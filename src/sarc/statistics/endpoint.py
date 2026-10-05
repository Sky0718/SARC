import math
from fractions import Fraction

import numpy as np

CONTRASTS = (
    "PRODIGY_joint",
    "DawnRank_joint",
    "PRODIGY_joint_minus_degree_update",
    "DawnRank_joint_minus_degree_update",
)

def require(condition, message):
    if (not condition):
        raise ValueError(message)

def as_fraction(value):
    require(
        (value is not None) and (not isinstance(value, bool)), "Missing or boolean endpoint"
    )
    if (isinstance(value, float)):
        require(math.isfinite(value), "Nonfinite endpoint")
    try:
        result = Fraction(value)
    except (ValueError, TypeError, ZeroDivisionError) as error:
        raise ValueError("Invalid rational endpoint") from error
    require(math.isfinite(float(result)), "Nonfinite converted endpoint")
    return result

def validate_scale(values, variance, standard_error):
    constant = all((value == values[0]) for (value) in (values))
    require(math.isfinite(standard_error), "Nonfinite standard error")
    if (constant):
        require((variance == 0) and (standard_error == 0), "Inconsistent constant coordinate")
    else:
        require((variance > 0) and (standard_error > 0), "Inconsistent zero standard error")
        require(
            (len(set(map(float, values))) == len(set(values))),
            "Distinct exact endpoints collapse in binary64",
        )
    return constant

def calculate(vectors, draw_values, names):
    require((len(vectors) >= 2), "At least two donors required")
    require((len(names) > 0), "Empty family")
    require(all((len(row) == len(names)) for (row) in (vectors)), "Endpoint width mismatch")
    exact = [[as_fraction(value) for (value) in (row)] for (row) in (vectors)]
    n = len(exact)
    require(all((len(row) == n) for (row) in (draw_values)), "Draw width mismatch")
    require((len(draw_values) >= 2), "At least two draws required")
    require(
        all(
            (isinstance(index, int) and (not isinstance(index, bool)) and (0 <= index < n))
            for (row) in (draw_values)
            for (index) in (row)
        ),
        "Invalid draw index",
    )
    draws = np.asarray(draw_values, dtype = np.int64)
    means = []
    errors = []
    constants = []
    variances = []
    resampled = np.empty((len(draws), len(names)), dtype = np.float64)
    for (column) in (range(len(names))):
        values = [row[column] for (row) in (exact)]
        mean = sum(values, Fraction()) / n
        variance = sum(((value - mean) ** 2 for (value) in (values)), Fraction()) / (n * (n - 1))
        standard_error = math.sqrt(float(variance))
        constant = validate_scale(values, variance, standard_error)
        denominator = math.lcm(*(value.denominator for (value) in (values)))
        integers = [(value.numerator * (denominator // value.denominator)) for (value) in (values)]
        require(
            (max(map(abs, integers)) * n <= np.iinfo(np.int64).max),
            "Integer accumulation exceeds int64",
        )
        sums = np.asarray(integers, dtype = np.int64)[draws].sum(axis = 1)
        if (constant):
            require(
                bool(np.all(sums == integers[0] * n)),
                "Nonzero exact deviation for constant coordinate",
            )
            resampled[:, column] = float(mean)
        else:
            resampled[:, column] = sums / (denominator * n)
        means.append(mean)
        errors.append(standard_error)
        variances.append(variance)
        constants.append(constant)
    centres = np.asarray(list(map(float, means)))
    active = [index for (index, constant) in (enumerate(constants)) if (not constant)]
    deviations = np.zeros_like(resampled)
    if (active):
        deviations[:, active] = (
            np.abs(resampled[:, active] - centres[active]) / np.asarray(errors)[active]
        )
    require(bool(np.all(np.isfinite(deviations))), "Nonfinite standardised deviation")
    maximum = deviations.max(axis = 1)
    critical = float(np.quantile(maximum, 0.95, method = "linear"))
    endpoints = []
    for (index, name) in (enumerate(names)):
        mean = float(means[index])
        standard_error = errors[index]
        coordinate_critical = float(np.quantile(deviations[:, index], 0.95, method = "linear"))
        endpoints.append(
            {
                "name": name,
                "mean_exact": str(means[index]),
                "mean": mean,
                "variance_of_mean_exact": str(variances[index]),
                "original_standard_error": standard_error,
                "constant": constants[index],
                "simultaneous_interval": [
                    mean - critical * standard_error,
                    mean + critical * standard_error,
                ],
                "unadjusted_critical_value": coordinate_critical,
                "unadjusted_interval": [
                    mean - coordinate_critical * standard_error,
                    mean + coordinate_critical * standard_error,
                ],
            }
        )
    return {
        "n_donors": n,
        "replicates": len(draws),
        "family_size": len(names),
        "endpoint_order": list(names),
        "critical_value": critical,
        "constant_coordinate_count": sum(constants),
        "endpoints": endpoints,
        "quantile_rule": "Linear interpolation at zero-based (B - 1) * 0.95",
        "fixed_original_standard_errors": True,
    }, {"means": resampled.tolist(), "maximum_standardised_deviation": maximum.tolist()}

def colorectal_family(paired, draws, donor_ids):
    require((len(draws) == 10000), "All 10000 saved draws required")
    require(
        (len(donor_ids) == len(set(donor_ids)) == 83),
        "Full distinct 83-donor population required",
    )
    require((donor_ids == sorted(donor_ids)), "Lexicographic donor order required")
    vectors = [{"donor_id": donor, "values": []} for (donor) in (donor_ids)]
    expected_members = None
    for (contrast) in (CONTRASTS):
        source = paired["update_comparator"][contrast]
        members = {
            donor: sorted(
                row["model_id"] for (row) in (source["models"]) if (row["donor_id"] == donor)
            )
            for (donor) in (donor_ids)
        }
        model_ids = [member for (group) in (members.values()) for (member) in (group)]
        require(
            (len(model_ids) == len(set(model_ids)) == 85), "Full 85-model population required"
        )
        require(
            (expected_members is None or members == expected_members),
            "Donor membership mismatch",
        )
        expected_members = members
        indexed = {row["donor_id"]: row for (row) in (source["donors"])}
        require(
            (len(indexed) == len(source["donors"]) == 83 and set(indexed) == set(donor_ids)),
            "Contrast donor mismatch",
        )
        for (vector) in (vectors):
            donor = vector["donor_id"]
            row = indexed[donor]
            require((row["models"] == len(members[donor])), "Within-donor count mismatch")
            require(
                (as_fraction(row["lower"]) <= as_fraction(row["upper"])), "Reversed endpoints"
            )
            vector["model_ids"] = members[donor]
            vector["values"].extend([row["lower"], row["upper"]])
        for (bound) in (("lower", "upper")):
            mean = sum((as_fraction(row[bound]) for (row) in (indexed.values())), Fraction()) / 83
            require(
                (mean == as_fraction(source["population"][bound])),
                "Population endpoint mismatch",
            )
    names = [
        name + "_" + bound
        for (name) in (("PRODIGY_joint", "DawnRank_joint", "PRODIGY_Theta", "DawnRank_Theta"))
        for (bound) in (("lower", "upper"))
    ]
    summary, bootstrap = calculate([row["values"] for (row) in (vectors)], draws, names)
    summary["models"] = 85
    summary["budget"] = 10
    summary["donor_ids"] = donor_ids
    summary["envelopes"] = {}
    for (index, contrast) in (enumerate(CONTRASTS)):
        lower, upper = summary["endpoints"][2 * index : 2 * index + 2]
        summary["envelopes"][contrast] = {
            "identification_exact": [lower["mean_exact"], upper["mean_exact"]],
            "simultaneous_outer": [
                lower["simultaneous_interval"][0],
                upper["simultaneous_interval"][1],
            ],
            "unadjusted_outer": [
                lower["unadjusted_interval"][0],
                upper["unadjusted_interval"][1],
            ],
        }
    return summary, vectors, bootstrap
