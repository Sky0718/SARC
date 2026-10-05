import math
from fractions import Fraction
from .matched_return import quantile, require

def uncertainty(vectors, indices, production = True):
    require(
        len(vectors) == 6
        and len({len(values) for (values) in (vectors.values())}) == 1,
        "FIXED_SIX_FAMILY_REQUIRED",
    )
    n = len(next(iter(vectors.values())))
    require(
        n >= 2
        and indices
        and all(
            (
                len(row) == n
                and all((type(index) is int and 0 <= index < n for (index) in (row)))
                for (row) in (indices)
            )
        ),
        "PAIRED_DRAW_GEOMETRY_INVALID",
    )
    if (production):
        require(n == 58 and len(indices) == 2000, "ORIGINAL_HELDOUT_DRAWS_REQUIRED")
    records, joint = ([], [])
    for (label, values) in (vectors.items()):
        if (any((value is None for (value) in (values)))):
            records.append(
                {
                    "contrast": label,
                    "status": "FULL_DONOR_CONTRAST_UNAVAILABLE",
                    "point": None,
                    "percentile_95": None,
                    "standard_error": None,
                }
            )
            continue
        require(
            all(
                (
                    isinstance(value, Fraction) or type(value) is int
                    for (value) in (values)
                )
            ),
            "EXACT_RATIONAL_DONOR_VALUES_REQUIRED",
        )
        exact_values = [Fraction(value) for (value) in (values)]
        exact_mean = sum(exact_values, Fraction(0)) / n
        exact_variance = sum(
            ((value - exact_mean) ** 2 for (value) in (exact_values)), Fraction(0)
        ) / (n * (n - 1))
        if (exact_variance == 0):
            point = float(exact_mean)
            draws = [point] * len(indices)
            scale = 0.0
        else:
            values = [float(value) for (value) in (exact_values)]
            point = sum(values) / n
            draws = [
                sum((values[index] for (index) in (row))) / n for (row) in (indices)
            ]
            scale = math.sqrt(float(exact_variance))
        records.append(
            {
                "contrast": label,
                "status": "AVAILABLE",
                "point": point,
                "percentile_95": [quantile(draws, 0.025), quantile(draws, 0.975)],
                "standard_error": scale,
            }
        )
        if (exact_variance > 0 and scale > 0 and math.isfinite(scale)):
            joint.append([abs(value - point) / scale for (value) in (draws)])
    available = len(joint) == 6
    critical = (
        quantile([max(row) for (row) in (zip(*joint))], 0.95) if (available) else None
    )
    for (row) in (records):
        row["simultaneous_sensitivity_95"] = (
            [
                row["point"] - critical * row["standard_error"],
                row["point"] + critical * row["standard_error"],
            ]
            if (available)
            else None
        )
    return {
        "status": "AVAILABLE_FIXED_SIX_FAMILY"
        if (available)
        else "JOINT_FAMILY_UNAVAILABLE_ZERO_SCALE_OR_MISSING",
        "draw_count": len(indices),
        "donor_count": n,
        "critical": critical,
        "interval_kind": "approximate_fixed_original_SE_maximum_deviation_sensitivity_not_exact_studentized_inference",
        "records": records,
    }
