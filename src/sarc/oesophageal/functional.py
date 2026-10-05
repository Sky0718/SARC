import csv
import math
from collections import defaultdict
from fractions import Fraction
import numpy as np
from .contracts import CONTRASTS, STATES, exact, require

CATEGORIES = {
    "core": ("BOTH_CORE", "ORGANOID_ONLY", "CELL_LINE_ONLY", "NOT_LISTED_IN_EITHER"),
    "role": ("ONCOGENE", "TSG", "DUAL_ONCOGENE_TSG", "OTHER_OR_UNSPECIFIED", "UNLISTED"),
}

def call_bounds(coefficients, calls):
    observed = Fraction()
    lower_unknown = Fraction()
    upper_unknown = Fraction()
    for (key, coefficient) in (coefficients.items()):
        state = calls[key]["state"]
        value = calls[key]["value"]
        if (state in ("MEASURED_POSITIVE", "MEASURED_NEGATIVE")):
            require(value in (0, 1), "Invalid measured classification")
            observed += coefficient * value
        else:
            require(value is None, "Unknown cell must not carry a numerical value")
            lower_unknown += min(Fraction(), coefficient)
            upper_unknown += max(Fraction(), coefficient)
    return {
        "observed": observed,
        "lower": observed + lower_unknown,
        "upper": observed + upper_unknown,
    }

def direction(lower, upper):
    if (upper < 0):
        return "NEGATIVE"
    if (lower > 0):
        return "POSITIVE"
    if (lower == upper == 0):
        return "ZERO"
    return "UNRESOLVED"

def average_rows(rows, fields):
    require(bool(rows), "No empty complete-cohort summaries")
    return {
        field: sum((row[field] for (row) in (rows)), Fraction()) / len(rows) for (field) in (fields)
    }

def evaluate_exact(allocation, scope, frozen, selected):
    calls = {(row["sample_id"], row["gene"]): row for (row) in (selected["calls"])}
    require(len(calls) == len(selected["calls"]), "Duplicate selected call")
    expected = {
        (model["sample_id"], gene)
        for (model) in (allocation["models"])
        for (state) in (STATES)
        for (gene) in (model["states"][state]["weights"])
    }
    require(set(calls) == expected, "Exact frozen selected-call support required")
    model_levels, model_contrasts, category_levels = [], [], []
    for (model) in (allocation["models"]):
        sample, donor = model["sample_id"], model["donor_id"]
        for (state) in (STATES):
            weights = {
                gene: exact(weight)
                for (gene, weight) in (model["states"][state]["weights"].items())
            }
            row = {
                "sample_id": sample,
                "donor_id": donor,
                "state": state,
                "R": sum(weights.values(), Fraction()),
                "positive": Fraction(),
                "zero": Fraction(),
                "absent": Fraction(),
                "null": Fraction(),
                "identity_conflict": Fraction(),
                "call_conflict": Fraction(),
            }
            key_map = {
                "MEASURED_POSITIVE": "positive",
                "MEASURED_NEGATIVE": "zero",
                "ABSENT": "absent",
                "EXPLICIT_NULL": "null",
                "IDENTITY_CONFLICT": "identity_conflict",
                "CONFLICTING_CALLS": "call_conflict",
            }
            for (gene, weight) in (weights.items()):
                row[key_map[calls[(sample, gene)]["state"]]] += weight
            row["unknown"] = (
                row["absent"] + row["null"] + row["identity_conflict"] + row["call_conflict"]
            )
            row["vacant"] = 10 - row["R"]
            row["lower"], row["upper"] = row["positive"], row["positive"] + row["unknown"]
            require(
                row["R"] == row["positive"] + row["zero"] + row["unknown"],
                "Absolute allocation accounting",
            )
            model_levels.append(row)
            for (partition, mapping) in (scope.get("categories", {}).items()):
                for (category) in (CATEGORIES[partition]):
                    mass = sum(
                        (
                            weight
                            for (gene, weight) in (weights.items())
                            if (mapping[gene] == category)
                        ),
                        Fraction(),
                    )
                    positive = sum(
                        (
                            weight
                            for (gene, weight) in (weights.items())
                            if (mapping[gene] == category
                            and calls[(sample, gene)]["state"] == "MEASURED_POSITIVE")
                        ),
                        Fraction(),
                    )
                    unknown = sum(
                        (
                            weight
                            for (gene, weight) in (weights.items())
                            if (mapping[gene] == category
                            and calls[(sample, gene)]["state"]
                            not in ("MEASURED_POSITIVE", "MEASURED_NEGATIVE"))
                        ),
                        Fraction(),
                    )
                    category_levels.append(
                        {
                            "sample_id": sample,
                            "donor_id": donor,
                            "state": state,
                            "partition": partition,
                            "category": category,
                            "R": mass,
                            "positive": positive,
                            "zero": mass - positive - unknown,
                            "unknown": unknown,
                            "pooled_positive_fraction": positive / mass if (mass) else None,
                        }
                    )
        for (contrast, signs) in (CONTRASTS.items()):
            coefficients = defaultdict(Fraction)
            for (state, sign) in (signs.items()):
                for (gene, weight) in (model["states"][state]["weights"].items()):
                    coefficients[(sample, gene)] += sign * exact(weight)
            bounds = call_bounds(coefficients, calls)
            model_contrasts.append(
                {
                    "sample_id": sample,
                    "donor_id": donor,
                    "contrast": contrast,
                    **bounds,
                    "direction": direction(bounds["lower"], bounds["upper"]),
                }
            )
    level_fields = [
        "R",
        "positive",
        "zero",
        "absent",
        "null",
        "identity_conflict",
        "call_conflict",
        "unknown",
        "vacant",
        "lower",
        "upper",
    ]
    donor_levels, donor_contrasts = [], []
    for (donor) in (frozen["donor_ids"]):
        for (state) in (STATES):
            rows = [
                row
                for (row) in (model_levels)
                if (row["donor_id"] == donor and row["state"] == state)
            ]
            donor_levels.append(
                {"donor_id": donor, "state": state, **average_rows(rows, level_fields)}
            )
        for (contrast) in (CONTRASTS):
            rows = [
                row
                for (row) in (model_contrasts)
                if (row["donor_id"] == donor and row["contrast"] == contrast)
            ]
            row = {
                "donor_id": donor,
                "contrast": contrast,
                **average_rows(rows, ["observed", "lower", "upper"]),
            }
            row["direction"] = direction(row["lower"], row["upper"])
            donor_contrasts.append(row)
    population_levels = [
        {
            "state": state,
            **average_rows(
                [row for (row) in (donor_levels) if (row["state"] == state)], level_fields
            ),
        }
        for (state) in (STATES)
    ]
    population_contrasts = []
    for (contrast) in (CONTRASTS):
        rows = [row for (row) in (donor_contrasts) if (row["contrast"] == contrast)]
        row = {
            "contrast": contrast,
            **average_rows(rows, ["observed", "lower", "upper"]),
            "donors": len(rows),
        }
        row["direction"] = direction(row["lower"], row["upper"])
        row["donor_directions"] = {
            sign: sum(value["direction"] == sign for (value) in (rows))
            for (sign) in (("NEGATIVE", "ZERO", "POSITIVE", "UNRESOLVED"))
        }
        population_contrasts.append(row)
    category_donors, category_population = [], []
    for (partition, mapping) in (scope.get("categories", {}).items()):
        for (category) in (CATEGORIES[partition]):
            for (state) in (STATES):
                for (donor) in (frozen["donor_ids"]):
                    rows = [
                        row
                        for (row) in (category_levels)
                        if (row["partition"] == partition
                        and row["category"] == category
                        and row["state"] == state
                        and row["donor_id"] == donor)
                    ]
                    average = average_rows(rows, ["R", "positive", "zero", "unknown"])
                    average["pooled_positive_fraction"] = (
                        average["positive"] / average["R"] if (average["R"]) else None
                    )
                    category_donors.append(
                        {
                            "donor_id": donor,
                            "state": state,
                            "partition": partition,
                            "category": category,
                            **average,
                        }
                    )
                rows = [
                    row
                    for (row) in (category_donors)
                    if (row["partition"] == partition
                    and row["category"] == category
                    and row["state"] == state)
                ]
                average = average_rows(rows, ["R", "positive", "zero", "unknown"])
                average["pooled_positive_fraction"] = (
                    average["positive"] / average["R"] if (average["R"]) else None
                )
                category_population.append(
                    {"state": state, "partition": partition, "category": category, **average}
                )
    unlinked = []
    for (unit) in (frozen["donor_ids"] + ["POPULATION"]):
        levels = {
            row["state"]: row
            for (row) in (population_levels if (unit == "POPULATION") else donor_levels)
            if (unit == "POPULATION" or row["donor_id"] == unit)
        }
        paired = next(
            row
            for (row) in (population_contrasts if (unit == "POPULATION") else donor_contrasts)
            if (row["contrast"] == "Delta" and (unit == "POPULATION" or row["donor_id"] == unit))
        )
        low = levels["G1H1"]["lower"] - levels["G0H0"]["upper"]
        high = levels["G1H1"]["upper"] - levels["G0H0"]["lower"]
        require(
            low <= paired["lower"] <= paired["upper"] <= high,
            "Paired shared-label interval must nest within unlinked arms",
        )
        unlinked.append(
            {
                "unit": unit,
                "unlinked_lower": low,
                "unlinked_upper": high,
                "shared_lower": paired["lower"],
                "shared_upper": paired["upper"],
                "width_reduction": (high - low) - (paired["upper"] - paired["lower"]),
                "unlinked_direction": direction(low, high),
                "shared_direction": paired["direction"],
            }
        )
    return {
        "model_levels": model_levels,
        "model_contrasts": model_contrasts,
        "donor_levels": donor_levels,
        "donor_contrasts": donor_contrasts,
        "population_levels": population_levels,
        "population_contrasts": population_contrasts,
        "model_category_levels": category_levels,
        "donor_category_levels": category_donors,
        "population_category_levels": category_population,
        "Delta_identity_gain": unlinked,
    }

def composition(donor_contrasts, donor_ids, output):
    require(
        donor_ids == sorted(donor_ids) and len(set(donor_ids)) == len(donor_ids),
        "Exact lexicographic donor order",
    )
    lookup = {(row["donor_id"], row["contrast"]): row for (row) in (donor_contrasts)}
    family = [("Delta", "lower"), ("Delta", "upper"), ("Theta", "lower"), ("Theta", "upper")]
    columns = [
        [lookup[(donor, contrast)][endpoint] for (donor) in (donor_ids)]
        for (contrast, endpoint) in (family)
    ]
    constant = [len(set(column)) == 1 for (column) in (columns)]
    n = len(donor_ids)
    require(n >= 2, "Composition diagnostic requires at least two donors")
    means_exact = [sum(column, Fraction()) / n for (column) in (columns)]
    variance_exact = [
        sum(((value - mean) ** 2 for (value) in (column)), Fraction()) / (n * (n - 1))
        for (column, mean) in (zip(columns, means_exact))
    ]
    se = np.array([math.sqrt(float(value)) for (value) in (variance_exact)])
    for (is_constant, value) in (zip(constant, se)):
        require(
            (value == 0 if (is_constant) else value > 0 and math.isfinite(value)),
            "Floating standard error invalid for exact variance class",
        )
    values = np.asarray([[float(columns[j][i]) for (j) in (range(4))] for (i) in (range(n))])
    means = np.asarray([float(value) for (value) in (means_exact)])
    generator = np.random.Generator(np.random.PCG64(20261005))
    draws = generator.integers(0, n, size = (10000, n), dtype = np.int64)
    np.save(output / "DONOR_DRAWS.npy", draws, allow_pickle = False)
    draw_means = values[draws].mean(axis = 1)
    t = np.zeros(10000)
    for (index) in (range(4)):
        if (not constant[index]):
            t = np.maximum(t, np.abs(draw_means[:, index] - means[index]) / se[index])
    ordered = np.sort(t)
    position = (10000 - 1) * 0.95
    lower = math.floor(position)
    q = float(ordered[lower] + (position - lower) * (ordered[lower + 1] - ordered[lower]))
    intervals = [
        {
            "endpoint": contrast + "_" + endpoint,
            "mean_exact": str(mean),
            "se_exact_squared": str(variance),
            "se": float(error),
            "constant": is_constant,
            "lower": float(number - q * error),
            "upper": float(number + q * error),
        }
        for ((contrast, endpoint), mean, variance, error, is_constant, number) in (zip(
            family, means_exact, variance_exact, se, constant, means
        ))
    ]
    np.save(output / "STUDENTISED_MAX.npy", t, allow_pickle = False)
    return {
        "family": [row["endpoint"] for (row) in (intervals)],
        "draws": 10000,
        "donor_ids": donor_ids,
        "seed": 20261005,
        "generator": "PCG64",
        "q95": q,
        "quantile_position_zero_based": position,
        "intervals": intervals,
        "outer_envelopes": {
            "Delta": [intervals[0]["lower"], intervals[1]["upper"]],
            "Theta": [intervals[2]["lower"], intervals[3]["upper"]],
        },
        "interpretation": "Conditional paired donor-composition diagnostic. This is not an independent validation, assay-error model, or cross-context test.",
    }

def serial(value):
    if (isinstance(value, Fraction)):
        return str(value)
    if (isinstance(value, dict)):
        return {key: serial(item) for (key, item) in (value.items())}
    if (isinstance(value, (list, tuple))):
        return [serial(item) for (item) in (value)]
    return value

def write_tables(output, results):
    for (name, rows) in (results.items()):
        if (not rows):
            continue
        with (output / (name + ".csv")).open("x", encoding = "utf-8", newline = "") as stream:
            writer = csv.DictWriter(stream, fieldnames = list(rows[0]), lineterminator = "\n")
            writer.writeheader()
            for (row) in (rows):
                writer.writerow(serial(row))
