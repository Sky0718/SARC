from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations
from .allocation import (
    BUDGETS,
    CELLS,
    CONTRASTS,
    METHODS,
    SEEDS,
    bound,
    distance,
    mean_weights,
    require,
    sign,
)

def reduce_values(rows, weights, fields):
    require(len(rows) == len(weights), "Reduction length mismatch")
    return {
        field: sum((row[field] * weight for (row, weight) in (zip(rows, weights))), Fraction())
        for (field) in (fields)
    }

def arm_mass(weights, sample, calls, budget):
    result = {
        key: Fraction()
        for (key) in (("returned", "known_positive", "known_zero", "explicit_null", "absent"))
    }
    for (gene, value) in (weights.items()):
        state, outcome = calls[(sample, gene)]
        field = (
            "known_positive"
            if (state == "MEASURED" and outcome == 1)
            else "known_zero"
            if (state == "MEASURED")
            else "explicit_null"
            if (state == "EXPLICIT_NULL")
            else "absent"
        )
        result[field] += value
        result["returned"] += value
    result["unknown"] = result["explicit_null"] + result["absent"]
    result["vacant"] = budget - result["returned"]
    require(
        result["returned"]
        == result["known_positive"] + result["known_zero"] + result["unknown"]
        and result["vacant"] >= 0,
        "Budget partition failure",
    )
    return result

def contrast(arms, vector, sample, calls):
    genes = set().union(*(set(arms[cell]) for (cell) in (vector)))
    coefficients = {
        (sample, gene): sum(
            (factor * arms[cell].get(gene, Fraction()) for (cell, factor) in (vector.items())),
            Fraction(),
        )
        for (gene) in (genes)
    }
    result = bound(coefficients, calls)
    independent_lower, independent_upper = Fraction(), Fraction()
    for (cell, factor) in (vector.items()):
        arm = bound({(sample, gene): value for (gene, value) in (arms[cell].items())}, calls)
        independent_lower += factor * arm["lower" if (factor > 0) else "upper"]
        independent_upper += factor * arm["upper" if (factor > 0) else "lower"]
    require(
        independent_lower <= result["lower"] <= result["upper"] <= independent_upper,
        "Shared bound outside marginal-arm envelope",
    )
    result.update(
        independent_lower = independent_lower,
        independent_upper = independent_upper,
        independent_width = independent_upper - independent_lower,
        independent_sign = sign(independent_lower, independent_upper),
        width_saved = independent_upper - independent_lower - result["width"],
    )
    return result

def direction_summary(rows):
    return dict(Counter(row["sign"] for (row) in (rows)))

def percentile(values, proportion):
    if (not values):
        return None
    ordered = sorted(values)
    index = (len(ordered) - 1) * proportion
    floor = index.numerator // index.denominator
    remainder = index - floor
    return (
        ordered[floor]
        if (floor == len(ordered) - 1)
        else ordered[floor] * (1 - remainder) + ordered[floor + 1] * remainder
    )

def distribution(values):
    valid = [value for (value) in (values) if (value is not None)]
    return {
        "total": len(values),
        "defined": len(valid),
        "undefined": len(values) - len(valid),
        "minimum": min(valid) if (valid) else None,
        "q25": percentile(valid, Fraction(1, 4)),
        "median": percentile(valid, Fraction(1, 2)),
        "q75": percentile(valid, Fraction(3, 4)),
        "maximum": max(valid) if (valid) else None,
    }

def overlap(first, second):
    genes = set(first) | set(second)
    intersection = sum(
        (min(first.get(gene, Fraction()), second.get(gene, Fraction())) for (gene) in (genes)),
        Fraction(),
    )
    union = sum(
        (max(first.get(gene, Fraction()), second.get(gene, Fraction())) for (gene) in (genes)),
        Fraction(),
    )
    return {
        "intersection": intersection,
        "union": union,
        "ratio": intersection / union if (union) else None,
    }

def action(arms, direction):
    values = {}
    for (gene) in (set(arms["G0H0"]) | set(arms["G1H1"])):
        delta = arms["G1H1"].get(gene, Fraction()) - arms["G0H0"].get(gene, Fraction())
        values[gene] = max(delta if (direction == "addition") else -delta, 0)
    return values

def run(models, allocations, calls):
    donors = defaultdict(list)
    for (model) in (models):
        donors[model["donor_id"]].append(model)
    require(
        len(models) == 85 and len(donors) == 83 and len(allocations) == 8160,
        "Full analysis population differs",
    )
    alphas = {row["model_id"]: Fraction(1, 83 * len(donors[row["donor_id"]])) for (row) in (models)}
    means = {}
    for (method) in (METHODS):
        for (model) in (models):
            for (budget) in (BUDGETS):
                means[(method, model["model_id"], budget)] = {
                    cell: mean_weights(
                        [
                            allocations[(method, model["model_id"], cell, seed, budget)]
                            for (seed) in (SEEDS if (method == "PRODIGY") else ("deterministic",))
                        ]
                    )
                    for (cell) in (CELLS)
                }
    levels, contrasts, seeds_out, distances, variations, cancellations, overlaps = (
        [],
        [],
        [],
        [],
        [],
        [],
        [],
    )
    contrast_fields = (
        "known",
        "lower",
        "upper",
        "width",
        "independent_lower",
        "independent_upper",
        "independent_width",
        "width_saved",
    )
    summaries = {
        "scope": {"models": 85, "donors": 83, "budgets": list(BUDGETS), "seeds": list(SEEDS)},
        "budget_contrasts": [],
        "evidence_allocation": [],
        "seed_joint": [],
        "seed_variability": [],
        "turnover_cancellation": [],
        "overlap": [],
        "paired_directions": [],
    }
    for (method) in (METHODS):
        for (budget) in (BUDGETS):
            per_model_contrasts, per_model_levels = {}, {}
            for (model) in (models):
                mid, sample = model["model_id"], model["sample_id"]
                arms = means[(method, mid, budget)]
                identity = {
                    "method": method,
                    "budget": budget,
                    "unit": "model",
                    "model_id": mid,
                    "sample_id": sample,
                    "donor_id": model["donor_id"],
                }
                for (cell) in (CELLS):
                    value = arm_mass(arms[cell], sample, calls, budget)
                    per_model_levels[(mid, cell)] = value
                    levels.append({**identity, "cell": cell, **value})
                for (name, vector) in (CONTRASTS.items()):
                    value = contrast(arms, vector, sample, calls)
                    per_model_contrasts[(mid, name)] = value
                    contrasts.append({**identity, "contrast": name, **value})
                for (name) in (("evidence_at_old_graph", "evidence_at_new_graph", "joint")):
                    vector = CONTRASTS[name]
                    old = next(cell for (cell, factor) in (vector.items()) if (factor == -1))
                    new = next(cell for (cell, factor) in (vector.items()) if (factor == 1))
                    value = distance(arms[old], arms[new])
                    distances.append(
                        {
                            **identity,
                            "seed": "mean",
                            "contrast": name,
                            **value,
                            "known_contrast": per_model_contrasts[(mid, name)]["known"],
                            "bound_sign": per_model_contrasts[(mid, name)]["sign"],
                        }
                    )
                    for (seed) in (SEEDS if (method == "PRODIGY") else ("deterministic",)):
                        first, second = (
                            allocations[(method, mid, cell, seed, budget)]
                            for (cell) in ((old, new))
                        )
                        distances.append(
                            {
                                **identity,
                                "seed": seed,
                                "contrast": name,
                                **distance(first, second),
                            }
                        )
            donor_contrasts = {}
            for (donor, group) in (sorted(donors.items())):
                weights = [Fraction(1, len(group))] * len(group)
                identity = {
                    "method": method,
                    "budget": budget,
                    "unit": "donor",
                    "model_id": "",
                    "sample_id": "",
                    "donor_id": donor,
                }
                for (cell) in (CELLS):
                    rows = [per_model_levels[(row["model_id"], cell)] for (row) in (group)]
                    levels.append(
                        {
                            **identity,
                            "cell": cell,
                            **reduce_values(rows, weights, rows[0].keys()),
                        }
                    )
                for (name) in (CONTRASTS):
                    rows = [per_model_contrasts[(row["model_id"], name)] for (row) in (group)]
                    value = reduce_values(rows, weights, contrast_fields)
                    value.update(
                        sign = sign(value["lower"], value["upper"]),
                        independent_sign = sign(
                            value["independent_lower"], value["independent_upper"]
                        ),
                    )
                    donor_contrasts[(donor, name)] = value
                    contrasts.append({**identity, "contrast": name, **value})
            for (cell) in (CELLS):
                rows = [per_model_levels[(model["model_id"], cell)] for (model) in (models)]
                value = reduce_values(
                    rows, [alphas[model["model_id"]] for (model) in (models)], rows[0].keys()
                )
                levels.append(
                    {
                        "method": method,
                        "budget": budget,
                        "unit": "population",
                        "model_id": "",
                        "sample_id": "",
                        "donor_id": "",
                        "cell": cell,
                        **value,
                    }
                )
            for (name) in (CONTRASTS):
                rows = [donor_contrasts[(donor, name)] for (donor) in (sorted(donors))]
                value = reduce_values(rows, [Fraction(1, 83)] * 83, contrast_fields)
                value.update(
                    sign = sign(value["lower"], value["upper"]),
                    independent_sign = sign(
                        value["independent_lower"], value["independent_upper"]
                    ),
                )
                contrasts.append(
                    {
                        "method": method,
                        "budget": budget,
                        "unit": "population",
                        "model_id": "",
                        "sample_id": "",
                        "donor_id": "",
                        "contrast": name,
                        **value,
                    }
                )
                summaries["budget_contrasts"].append(
                    {
                        "method": method,
                        "budget": budget,
                        "contrast": name,
                        **value,
                        "donor_signs": direction_summary(rows),
                        "donors_shared_identified_independent_unresolved": sum(
                            row["sign"] != "unresolved"
                            and row["independent_sign"] == "unresolved"
                            for (row) in (rows)
                        ),
                        "donors_width_strictly_reduced": sum(
                            row["width_saved"] > 0 for (row) in (rows)
                        ),
                    }
                )
            for (name) in (("evidence_at_old_graph", "evidence_at_new_graph")):
                rows = [
                    row
                    for (row) in (distances)
                    if (row["method"] == method
                    and row["budget"] == budget
                    and row["contrast"] == name
                    and row["seed"] == "mean")
                ]
                summaries["evidence_allocation"].append(
                    {
                        "method": method,
                        "budget": budget,
                        "contrast": name,
                        "models": len(rows),
                        "allocation_changed_models": sum(row["l1"] > 0 for (row) in (rows)),
                        "allocation_changed_with_zero_known_contrast": sum(
                            row["l1"] > 0 and row["known_contrast"] == 0 for (row) in (rows)
                        ),
                        "weighted_l1": sum(
                            (alphas[row["model_id"]] * row["l1"] for (row) in (rows)), Fraction()
                        ),
                        "l1_distribution": distribution([row["l1"] for (row) in (rows)]),
                        "functional_signs": dict(Counter(row["bound_sign"] for (row) in (rows))),
                    }
                )
    for (budget) in (BUDGETS):
        for (seed) in (SEEDS):
            donor_values = defaultdict(list)
            for (model) in (models):
                arms = {
                    cell: allocations[("PRODIGY", model["model_id"], cell, seed, budget)]
                    for (cell) in (CELLS)
                }
                value = contrast(arms, CONTRASTS["joint"], model["sample_id"], calls)
                donor_values[model["donor_id"]].append(value)
                seeds_out.append(
                    {
                        "budget": budget,
                        "seed": seed,
                        "unit": "model",
                        "model_id": model["model_id"],
                        "donor_id": model["donor_id"],
                        **value,
                    }
                )
            reduced = {}
            for (donor, rows) in (sorted(donor_values.items())):
                value = reduce_values(
                    rows, [Fraction(1, len(rows))] * len(rows), contrast_fields
                )
                value.update(
                    sign = sign(value["lower"], value["upper"]),
                    independent_sign = sign(
                        value["independent_lower"], value["independent_upper"]
                    ),
                )
                reduced[donor] = value
                seeds_out.append(
                    {
                        "budget": budget,
                        "seed": seed,
                        "unit": "donor",
                        "model_id": "",
                        "donor_id": donor,
                        **value,
                    }
                )
            value = reduce_values(
                list(reduced.values()), [Fraction(1, 83)] * 83, contrast_fields
            )
            value.update(
                sign = sign(value["lower"], value["upper"]),
                independent_sign = sign(value["independent_lower"], value["independent_upper"]),
            )
            dawn = {
                row["donor_id"]: row
                for (row) in (contrasts)
                if (row["method"] == "DawnRank"
                and row["budget"] == budget
                and row["contrast"] == "joint"
                and row["unit"] == "donor")
            }
            summaries["seed_joint"].append(
                {
                    "budget": budget,
                    "seed": seed,
                    **value,
                    "donor_signs": direction_summary(list(reduced.values())),
                    "opposing_DawnRank_donors": sum(
                        {row["sign"], dawn[donor]["sign"]} == {"positive", "negative"}
                        for (donor, row) in (reduced.items())
                    ),
                }
            )
            seeds_out.append(
                {
                    "budget": budget,
                    "seed": seed,
                    "unit": "population",
                    "model_id": "",
                    "donor_id": "",
                    **value,
                }
            )
        for (model) in (models):
            mid = model["model_id"]
            for (cell) in (CELLS):
                for (first, second) in (combinations(SEEDS, 2)):
                    value = distance(
                        allocations[("PRODIGY", mid, cell, first, budget)],
                        allocations[("PRODIGY", mid, cell, second, budget)],
                    )
                    variations.append(
                        {
                            "budget": budget,
                            "model_id": mid,
                            "donor_id": model["donor_id"],
                            "cell": cell,
                            "seed_a": first,
                            "seed_b": second,
                            **value,
                        }
                    )
            paired = [
                distance(
                    allocations[("PRODIGY", mid, "G0H0", seed, budget)],
                    allocations[("PRODIGY", mid, "G1H1", seed, budget)],
                )
                for (seed) in (SEEDS)
            ]
            mean = reduce_values(paired, [Fraction(1, 5)] * 5, ("l1", "addition", "removal"))
            net = distance(
                means[("PRODIGY", mid, budget)]["G0H0"], means[("PRODIGY", mid, budget)]["G1H1"]
            )
            require(
                mean["l1"] >= net["l1"]
                and mean["addition"] >= net["addition"]
                and mean["removal"] >= net["removal"],
                "Jensen cancellation identity violated",
            )
            cancellations.append(
                {
                    "budget": budget,
                    "model_id": mid,
                    "donor_id": model["donor_id"],
                    "mean_seed_l1": mean["l1"],
                    "net_mean_l1": net["l1"],
                    "cancelled_l1": mean["l1"] - net["l1"],
                    "cancelled_fraction": (mean["l1"] - net["l1"]) / mean["l1"]
                    if (mean["l1"])
                    else None,
                    "mean_seed_addition": mean["addition"],
                    "net_mean_addition": net["addition"],
                    "mean_seed_removal": mean["removal"],
                    "net_mean_removal": net["removal"],
                }
            )
        for (cell) in (CELLS):
            rows = [
                row for (row) in (variations) if (row["budget"] == budget and row["cell"] == cell)
            ]
            summaries["seed_variability"].append(
                {
                    "budget": budget,
                    "cell": cell,
                    "model_seed_pairs": len(rows),
                    "models_with_seed_allocation_change": len(
                        {row["model_id"] for (row) in (rows) if (row["l1"] > 0)}
                    ),
                    "changed_pairs": sum(row["l1"] > 0 for (row) in (rows)),
                    "donor_weighted_mean_pair_l1": sum(
                        (alphas[row["model_id"]] * row["l1"] / 10 for (row) in (rows)), Fraction()
                    ),
                    "pair_l1_distribution": distribution([row["l1"] for (row) in (rows)]),
                }
            )
        rows = [row for (row) in (cancellations) if (row["budget"] == budget)]
        summary = reduce_values(
            rows,
            [alphas[row["model_id"]] for (row) in (rows)],
            (
                "mean_seed_l1",
                "net_mean_l1",
                "cancelled_l1",
                "mean_seed_addition",
                "net_mean_addition",
                "mean_seed_removal",
                "net_mean_removal",
            ),
        )
        summary.update(
            budget = budget, models_with_cancellation = sum(row["cancelled_l1"] > 0 for (row) in (rows))
        )
        summaries["turnover_cancellation"].append(summary)
        for (direction) in (("addition", "removal")):
            model_overlap = {}
            for (model) in (models):
                value = overlap(
                    action(means[("PRODIGY", model["model_id"], budget)], direction),
                    action(means[("DawnRank", model["model_id"], budget)], direction),
                )
                model_overlap[model["model_id"]] = value
                overlaps.append(
                    {
                        "budget": budget,
                        "direction": direction,
                        "unit": "model",
                        "model_id": model["model_id"],
                        "donor_id": model["donor_id"],
                        **value,
                    }
                )
            donor_overlap = []
            for (donor, group) in (sorted(donors.items())):
                value = reduce_values(
                    [model_overlap[row["model_id"]] for (row) in (group)],
                    [Fraction(1, len(group))] * len(group),
                    ("intersection", "union"),
                )
                value["ratio"] = (
                    value["intersection"] / value["union"] if (value["union"]) else None
                )
                donor_overlap.append(value)
                overlaps.append(
                    {
                        "budget": budget,
                        "direction": direction,
                        "unit": "donor",
                        "model_id": "",
                        "donor_id": donor,
                        **value,
                    }
                )
            value = reduce_values(
                donor_overlap, [Fraction(1, 83)] * 83, ("intersection", "union")
            )
            value["ratio"] = value["intersection"] / value["union"] if (value["union"]) else None
            overlaps.append(
                {
                    "budget": budget,
                    "direction": direction,
                    "unit": "population",
                    "model_id": "",
                    "donor_id": "",
                    **value,
                }
            )
            summaries["overlap"].append(
                {
                    "budget": budget,
                    "direction": direction,
                    "population": value,
                    "donor_ratio_distribution": distribution(
                        [row["ratio"] for (row) in (donor_overlap)]
                    ),
                    "zero_shared_positive_union_donors": sum(
                        row["intersection"] == 0 and row["union"] > 0 for (row) in (donor_overlap)
                    ),
                }
            )
        pro = {
            row["donor_id"]: row
            for (row) in (contrasts)
            if (row["method"] == "PRODIGY"
            and row["budget"] == budget
            and row["contrast"] == "joint"
            and row["unit"] == "donor")
        }
        dawn = {
            row["donor_id"]: row
            for (row) in (contrasts)
            if (row["method"] == "DawnRank"
            and row["budget"] == budget
            and row["contrast"] == "joint"
            and row["unit"] == "donor")
        }
        summaries["paired_directions"].append(
            {
                "budget": budget,
                "complete_donors": len(pro),
                "matrix": dict(
                    Counter(
                        pro[donor]["sign"] + "_" + dawn[donor]["sign"]
                        for (donor) in (sorted(donors))
                    )
                ),
                "opposing": sum(
                    {pro[donor]["sign"], dawn[donor]["sign"]} == {"positive", "negative"}
                    for (donor) in (sorted(donors))
                ),
            }
        )
    return {
        "summary": summaries,
        "arm_levels": levels,
        "contrasts": contrasts,
        "seed_joint": seeds_out,
        "allocation_distances": distances,
        "seed_variability": variations,
        "turnover_cancellation": cancellations,
        "donor_overlap": overlaps,
    }
