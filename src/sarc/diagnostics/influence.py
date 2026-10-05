import math
from collections import Counter, defaultdict
from fractions import Fraction

CELLS = ("G0H0", "G1H0", "G0H1", "G1H1")

SEEDS = (104729, 130363, 155921, 196613, 228017)

MASS = (
    "returned",
    "known_positive",
    "known_zero",
    "explicit_null",
    "absent",
    "unknown",
    "vacant",
)

def need(condition, message):
    if (not condition):
        raise ValueError(message)

def labels_from(carriers):
    calls, origins = {}, defaultdict(list)
    identities = set()
    for (name, carrier) in (carriers):
        binding = carrier["source_binding"]
        need(
            isinstance(binding.get("source_id"), str) and binding["source_id"],
            "Explicit label source ID required",
        )
        identities.add(binding["source_id"])
        local = set()
        for (row) in (carrier["calls"]):
            key = (row["sample_id"], row["gene"])
            state, value = row["state"], row["value"]
            need(key not in local, "Duplicate label within carrier")
            local.add(key)
            need(state in ("MEASURED", "EXPLICIT_NULL", "ABSENT"), "Unknown assay state")
            need(
                (state == "MEASURED" and type(value) is int and value in (0, 1))
                or (state != "MEASURED" and value is None),
                "Invalid assay value/state",
            )
            need(
                key not in calls or calls[key] == (state, value),
                "Conflicting shared Table5 labels",
            )
            calls[key] = (state, value)
            origins[key].append(name)
    need(len(identities) == 1, "Selected calls do not refer to the same original Table5 source")
    return calls, origins

def sign(lower, upper):
    return (
        "STRICTLY_NEGATIVE"
        if (upper < 0)
        else "STRICTLY_POSITIVE"
        if (lower > 0)
        else "IDENTIFIED_ZERO"
        if (lower == upper == 0)
        else "UNRESOLVED"
    )

def bounds(coefficients, calls):
    known, negative, positive = Fraction(), Fraction(), Fraction()
    for (key, coefficient) in (coefficients.items()):
        state, value = calls[key]
        if (state == "MEASURED"):
            known += coefficient * value
        else:
            negative += min(coefficient, 0)
            positive += max(coefficient, 0)
    lower, upper = known + negative, known + positive
    return {
        "known": known,
        "unknown_negative": negative,
        "unknown_positive": positive,
        "lower": lower,
        "upper": upper,
        "width": upper - lower,
        "sign": sign(lower, upper),
    }

def mass(weights, sample, calls):
    result = dict.fromkeys(MASS, Fraction())
    counts = Counter()
    for (gene, weight) in (weights.items()):
        need(weight > 0, "Only positive allocations enter absolute mass")
        state, value = calls[(sample, gene)]
        group = (
            "known_positive"
            if (state == "MEASURED" and value == 1)
            else "known_zero"
            if (state == "MEASURED")
            else "explicit_null"
            if (state == "EXPLICIT_NULL")
            else "absent"
        )
        result["returned"] += weight
        result[group] += weight
        counts[group] += 1
    result["unknown"] = result["explicit_null"] + result["absent"]
    result["vacant"] = 10 - result["returned"]
    need(
        result["returned"]
        == result["known_positive"] + result["known_zero"] + result["unknown"],
        "Absolute mass partition failure",
    )
    need(
        result["vacant"] >= 0 and result["returned"] + result["vacant"] == 10,
        "Capacity partition failure",
    )
    return result, counts

def weighted(records):
    return {
        field: sum((weight * row[field] for (weight, row) in (records)), Fraction())
        for (field) in (MASS)
    }

def witness(rows, indices, gain):
    selected = [rows[index] for (index) in (indices)]
    need(
        sum((row["upward_gain"] for (row) in (selected)), Fraction()) == gain,
        "Flip witness gain differs",
    )
    allocated_mass = (
        sum((row["joint_arm_allocated_mass"] for (row) in (selected)), Fraction())
        if (all("joint_arm_allocated_mass" in row for (row) in (selected)))
        else None
    )
    return {
        "flip_count": len(selected),
        "contrast_leverage_burden": gain,
        "joint_arm_allocated_mass_of_witness": allocated_mass,
        "distinct_genes": len({row["gene"] for (row) in (selected)}),
        "distinct_donors": len({row["donor_id"] for (row) in (selected)}),
        "sample_gene_ids": [[row["sample_id"], row["gene"]] for (row) in (selected)],
    }

def minimum_burden(rows, threshold, strict):
    denominator = math.lcm(
        threshold.denominator, *(row["upward_gain"].denominator for (row) in (rows))
    )
    target = int(threshold * denominator) + int(strict)
    values = [int(row["upward_gain"] * denominator) for (row) in (rows)]
    total = sum(values)
    if (target <= 0):
        return {"status": "EXACT", **witness(rows, [], Fraction())}
    if (target > total):
        return {"status": "UNREACHABLE"}
    if (total > 2000000 or len(values) * (total + 1) > 300000000):
        return {
            "status": "EXACT_STATE_SPACE_LIMIT",
            "scaled_total": total,
            "denominator": denominator,
        }
    history = [1]
    for (value) in (values):
        history.append(history[-1] | (history[-1] << value))
    shifted = history[-1] >> target
    if (not shifted):
        return {"status": "UNREACHABLE"}
    lowbit = shifted & -shifted
    selected_total = target + lowbit.bit_length() - 1
    remaining, selected = selected_total, []
    for (index) in (reversed(range(len(values)))):
        if ((history[index] >> remaining) & 1):
            continue
        selected.append(index)
        remaining -= values[index]
        need(remaining >= 0, "Subset reconstruction underflow")
    need(remaining == 0, "Subset reconstruction failed")
    return {
        "status": "EXACT",
        "integer_denominator": denominator,
        **witness(rows, sorted(selected), Fraction(selected_total, denominator)),
    }

def flip_analysis(coefficients, calls, sample_info, allocation_mass = None):
    base = bounds(coefficients, calls)
    need(
        base["upper"] < 0,
        "Primary flip diagnostic requires the accepted strictly negative joint contrast",
    )
    cells, beneficial = [], []
    for ((sample, gene), coefficient) in (sorted(coefficients.items())):
        state, value = calls[(sample, gene)]
        if (state != "MEASURED" or not coefficient):
            continue
        delta = coefficient * (1 - 2 * value)
        row = {
            "sample_id": sample,
            "gene": gene,
            "donor_id": sample_info[sample]["donor_id"],
            "label": value,
            "coefficient": coefficient,
            "flip_delta": delta,
            "upward_gain": max(delta, 0),
        }
        if (allocation_mass is not None):
            row["joint_arm_allocated_mass"] = allocation_mass[(sample, gene)]
        cells.append(row)
        if (delta > 0):
            beneficial.append(row)
    beneficial.sort(key = lambda row: (-row["upward_gain"], row["sample_id"], row["gene"]))
    requirements = {
        "lose_strict_negative": (-base["upper"], False),
        "positive_completion_possible": (-base["upper"], True),
        "strictly_positive_identified": (-base["lower"], True),
    }
    solutions, curves, gain = {}, [], Fraction()
    for (count) in (range(len(beneficial) + 1)):
        if (count):
            gain += beneficial[count - 1]["upward_gain"]
        curves.append(
            {
                "flip_count": count,
                "maximum_upward_gain": gain,
                "lower_after": base["lower"] + gain,
                "upper_after": base["upper"] + gain,
                "sign_after": sign(base["lower"] + gain, base["upper"] + gain),
            }
        )
        for (name, (threshold, strict)) in (requirements.items()):
            if (name not in solutions and (gain > threshold if (strict) else gain >= threshold)):
                solutions[name] = {
                    "minimum_count": {
                        "status": "EXACT",
                        **witness(beneficial, list(range(count)), gain),
                    }
                }
    for (name, (threshold, strict)) in (requirements.items()):
        solutions.setdefault(name, {"minimum_count": {"status": "UNREACHABLE"}})
        solutions[name]["threshold_gain"] = threshold
        solutions[name]["strict_threshold"] = strict
        solutions[name]["minimum_burden"] = minimum_burden(beneficial, threshold, strict)
    return (
        {
            "base": base,
            "measured_nonzero_coefficient_calls": len(cells),
            "beneficial_upward_flips": len(beneficial),
            "solutions": solutions,
        },
        cells,
        curves,
    )

def gene_analysis(coefficients, calls, genes, sample_info):
    base, records = bounds(coefficients, calls), []
    for (gene) in (sorted(genes)):
        subset = {key: value for (key, value) in (coefficients.items()) if (key[1] == gene)}
        value = bounds(subset, calls)
        remaining = {key: value for (key, value) in (coefficients.items()) if (key[1] != gene)}
        omitted = bounds(remaining, calls)
        need(
            omitted["known"] == base["known"] - value["known"]
            and omitted["lower"] == base["lower"] - value["lower"]
            and omitted["upper"] == base["upper"] - value["upper"],
            "Disjoint gene bounds fail exact additivity",
        )
        nonzero = {key: value for (key, value) in (subset.items()) if (value)}
        record = {
            "gene": gene,
            "nonzero_samples": len(nonzero),
            "nonzero_donors": len({sample_info[key[0]]["donor_id"] for (key) in (nonzero)}),
            "measured_positive_calls": sum(calls[key] == ("MEASURED", 1) for (key) in (nonzero)),
            "measured_zero_calls": sum(calls[key] == ("MEASURED", 0) for (key) in (nonzero)),
            "unknown_calls": sum(calls[key][0] != "MEASURED" for (key) in (nonzero)),
            "known": value["known"],
            "known_positive_coefficients": sum(
                (
                    coefficient
                    for (key, coefficient) in (nonzero.items())
                    if (calls[key] == ("MEASURED", 1) and coefficient > 0)
                ),
                Fraction(),
            ),
            "known_negative_coefficients": sum(
                (
                    coefficient
                    for (key, coefficient) in (nonzero.items())
                    if (calls[key] == ("MEASURED", 1) and coefficient < 0)
                ),
                Fraction(),
            ),
            "unknown_negative": value["unknown_negative"],
            "unknown_positive": value["unknown_positive"],
            "lower": value["lower"],
            "upper": value["upper"],
            "width": value["width"],
            "absolute_coefficient_mass": sum(
                (abs(value) for (value) in (subset.values())), Fraction()
            ),
            "absolute_known_call_mass": sum(
                (
                    abs(coefficient)
                    for (key, coefficient) in (subset.items())
                    if (calls[key] == ("MEASURED", 1))
                ),
                Fraction(),
            ),
            "absolute_net_known_gene_contribution": abs(value["known"]),
            "omitted_lower": omitted["lower"],
            "omitted_upper": omitted["upper"],
            "omitted_known": omitted["known"],
            "omitted_sign": omitted["sign"],
        }
        records.append(record)
    summaries = {}
    for (field) in ((
        "absolute_net_known_gene_contribution",
        "absolute_known_call_mass",
        "absolute_coefficient_mass",
    )):
        ranking = sorted(records, key = lambda row: (-row[field], row["gene"]))
        total = sum((row[field] for (row) in (records)), Fraction())
        cumulative = Fraction()
        for (position, row) in (enumerate(ranking, 1)):
            cumulative += row[field]
            row[field + "_rank"] = position
            row[field + "_cumulative_share"] = cumulative / total if (total) else None
        summaries[field] = {
            "denominator": total,
            "top": {
                str(k): {
                    "gene_ids": [row["gene"] for (row) in (ranking[:k])],
                    "mass": sum((row[field] for (row) in (ranking[:k])), Fraction()),
                    "share": sum((row[field] for (row) in (ranking[:k])), Fraction()) / total
                    if (total)
                    else None,
                }
                for (k) in ((1, 5, 10))
            },
        }
    return {
        "genes": len(records),
        "nonzero_joint_genes": sum(row["nonzero_samples"] > 0 for (row) in (records)),
        "concentration": summaries,
        "omission_sign_counts": dict(Counter(row["omitted_sign"] for (row) in (records))),
        "omissions_losing_negative_identification": [
            row["gene"] for (row) in (records) if (row["omitted_upper"] >= 0)
        ],
    }, records

def analyse_method(method, models, calls, origins):
    cell_rows, ledger, per_model = [], [], {}
    coefficients, genes, sample_info, allocation_mass = {}, set(), {}, {}
    for (model) in (models):
        sample = model["sample_id"]
        sample_info[sample] = model
        arms = model["methods"][method]
        universe = set().union(*(set(value) for (value) in (arms.values())))
        genes.update(universe)
        for (cell, weights) in (arms.items()):
            value, counts = mass(weights, sample, calls)
            per_model[(sample, cell)] = value
            cell_rows.append(
                {
                    "method": method,
                    "level": "model",
                    "model_id": model["model_id"],
                    "sample_id": sample,
                    "donor_id": model["donor_id"],
                    "cell": cell,
                    "models": 1,
                    "donors": 1,
                    "original_weight": model["alpha"],
                    **value,
                    **{
                        key + "_sample_gene_count": counts[key]
                        for (key) in (MASS)
                        if (key not in ("returned", "unknown", "vacant"))
                    },
                }
            )
        for (gene) in (sorted(universe)):
            coefficient = model["alpha"] * (
                arms["G1H1"].get(gene, Fraction()) - arms["G0H0"].get(gene, Fraction())
            )
            coefficients[(sample, gene)] = coefficient
            allocation_mass[(sample, gene)] = model["alpha"] * (
                arms["G1H1"].get(gene, Fraction()) + arms["G0H0"].get(gene, Fraction())
            )
            state, label = calls[(sample, gene)]
            ledger.append(
                {
                    "method": method,
                    "model_id": model["model_id"],
                    "sample_id": sample,
                    "donor_id": model["donor_id"],
                    "gene": gene,
                    "donor_weight": model["alpha"],
                    **{cell: arms[cell].get(gene, Fraction()) for (cell) in (CELLS)},
                    "joint_coefficient": coefficient,
                    "label_state": state,
                    "label_value": label,
                    "label_carriers": origins[(sample, gene)],
                    "known_contribution": coefficient * label if (state == "MEASURED") else 0,
                    "unknown_lower_contribution": min(coefficient, 0)
                    if (state != "MEASURED")
                    else 0,
                    "unknown_upper_contribution": max(coefficient, 0)
                    if (state != "MEASURED")
                    else 0,
                }
            )
    donor_ids = sorted({model["donor_id"] for (model) in (models)})
    for (donor) in (donor_ids):
        included = [model for (model) in (models) if (model["donor_id"] == donor)]
        for (cell) in (CELLS):
            value = weighted(
                [
                    (Fraction(1, len(included)), per_model[(model["sample_id"], cell)])
                    for (model) in (included)
                ]
            )
            cell_rows.append(
                {
                    "method": method,
                    "level": "donor",
                    "model_id": "",
                    "sample_id": "",
                    "donor_id": donor,
                    "cell": cell,
                    "models": len(included),
                    "donors": 1,
                    "original_weight": Fraction(1, len(donor_ids)),
                    **value,
                }
            )
    population = {}
    for (cell) in (CELLS):
        value = weighted(
            [(model["alpha"], per_model[(model["sample_id"], cell)]) for (model) in (models)]
        )
        population[cell] = value
        cell_rows.append(
            {
                "method": method,
                "level": "population",
                "model_id": "",
                "sample_id": "",
                "donor_id": "",
                "cell": cell,
                "models": len(models),
                "donors": len(donor_ids),
                "original_weight": 1,
                **value,
            }
        )
    joint = bounds(coefficients, calls)
    need(
        joint["known"]
        == population["G1H1"]["known_positive"] - population["G0H0"]["known_positive"],
        "Absolute levels do not reproduce known joint contrast",
    )
    flips, flip_rows, curves = flip_analysis(coefficients, calls, sample_info, allocation_mass)
    genes_summary, gene_rows = gene_analysis(coefficients, calls, genes, sample_info)
    for (collection) in ((flip_rows, curves, gene_rows)):
        for (row) in (collection):
            row["method"] = method
    active = [row for (row) in (ledger) if (row["joint_coefficient"])]
    masks = dict(
        Counter(
            row["label_state"]
            + ("_" + str(row["label_value"]) if (row["label_state"] == "MEASURED") else "")
            for (row) in (active)
        )
    )
    return {
        "models": len(models),
        "donors": len(donor_ids),
        "four_cell_levels": population,
        "positive_weight_sample_gene_union": len(ledger),
        "all_arm_label_states": dict(
            Counter(
                row["label_state"]
                + ("_" + str(row["label_value"]) if (row["label_state"] == "MEASURED") else "")
                for (row) in (ledger)
            )
        ),
        "joint_nonzero_sample_genes": len(active),
        "joint_label_states": masks,
        "joint": joint,
        "label_flip": flips,
        "gene_influence": genes_summary,
    }, {
        "cell_levels": cell_rows,
        "sample_gene_ledger": ledger,
        "flip_cells": flip_rows,
        "flip_curves": curves,
        "gene_influence": gene_rows,
    }
