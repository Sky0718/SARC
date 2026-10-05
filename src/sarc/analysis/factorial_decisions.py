from collections import Counter
from fractions import Fraction

SEEDS = [104729, 130363, 155921, 196613, 228017]
ORDER = [
    "joint",
    "graph_H0",
    "graph_H1",
    "joint_minus_graph_H0",
    "joint_minus_graph_H1",
    "interaction",
]
ALIASES = {
    "graph_H0": "graph_H0",
    "graph_H1": "graph_H1",
    "joint_minus_graph_H0": "evidence_G1",
    "joint_minus_graph_H1": "evidence_G0",
    "interaction": "interaction",
}
CELLS = {
    "00": ("11.0", "H0"),
    "01": ("11.0", "H1"),
    "10": ("11.5", "H0"),
    "11": ("11.5", "H1"),
}
MASKS = ["MEASURED_ZERO", "MEASURED_ONE", "EXPLICIT_NULL", "ABSENT"]
SIGNS = [
    "STRICTLY_POSITIVE",
    "STRICTLY_NEGATIVE",
    "IDENTIFIED_ZERO",
    "SIGN_UNRESOLVED",
    "NATIVE_UNAVAILABLE",
]
FIELDS = ["known_contribution", "lower", "upper"]
LINEAGE = {
    "joint": "Mean of five accepted allocation differences W11-W00; both marginal algebraic paths checked",
    **ALIASES,
}

def need(value, message):
    if (not value):
        raise ValueError(message)

def fraction(value):
    need(isinstance(value, str), "Fractions must be canonical strings")
    result = Fraction(value)
    need(str(result) == value, "Noncanonical fraction")
    return result

def terms(mapping):
    return [
        {"sample_id": sample, "gene": gene, "coefficient": str(value)}
        for (((sample, gene), value)) in (sorted(mapping.items()))
        if (value)
    ]

def unterms(rows):
    result = {}
    for (row) in (rows):
        key = (row["sample_id"], row["gene"])
        need(
            all((isinstance(value, str) and value for (value) in (key)))
            and key not in result,
            "Coefficient identity invalid or repeated",
        )
        result[key] = fraction(row["coefficient"])
        need(result[key] != 0, "Zero coefficient must be cancelled")
    return result

def combine(parts):
    result = {}
    for (scale, mapping) in (parts):
        for (key, value) in (mapping.items()):
            result[key] = result.get(key, Fraction()) + scale * value
    return {key: value for ((key, value)) in (result.items()) if (value)}

def coefficient_record(mapping, missing):
    return {
        "available": not missing,
        "coefficients": terms(mapping) if (not missing) else None,
        "missing": missing,
    }

def validate_protocol(protocol, parent):
    need(
        protocol["schema"] == "factorial_decision_closure_protocol_v1"
        and protocol["status"] == "FROZEN_BEFORE_NEW_DERIVATION_EXPLORATORY",
        "Frozen exploratory protocol required",
    )
    need(
        parent["schema"] == "common_domain_factorial_coefficient_freeze_v1"
        and parent["prelabel_freeze"] is True
        and (parent["accounting_complete"] is True),
        "Accepted coefficient freeze required",
    )
    scope = parent["scope"]
    need(
        parent["data_class"] == scope["data_class"]
        and parent["data_class"] in ("SYNTHETIC", "PRODUCTION"),
        "Data class",
    )
    need(
        scope["budget"] == protocol["scope"]["budget"] == 10
        and scope["master_seeds"] == protocol["scope"]["seeds"] == SEEDS
        and (
            scope["transition"] == protocol["scope"]["transition"] == ["11.0", "11.5"]
        ),
        "Frozen budget, seeds or transition differs",
    )
    definitions = {
        "joint": "W11-W00",
        "graph_H0": "W10-W00",
        "graph_H1": "W11-W01",
        "joint_minus_graph_H0": "W11-W10",
        "joint_minus_graph_H1": "W01-W00",
        "interaction": "W11-W01-W10+W00",
    }
    need(
        {name: value["definition"] for ((name, value)) in (protocol["contrasts"].items())}
        == definitions,
        "Frozen contrast definitions differ",
    )
    models = scope["models"]
    need(
        models
        and all(
            (
                all(
                    (
                        isinstance(model[key], str) and model[key]
                        for (key) in (("model_id", "donor_id", "sample_id"))
                    )
                )
                for (model) in (models)
            )
        ),
        "Complete model identities",
    )
    need(
        len({model["model_id"] for (model) in (models)}) == len(models),
        "Repeated model",
    )
    need(
        len({model["sample_id"] for (model) in (models)}) == len(models),
        "Sample identities must be disjoint across models and donors",
    )
    donors = len({model["donor_id"] for (model) in (models)})
    counts = {
        "models": len(models),
        "donors": donors,
        "cells": len(models) * 4,
        "seed_states": len(models) * 20,
    }
    need(
        counts
        == parent["counts"]
        == scope["expected"]
        == {key: protocol["scope"][key] for (key) in (counts)},
        "Full population coverage differs",
    )
    need(donors >= 6, "All k=1..5 sensitivities require at least six donors")
    if (parent["data_class"] == "PRODUCTION"):
        need((len(models), donors) == (85, 83), "All85models and83donors required")
    else:
        need(
            all(
                (
                    all(
                        (
                            model[key].startswith("SYN_")
                            for (key) in (("model_id", "donor_id", "sample_id"))
                        )
                    )
                    for (model) in (models)
                )
            ),
            "Synthetic identity prefix",
        )
    return (sorted(models, key = lambda model: model["model_id"]), counts)

def freeze(protocol, parent):
    models, counts = validate_protocol(protocol, parent)
    mids = {model["model_id"] for (model) in (models)}
    old = {}
    for (row) in (parent["model_contrasts"]):
        key = (row["model_id"], row["contrast"])
        need(
            key not in old and key[0] in mids and (key[1] in set(ALIASES.values())),
            "Duplicate or extra inherited model contrast",
        )
        old[key] = row
    need(len(old) == len(models) * 5, "Missing inherited model contrast")
    allocations = {}
    inverse = {value: key for ((key, value)) in (CELLS.items())}
    for (row) in (parent["allocations"]):
        need(
            (row["graph"], row["evidence"]) in inverse
            and type(row["master_seed"]) is int,
            "Unknown allocation state",
        )
        key = (
            row["model_id"],
            inverse[row["graph"], row["evidence"]],
            row["master_seed"],
        )
        need(
            key[0] in mids and key[2] in SEEDS and (key not in allocations),
            "Duplicate or extra allocation",
        )
        status = row["status"]
        need(
            status in ("SUCCESS", "SUCCESS_EMPTY", "UNAVAILABLE", "FAILED"),
            "Explicit allocation status required",
        )
        if (status in ("UNAVAILABLE", "FAILED")):
            need(
                row["weights"] is None
                and row["capacity"] is None
                and row.get("reason"),
                "Unavailable allocation must preserve missingness",
            )
            weights = None
        else:
            weights = {}
            for (gene, value) in (row["weights"]):
                need(
                    isinstance(gene, str) and gene and (gene not in weights),
                    "Allocation gene identity",
                )
                weights[gene] = fraction(value)
                need(0 <= weights[gene] <= 1, "Allocation weight range")
            need(
                sum(weights.values(), Fraction()) == fraction(row["capacity"])
                and 0 <= fraction(row["capacity"]) <= 10,
                "Allocation capacity differs",
            )
            if (status == "SUCCESS_EMPTY"):
                need(
                    not weights and row.get("reason"),
                    "Successful empty state must remain explicit",
                )
        allocations[key] = (row, weights)
    need(
        len(allocations) == counts["seed_states"],
        "Missing complete allocation population",
    )
    output = []
    for (model) in (models):
        mid, sample = (model["model_id"], model["sample_id"])
        capacities, seed_capacities, seed_joint = ({}, {}, [])
        for (cell) in (CELLS):
            records = [allocations[mid, cell, seed][0] for (seed) in (SEEDS)]
            seed_capacities[cell] = [
                {
                    key: row.get(key)
                    for (key) in (("master_seed", "status", "reason", "capacity"))
                }
                for (row) in (records)
            ]
            capacities[cell] = (
                str(
                    sum((fraction(row["capacity"]) for (row) in (records)), Fraction())
                    / 5
                )
                if (all((row["capacity"] is not None for (row) in (records))))
                else None
            )
        for (seed) in (SEEDS):
            selected = [
                (cell, allocations[mid, cell, seed][1]) for (cell) in (("11", "00"))
            ]
            missing = [cell for ((cell, weights)) in (selected) if (weights is None)]
            mapping = (
                {}
                if (missing)
                else combine(
                    [
                        (
                            Fraction(1 if (cell == "11") else -1),
                            {
                                (sample, gene): value
                                for ((gene, value)) in (weights.items())
                            },
                        )
                        for ((cell, weights)) in (selected)
                    ]
                )
            )
            seed_joint.append(
                {"master_seed": seed, **coefficient_record(mapping, missing)}
            )
        missing = [
            row["master_seed"] for (row) in (seed_joint) if (not row["available"])
        ]
        joint = (
            {}
            if (missing)
            else combine(
                [
                    (Fraction(1, 5), unterms(row["coefficients"]))
                    for (row) in (seed_joint)
                ]
            )
        )
        contrasts = {"joint": coefficient_record(joint, missing)}
        for (name, alias) in (ALIASES.items()):
            row = old[mid, alias]
            need(
                row["sample_id"] == sample
                and row["donor_id"] == model["donor_id"]
                and (type(row["available"]) is bool),
                "Inherited model identity differs",
            )
            mapping = unterms(row["coefficients"]) if (row["available"]) else {}
            need(
                all((key[0] == sample for (key) in (mapping))),
                "Model coefficients contain another sample",
            )
            need(
                row["available"]
                or (row["coefficients"] is None and row.get("missing")),
                "Unavailable contrast must preserve missing states",
            )
            contrasts[name] = coefficient_record(
                mapping, [] if (row["available"]) else row["missing"]
            )
        paths_match = None
        if (all((contrasts[name]["available"] for (name) in (ORDER)))):
            path0 = combine(
                [
                    (Fraction(1), unterms(contrasts[name]["coefficients"]))
                    for (name) in (("graph_H0", "joint_minus_graph_H0"))
                ]
            )
            path1 = combine(
                [
                    (Fraction(1), unterms(contrasts[name]["coefficients"]))
                    for (name) in (("graph_H1", "joint_minus_graph_H1"))
                ]
            )
            need(
                joint == path0 == path1,
                "Joint algebraic paths differ from accepted allocations",
            )
            interaction = combine(
                [
                    (Fraction(1), unterms(contrasts["graph_H1"]["coefficients"])),
                    (Fraction(-1), unterms(contrasts["graph_H0"]["coefficients"])),
                ]
            )
            need(
                interaction == unterms(contrasts["interaction"]["coefficients"]),
                "Accepted interaction identity differs",
            )
            paths_match = True
        output.append(
            {key: model[key] for (key) in (("model_id", "donor_id", "sample_id"))}
            | {
                "capacities": capacities,
                "seed_capacities": seed_capacities,
                "seed_joint": seed_joint,
                "contrasts": contrasts,
                "joint_algebraic_paths_match": paths_match,
            }
        )
    requested = sorted(
        {
            (row["sample_id"], row["gene"])
            for (model) in (output)
            for (contrast) in (model["contrasts"].values())
            if (contrast["available"])
            for (row) in (contrast["coefficients"])
        }
    )
    accepted = {
        (row["sample_id"], row["gene"]) for (row) in (parent["requested_labels"])
    }
    need(
        len(accepted) == len(parent["requested_labels"]) and set(requested) <= accepted,
        "New identities outside accepted selected-call set",
    )
    return {
        "schema": "factorial_decision_coefficient_freeze_v1",
        "data_class": parent["data_class"],
        "counts": counts,
        "budget": 10,
        "master_seeds": SEEDS,
        "contrast_order": ORDER,
        "requested_labels": [
            {"sample_id": sample, "gene": gene} for ((sample, gene)) in (requested)
        ],
        "models": output,
        "prelabel_freeze": True,
        "inherited_input_bindings": parent.get("input_bindings", {}),
        "coefficient_lineage": LINEAGE,
    }

def sign(lower, upper):
    return (
        "STRICTLY_POSITIVE"
        if (lower > 0)
        else "STRICTLY_NEGATIVE"
        if (upper < 0)
        else "IDENTIFIED_ZERO"
        if (lower == upper == 0)
        else "SIGN_UNRESOLVED"
    )

def bound_record(known, lower, upper, masks):
    need(lower <= known <= upper, "Invalid exact identification interval")
    return {
        "known_contribution": str(known),
        "lower": str(lower),
        "upper": str(upper),
        "per_offered_slot": {
            key: str(value / 10)
            for ((key, value)) in (zip(FIELDS, (known, lower, upper)))
        },
        "masks": {mask: masks.get(mask, 0) for (mask) in (MASKS)},
        "sign": sign(lower, upper),
    }

def bounds(mapping, calls):
    known, low, high = (Fraction(), Fraction(), Fraction())
    masks = Counter()
    for (key, value) in (mapping.items()):
        state, call = calls[key]
        mask = (
            "MEASURED_ONE"
            if (state == "MEASURED" and call == 1)
            else "MEASURED_ZERO"
            if (state == "MEASURED")
            else state
        )
        masks[mask] += 1
        if (state == "MEASURED"):
            known += value * call
        else:
            low += min(0, value)
            high += max(0, value)
    return bound_record(known, known + low, known + high, masks)

def disjoint_bounds(parts):
    known, lower, upper = (Fraction(), Fraction(), Fraction())
    masks = Counter()
    for (scale, value) in (parts):
        if (not scale):
            continue
        known += scale * fraction(value["known_contribution"])
        left, right = (
            scale * fraction(value["lower"]),
            scale * fraction(value["upper"]),
        )
        lower += min(left, right)
        upper += max(left, right)
        masks.update(value["masks"])
    return bound_record(known, lower, upper, masks)

def comparisons(contrasts):
    output = {}
    for (graph) in (("graph_H0", "graph_H1")):
        joint = contrasts["joint"]["value"]
        other = contrasts[graph]["value"]
        js, gs = (
            joint["sign"] if (joint) else "NATIVE_UNAVAILABLE",
            other["sign"] if (other) else "NATIVE_UNAVAILABLE",
        )
        classification = "UNRESOLVED"
        strict = {"STRICTLY_POSITIVE", "STRICTLY_NEGATIVE"}
        if (js in strict and gs in strict):
            classification = (
                "SAME_ROBUST_DIRECTION" if (js == gs) else "ROBUST_OPPOSITE_DIRECTION"
            )
        elif (js == gs == "IDENTIFIED_ZERO"):
            classification = "BOTH_IDENTIFIED_ZERO"
        elif (js == "IDENTIFIED_ZERO" and gs != "NATIVE_UNAVAILABLE"):
            classification = "JOINT_IDENTIFIED_ZERO"
        elif (gs == "IDENTIFIED_ZERO" and js != "NATIVE_UNAVAILABLE"):
            classification = "GRAPH_IDENTIFIED_ZERO"
        output[graph] = {
            "classification": classification,
            "joint_sign": js,
            "graph_sign": gs,
            "paired_difference": contrasts["joint_minus_" + graph]["value"],
        }
    return output

def mean_capacities(rows):
    return {
        cell: str(
            sum((fraction(row["capacities"][cell]) for (row) in (rows)), Fraction())
            / len(rows)
        )
        if (all((row["capacities"][cell] is not None for (row) in (rows))))
        else None
        for (cell) in (CELLS)
    }

def evaluate(frozen, labels, parent_readout = None):
    need(
        frozen["schema"] == "factorial_decision_coefficient_freeze_v1"
        and frozen["prelabel_freeze"] is True,
        "New coefficient freeze required before labels",
    )
    need(
        labels["schema"] == "common_domain_factorial_selected_labels_v1"
        and labels["data_class"] == frozen["data_class"]
        and (labels["label_contract"] == "final_Table5_binary_original_mapping"),
        "Accepted selected label semantics differ",
    )
    all_calls = {}
    for (row) in (labels["calls"]):
        key = (row["sample_id"], row["gene"])
        need(
            key not in all_calls
            and all((isinstance(part, str) and part for (part) in (key))),
            "Repeated or invalid selected call",
        )
        state, value = (row["state"], row["value"])
        need(state in ("MEASURED", "EXPLICIT_NULL", "ABSENT"), "Unknown mask state")
        need(
            type(value) is int and value in (0, 1)
            if (state == "MEASURED")
            else value is None,
            "Unknown/null and measured0/1 must remain distinct",
        )
        all_calls[key] = (state, value)
    requested = {
        (row["sample_id"], row["gene"]) for (row) in (frozen["requested_labels"])
    }
    need(
        len(requested) == len(frozen["requested_labels"])
        and requested <= set(all_calls),
        "Missing selected-call coverage",
    )
    calls = {key: all_calls[key] for (key) in (requested)}
    output = {
        "schema": "factorial_decision_readout_v1",
        "data_class": frozen["data_class"],
        "counts": frozen["counts"],
        "budget": 10,
        "master_seeds": SEEDS,
        "contrast_order": ORDER,
        "selected_calls": [
            {
                "sample_id": sample,
                "gene": gene,
                "state": calls[sample, gene][0],
                "value": calls[sample, gene][1],
            }
            for ((sample, gene)) in (sorted(calls))
        ],
        "models": [],
        "donors": [],
        "population": None,
        "leave_one_donor_out": [],
        "worst_case_removals": [],
        "donor_sign_counts": {},
    }
    roster = frozen["models"]
    need(
        len({model["sample_id"] for (model) in (roster)}) == len(roster),
        "Disjoint support requirement failed",
    )
    by_donor = {}
    for (model) in (roster):
        by_donor.setdefault(model["donor_id"], []).append(model)
        contrasts = {
            name: {
                "available": model["contrasts"][name]["available"],
                "value": bounds(
                    unterms(model["contrasts"][name]["coefficients"]), calls
                )
                if (model["contrasts"][name]["available"])
                else None,
            }
            for (name) in (ORDER)
        }
        output["models"].append(
            {
                key: model[key]
                for (key) in (("model_id", "donor_id", "sample_id", "capacities"))
            }
            | {"contrasts": contrasts, "comparisons": comparisons(contrasts)}
        )
    donor_maps = {}
    for (donor, models) in (sorted(by_donor.items())):
        contrasts, donor_maps[donor] = ({}, {})
        for (name) in (ORDER):
            available = all(
                (model["contrasts"][name]["available"] for (model) in (models))
            )
            mapping = (
                combine(
                    [
                        (
                            Fraction(1, len(models)),
                            unterms(model["contrasts"][name]["coefficients"]),
                        )
                        for (model) in (models)
                    ]
                )
                if (available)
                else None
            )
            donor_maps[donor][name] = mapping
            contrasts[name] = {
                "available": available,
                "value": bounds(mapping, calls) if (available) else None,
            }
        output["donors"].append(
            {
                "donor_id": donor,
                "model_ids": sorted((model["model_id"] for (model) in (models))),
                "capacities": mean_capacities(models),
                "contrasts": contrasts,
                "comparisons": comparisons(contrasts),
            }
        )
    donors = output["donors"]
    count = len(donors)
    need(
        count == frozen["counts"]["donors"]
        and len(roster) == frozen["counts"]["models"]
        and (count >= 6),
        "Full readout population required",
    )
    population = {}
    for (name) in (ORDER):
        available = all((row["contrasts"][name]["available"] for (row) in (donors)))
        mapping = (
            combine(
                [
                    (Fraction(1, count), donor_maps[row["donor_id"]][name])
                    for (row) in (donors)
                ]
            )
            if (available)
            else None
        )
        population[name] = {
            "available": available,
            "value": bounds(mapping, calls) if (available) else None,
        }
        output["donor_sign_counts"][name] = {
            state: sum(
                (
                    (
                        row["contrasts"][name]["value"]["sign"]
                        if (row["contrasts"][name]["available"])
                        else "NATIVE_UNAVAILABLE"
                    )
                    == state
                    for (row) in (donors)
                )
            )
            for (state) in (SIGNS)
        }
    output["population"] = {
        "donor_ids": [row["donor_id"] for (row) in (donors)],
        "capacities": mean_capacities(donors),
        "contrasts": population,
        "comparisons": comparisons(population),
    }
    for (removed) in (donors):
        remaining = [
            row for (row) in (donors) if (row["donor_id"] != removed["donor_id"])
        ]
        contrasts = {}
        for (name) in (ORDER):
            available = all(
                (row["contrasts"][name]["available"] for (row) in (remaining))
            )
            value = (
                disjoint_bounds(
                    [
                        (Fraction(1, count - 1), row["contrasts"][name]["value"])
                        for (row) in (remaining)
                    ]
                )
                if (available)
                else None
            )
            full = population[name]["value"]
            difference = (
                disjoint_bounds(
                    [
                        (
                            Fraction(1, count * (count - 1)),
                            row["contrasts"][name]["value"],
                        )
                        for (row) in (remaining)
                    ]
                    + [(Fraction(-1, count), removed["contrasts"][name]["value"])]
                )
                if (full is not None)
                else None
            )
            shifts = (
                {
                    key: str(fraction(value[key]) - fraction(full[key]))
                    for (key) in (FIELDS)
                }
                if (value is not None and full is not None)
                else None
            )
            contrasts[name] = {
                "available": available,
                "value": value,
                "paired_difference_from_full": difference,
                "endpoint_shifts_from_full": shifts,
            }
        output["leave_one_donor_out"].append(
            {
                "removed_donor_id": removed["donor_id"],
                "remaining_donors": count - 1,
                "contrasts": contrasts,
                "comparisons": comparisons(contrasts),
            }
        )
    for (k) in (range(1, 6)):
        contrasts = {}
        for (name) in (ORDER):
            available = population[name]["available"]
            record = {"available": available}
            for (label, field, descending) in ([
                ("minimum_lower", "lower", True),
                ("maximum_upper", "upper", False),
                ("minimum_known", "known_contribution", True),
                ("maximum_known", "known_contribution", False),
            ]):
                if (not available):
                    record[label] = None
                    continue
                selected = sorted(
                    donors,
                    key = lambda row: (
                        (-1 if (descending) else 1)
                        * fraction(row["contrasts"][name]["value"][field]),
                        row["donor_id"],
                    ),
                )[:k]
                removed_ids = sorted((row["donor_id"] for (row) in (selected)))
                remaining = [
                    row for (row) in (donors) if (row["donor_id"] not in removed_ids)
                ]
                record[label] = {
                    "removed_donor_ids": removed_ids,
                    "value": disjoint_bounds(
                        [
                            (Fraction(1, count - k), row["contrasts"][name]["value"])
                            for (row) in (remaining)
                        ]
                    ),
                }
            contrasts[name] = record
        output["worst_case_removals"].append(
            {"k": k, "remaining_donors": count - k, "contrasts": contrasts}
        )
    if (parent_readout is not None):
        validate_reuse(output, parent_readout)
    return output

def validate_reuse(output, parent):
    need(
        parent["schema"] == "common_domain_factorial_readout_v1"
        and parent["data_class"] == output["data_class"]
        and (parent["counts"] == output["counts"]),
        "Accepted readout identity differs",
    )
    for (section, old_section) in ([
        ("models", "model_contrasts"),
        ("population", "population_contrasts"),
    ]):
        old = {}
        for (row) in (parent[old_section]):
            key = (row.get("model_id"), row["contrast"])
            need(key not in old, "Duplicate accepted readout row")
            old[key] = row
        rows = output["models"] if (section == "models") else [output["population"]]
        need(len(old) == len(rows) * 5, "Accepted readout contrast coverage differs")
        for (row) in (rows):
            for (name, alias) in (ALIASES.items()):
                previous = old[row.get("model_id"), alias]
                current = row["contrasts"][name]
                need(
                    previous["available"] == current["available"],
                    "Reuse availability differs",
                )
                if (current["available"]):
                    for (key) in (FIELDS + ["per_offered_slot", "sign"]):
                        need(
                            previous["value"][key] == current["value"][key],
                            "Accepted result reuse differs: " + key,
                        )
