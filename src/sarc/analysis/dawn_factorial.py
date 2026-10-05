import math
from collections import Counter
from fractions import Fraction

CELLS = ["G0_H0", "G1_H0", "G0_H1", "G1_H1"]
SEEDS = [104729, 130363, 155921, 196613, 228017]
CONTRASTS = {
    "graph_H0": {"G1_H0": 1, "G0_H0": -1},
    "graph_H1": {"G1_H1": 1, "G0_H1": -1},
    "evidence_G0": {"G0_H1": 1, "G0_H0": -1},
    "evidence_G1": {"G1_H1": 1, "G1_H0": -1},
    "interaction": {"G1_H1": 1, "G1_H0": -1, "G0_H1": -1, "G0_H0": 1},
    "joint": {"G1_H1": 1, "G0_H0": -1},
}
SIGNS = [
    "STRICTLY_POSITIVE",
    "STRICTLY_NEGATIVE",
    "IDENTIFIED_ZERO",
    "SIGN_UNRESOLVED",
    "NATIVE_UNAVAILABLE",
]

def need(value, message):
    if (not value):
        raise ValueError(message)

def fraction(value):
    need(isinstance(value, str), "Fraction must be string")
    answer = Fraction(value)
    need(str(answer) == value, "Noncanonical fraction")
    return answer

def terms(mapping):
    return [
        {"sample_id": sample, "gene": gene, "coefficient": str(value)}
        for (((sample, gene), value)) in (sorted(mapping.items()))
        if (value)
    ]

def unterms(values):
    output = {}
    for (row) in (values):
        key = (row["sample_id"], row["gene"])
        need(key not in output, "Repeated coefficient identity")
        output[key] = fraction(row["coefficient"])
        need(output[key] != 0, "Zero coefficient must cancel")
    return output

def combine(parts):
    output = {}
    for (scale, part) in (parts):
        for (key, value) in (part.items()):
            output[key] = output.get(key, Fraction()) + scale * value
    return {key: value for ((key, value)) in (output.items()) if (value)}

def allocate(scores, eligible):
    need(
        len({gene for ((gene, score, percentile)) in (scores)}) == len(scores),
        "Repeated native gene",
    )
    need(
        all(
            (
                math.isfinite(score)
                and score >= 0
                and math.isfinite(percentile)
                and (0 < percentile < 100)
                for ((gene, score, percentile)) in (scores)
            )
        ),
        "Native score or percentile invalid",
    )
    for (left, right) in (zip(scores, scores[1:])):
        need(
            left[1] >= right[1] and left[2] >= right[2], "Native ranking not descending"
        )
        need(
            (left[1] == right[1]) == (left[2] == right[2]),
            "Score and native percentile tie partitions differ",
        )
    eligible = set(eligible)
    selected = [row for (row) in (scores) if (row[0] in eligible)]
    weights, offset = ({}, 0)
    while (offset < len(selected) and offset < 10):
        end = offset + 1
        while (end < len(selected) and selected[end][2] == selected[offset][2]):
            end += 1
        weight = Fraction(min(end - offset, 10 - offset), end - offset)
        for (gene, score, percentile) in (selected[offset:end]):
            weights[gene] = weight
        offset = end
    need(
        sum(weights.values(), Fraction()) == min(10, len(selected)),
        "Offered budget allocation mismatch",
    )
    return (weights, len(selected))

def coefficient_record(mapping, missing):
    if (missing):
        return {
            "available": False,
            "coefficients": None,
            "missing": missing,
            "positive_mass": None,
            "negative_mass": None,
            "half_l1": None,
            "capacity_contrast": None,
        }
    positive = sum((value for (value) in (mapping.values()) if (value > 0)), Fraction())
    negative = -sum(
        (value for (value) in (mapping.values()) if (value < 0)), Fraction()
    )
    return {
        "available": True,
        "coefficients": terms(mapping),
        "missing": [],
        "positive_mass": str(positive),
        "negative_mass": str(negative),
        "half_l1": str((positive + negative) / 2),
        "capacity_contrast": str(positive - negative),
    }

def validate_scope(scope, data_class):
    need(
        data_class in ("PRODUCTION", "SYNTHETIC") and scope["data_class"] == data_class,
        "Data class",
    )
    need(
        scope["budget"] == 10
        and scope["master_seeds"] == SEEDS
        and (scope["label_contract"] == "final_Table5_binary_original_mapping"),
        "Frozen query contract",
    )
    models = scope["models"]
    need(
        models
        and len({m["sample_id"] for (m) in (models)}) == len(models)
        and (len({m["model_id"] for (m) in (models)}) == len(models)),
        "Complete disjoint model identities",
    )
    for (model) in (models):
        need(
            len(model["root_genes"]) == len(set(model["root_genes"]))
            and len(model["eligible_genes"]) == len(set(model["eligible_genes"]))
            and (set(model["eligible_genes"]) <= set(model["root_genes"])),
            "Query eligibility invalid",
        )
    donors = len({model["donor_id"] for (model) in (models)})
    if (data_class == "PRODUCTION"):
        need((len(models), donors) == (85, 83), "Full85model83donor scope required")
    else:
        need(
            all(
                (
                    all(
                        (
                            model[key].startswith("SYN_")
                            for (key) in (("model_id", "sample_id", "donor_id"))
                        )
                    )
                    for (model) in (models)
                )
            ),
            "Fixture identities required",
        )
    need(donors > 1, "Donor leaveout requires at least two donors")
    return {
        "models": len(models),
        "donors": donors,
        "physical_model_cells": len(models) * 4,
        "logical_seed_references": len(models) * 20,
    }

def build_freeze(scope, raw_cells, data_class, input_bindings = None):
    counts = validate_scope(scope, data_class)
    samples = {model["sample_id"] for (model) in (scope["models"])}
    need(
        set(raw_cells) == set(CELLS)
        and all((set(raw_cells[cell]) == samples for (cell) in (CELLS))),
        "Full model-cell coverage",
    )
    output, requests = ([], set())
    for (model) in (scope["models"]):
        sample = model["sample_id"]
        allocations, maps = ({}, {})
        for (cell) in (CELLS):
            native = raw_cells[cell][sample]
            status = native["status"]
            need(
                status in ("SUCCESS", "SUCCESS_EMPTY", "FAILED", "UNAVAILABLE"),
                "Unknown or unresolved alias state",
            )
            scores = native["scores"]
            if (status in ("FAILED", "UNAVAILABLE")):
                need(
                    not scores and native.get("reason"),
                    "Unavailable requires reason and no rank",
                )
                weight, eligible_count, capacity, vacancy = (None, None, None, None)
            else:
                need(
                    (status == "SUCCESS") == bool(scores),
                    "Native status and rank availability disagree",
                )
                need(
                    set((row[0] for (row) in (scores))) == set(model["root_genes"]),
                    "Dawn must retain every bound common mutation root",
                )
                weight, eligible_count = allocate(scores, model["eligible_genes"])
                capacity = str(sum(weight.values(), Fraction()))
                vacancy = str(10 - fraction(capacity))
                requests.update(((sample, gene) for (gene) in (weight)))
            maps[cell] = (
                None
                if (weight is None)
                else {(sample, gene): value for ((gene, value)) in (weight.items())}
            )
            allocations[cell] = {
                key: native.get(key)
                for (key) in (
                    (
                        "status",
                        "reason",
                        "iterations",
                        "terminal_update_norm",
                        "termination_reason",
                        "fixed_point_error_certified",
                    )
                )
            }
            allocations[cell].update(
                weights = None
                if (weight is None)
                else [[gene, str(value)] for ((gene, value)) in (sorted(weight.items()))],
                capacity = capacity,
                vacancies = vacancy,
                emitted_count = len(scores) if (weight is not None) else None,
                eligible_emitted_count = eligible_count,
                selected_gene_count = len(weight) if (weight is not None) else None,
                seed_references = SEEDS,
                seed_semantics = "Five deterministic aliases, not independent replicates",
            )
        contrasts = {}
        for (name, specification) in (CONTRASTS.items()):
            missing = [cell for (cell) in (specification) if (maps[cell] is None)]
            mapping = (
                {}
                if (missing)
                else combine(
                    [
                        (Fraction(scale), maps[cell])
                        for ((cell, scale)) in (specification.items())
                    ]
                )
            )
            contrasts[name] = coefficient_record(mapping, missing)
        output.append(
            {key: model[key] for (key) in (("model_id", "sample_id", "donor_id"))}
            | {"cells": allocations, "contrasts": contrasts}
        )
    return {
        "schema": "dawn_factorial_allocation_freeze_v1",
        "data_class": data_class,
        "prelabel_freeze": True,
        "counts": counts,
        "budget": 10,
        "master_seeds": SEEDS,
        "contrast_definitions": CONTRASTS,
        "models": output,
        "requested_labels": [
            {"sample_id": sample, "gene": gene} for ((sample, gene)) in (sorted(requests))
        ],
        "input_bindings": input_bindings or {},
        "label_contract": scope["label_contract"],
        "table5_source": scope.get("table5_source"),
        "exposure": "Previously exposed assay; post-hoc exploratory second-method extension, not independent confirmation",
        "unknown_identity": "Literal sample_id and gene shared across all cells and deterministic seed aliases",
        "rank_precision": "Native %.17g score and native percentile parsed as binary64; no rounding or tolerance ties",
        "eligibility": "Inherited frozen molecular query eligibility before budget allocation; no assay-label filtering or refill",
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

def bounds(mapping, calls):
    known, lower, upper = (Fraction(), Fraction(), Fraction())
    masks = Counter()
    for (key, coefficient) in (mapping.items()):
        state, value = calls[key]
        masks[
            "MEASURED_ONE"
            if (state == "MEASURED" and value)
            else "MEASURED_ZERO"
            if (state == "MEASURED")
            else state
        ] += 1
        if (state == "MEASURED"):
            known += coefficient * value
        else:
            lower += min(coefficient, 0)
            upper += max(coefficient, 0)
    lower += known
    upper += known
    return {
        "known_contribution": str(known),
        "lower": str(lower),
        "upper": str(upper),
        "per_offered_slot": {
            "known_contribution": str(known / 10),
            "lower": str(lower / 10),
            "upper": str(upper / 10),
        },
        "sign": sign(lower, upper),
        "masks": dict(sorted(masks.items())),
    }

def cell_mass(weights, sample, calls):
    if (weights is None):
        return None
    masses = {
        key: Fraction()
        for (key) in (
            (
                "returned",
                "assessed",
                "positive",
                "negative",
                "unknown",
                "explicit_null",
                "absent",
            )
        )
    }
    for (gene, value) in (weights):
        weight = fraction(value)
        state, label = calls[sample, gene]
        masses["returned"] += weight
        if (state == "MEASURED"):
            masses["assessed"] += weight
            masses["positive" if (label) else "negative"] += weight
        else:
            masses["unknown"] += weight
            masses["explicit_null" if (state == "EXPLICIT_NULL") else "absent"] += (
                weight
            )
    masses["vacancies"] = 10 - masses["returned"]
    return {key: str(value) for ((key, value)) in (masses.items())}

def average_mass(records):
    if (any((record is None for (record) in (records)))):
        return None
    return {
        key: str(
            sum((fraction(record[key]) for (record) in (records)), Fraction())
            / len(records)
        )
        for (key) in (records[0])
    }

def evaluated_contrast(record, calls):
    return {
        "available": record["available"],
        "value": bounds(unterms(record["coefficients"]), calls)
        if (record["available"])
        else None,
        **{
            key: record[key]
            for (key) in (
                ("positive_mass", "negative_mass", "half_l1", "capacity_contrast")
            )
        },
        "missing": record["missing"],
    }

def aggregate(records, calls, ids = None):
    contrasts = {}
    for (name) in (CONTRASTS):
        missing = [
            index
            for ((index, row)) in (enumerate(records))
            if (not row["contrasts"][name]["available"])
        ]
        mapping = (
            {}
            if (missing)
            else combine(
                [
                    (
                        Fraction(1, len(records)),
                        unterms(row["contrasts"][name]["coefficients"]),
                    )
                    for (row) in (records)
                ]
            )
        )
        contrasts[name] = coefficient_record(mapping, missing)
    return contrasts

def evaluate(frozen, selected):
    need(
        frozen.get("table5_source")
        and frozen["table5_source"] == selected.get("source_binding"),
        "Allocation and labels require the same assay source",
    )
    need(
        frozen["schema"] == "dawn_factorial_allocation_freeze_v1"
        and frozen["prelabel_freeze"],
        "Prior allocation freeze required",
    )
    need(
        selected["label_contract"] == frozen["label_contract"]
        and selected["data_class"] == frozen["data_class"],
        "Label contract differs",
    )
    calls = {}
    for (row) in (selected["calls"]):
        key = (row["sample_id"], row["gene"])
        state, value = (row["state"], row["value"])
        need(
            key not in calls and state in ("MEASURED", "EXPLICIT_NULL", "ABSENT"),
            "Duplicate or malformed selected label",
        )
        need(
            type(value) is int and value in (0, 1)
            if (state == "MEASURED")
            else value is None,
            "Measured/unknown semantics differ",
        )
        calls[key] = (state, value)
    need(
        set(calls)
        == {(row["sample_id"], row["gene"]) for (row) in (frozen["requested_labels"])},
        "Selected-only label coverage differs",
    )
    output = {
        "schema": "dawn_factorial_readout_v1",
        "data_class": frozen["data_class"],
        "counts": frozen["counts"],
        "budget": 10,
        "master_seeds": SEEDS,
        "models": [],
        "donors": [],
        "population": {},
        "leave_one_donor_out": [],
        "donor_sign_counts": {},
        "bounds_type": "Sharp shared-binary-unknown identification intervals, not confidence intervals or p-values",
        "scientifically_complete": all(
            (
                cell["status"] in ("SUCCESS", "SUCCESS_EMPTY")
                for (model) in (frozen["models"])
                for (cell) in (model["cells"].values())
            )
        ),
        "seed_semantics": "Five deterministic labels reference one physical observation per model and cell; no seed precision",
        "exposure": frozen["exposure"],
        "selected_call_count": len(calls),
    }
    donor_models, model_out = ({}, {})
    for (model) in (frozen["models"]):
        donor_models.setdefault(model["donor_id"], []).append(model)
        cell_values = {
            cell: {
                "status": row["status"],
                "mass": cell_mass(row["weights"], model["sample_id"], calls),
            }
            for ((cell, row)) in (model["cells"].items())
        }
        record = {
            key: model[key] for (key) in (("model_id", "sample_id", "donor_id"))
        } | {
            "cells": cell_values,
            "contrasts": {
                name: evaluated_contrast(row, calls)
                for ((name, row)) in (model["contrasts"].items())
            },
        }
        output["models"].append(record)
        model_out[model["model_id"]] = record
    donor_freezes = []
    for (donor, models) in (sorted(donor_models.items())):
        contrasts = aggregate(models, calls)
        donor_freezes.append({"donor_id": donor, "contrasts": contrasts})
        output["donors"].append(
            {
                "donor_id": donor,
                "model_ids": [model["model_id"] for (model) in (models)],
                "cells": {
                    cell: {
                        "mass": average_mass(
                            [
                                model_out[model["model_id"]]["cells"][cell]["mass"]
                                for (model) in (models)
                            ]
                        )
                    }
                    for (cell) in (CELLS)
                },
                "contrasts": {
                    name: evaluated_contrast(row, calls)
                    for ((name, row)) in (contrasts.items())
                },
            }
        )
    population = aggregate(donor_freezes, calls)
    output["population"] = {
        "donor_ids": [row["donor_id"] for (row) in (donor_freezes)],
        "cells": {
            cell: {
                "mass": average_mass(
                    [row["cells"][cell]["mass"] for (row) in (output["donors"])]
                )
            }
            for (cell) in (CELLS)
        },
        "contrasts": {
            name: evaluated_contrast(row, calls) for ((name, row)) in (population.items())
        },
    }
    for (name) in (CONTRASTS):
        output["donor_sign_counts"][name] = dict.fromkeys(SIGNS, 0)
        for (row) in (output["donors"]):
            contrast = row["contrasts"][name]
            output["donor_sign_counts"][name][
                contrast["value"]["sign"]
                if (contrast["available"])
                else "NATIVE_UNAVAILABLE"
            ] += 1
    for (removed) in (donor_freezes):
        remaining = [
            row for (row) in (donor_freezes) if (row["donor_id"] != removed["donor_id"])
        ]
        contrasts = aggregate(remaining, calls)
        output["leave_one_donor_out"].append(
            {
                "removed_donor_id": removed["donor_id"],
                "remaining_donors": len(remaining),
                "contrasts": {
                    name: evaluated_contrast(row, calls)
                    for ((name, row)) in (contrasts.items())
                },
            }
        )
    output["leaveout_summary"] = {}
    for (name) in (CONTRASTS):
        contrasts = [
            row["contrasts"][name] for (row) in (output["leave_one_donor_out"])
        ]
        available = [row["value"] for (row) in (contrasts) if (row["available"])]
        output["leaveout_summary"][name] = {
            "available": len(available),
            "unavailable": len(contrasts) - len(available),
            "sign_counts": dict(Counter((row["sign"] for (row) in (available)))),
            "ranges": {
                field: [
                    str(min((fraction(row[field]) for (row) in (available)))),
                    str(max((fraction(row[field]) for (row) in (available)))),
                ]
                for (field) in (("known_contribution", "lower", "upper"))
            }
            if (available)
            else None,
        }
    output["native_diagnostics"] = {
        "status_counts": dict(
            Counter(
                (
                    row["status"]
                    for (model) in (frozen["models"])
                    for (row) in (model["cells"].values())
                )
            )
        ),
        "termination_counts": dict(
            Counter(
                (
                    row["termination_reason"]
                    for (model) in (frozen["models"])
                    for (row) in (model["cells"].values())
                )
            )
        ),
        "fixed_point_error_certified": all(
            (
                row["fixed_point_error_certified"]
                for (model) in (frozen["models"])
                for (row) in (model["cells"].values())
            )
        ),
        "interpretation": "Finite-iteration native implementation sensitivity; threshold stopping alone does not certify fixed-point error",
    }
    return output
