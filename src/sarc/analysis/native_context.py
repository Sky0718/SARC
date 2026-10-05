import copy
import math
from collections import Counter
from fractions import Fraction

SEEDS = [104729, 130363, 155921, 196613, 228017]
ARMS = ["N0", "N1", "C0", "C1"]
SIGNS = ["STRICTLY_POSITIVE", "STRICTLY_NEGATIVE", "IDENTIFIED_ZERO", "SIGN_UNRESOLVED"]
CONTRASTS = {
    "native_update": {"N1": 1, "N0": -1},
    "common_update": {"C1": 1, "C0": -1},
    "update_difference": {"N1": 1, "N0": -1, "C1": -1, "C0": 1},
}
NEW_CONTRASTS = ["native_update", "update_difference"]

def need(condition, message):
    if (not condition):
        raise ValueError(message)

def fraction(value):
    need(isinstance(value, str), "Rational must be a string")
    answer = Fraction(value)
    need(str(answer) == value, "Noncanonical rational")
    return answer

def terms(mapping):
    return [
        {"sample_id": sample, "gene": gene, "coefficient": str(value)}
        for (((sample, gene), value)) in (sorted(mapping.items()))
        if (value)
    ]

def unterms(records):
    result = {}
    for (row) in (records):
        key = (row["sample_id"], row["gene"])
        need(key not in result, "Repeated coefficient identity")
        result[key] = fraction(row["coefficient"])
        need(result[key] != 0, "Zero coefficient must cancel")
    return result

def combine(parts):
    output = {}
    for (scale, mapping) in (parts):
        for (key, value) in (mapping.items()):
            output[key] = output.get(key, Fraction()) + scale * value
    return {key: value for ((key, value)) in (output.items()) if (value)}

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

def unique_strings(values, name):
    need(
        isinstance(values, list)
        and all((isinstance(value, str) and value for (value) in (values))),
        name + " identifiers invalid",
    )
    need(len(values) == len(set(values)), name + " duplicate identifiers")
    return set(values)

def validate_protocol(protocol, data_class):
    need(data_class in ("PRODUCTION", "SYNTHETIC"), "Data class invalid")
    need(
        protocol["schema"] == "native_context_discriminator_protocol_v1",
        "Protocol schema differs",
    )
    need(
        protocol["budget"] == 10 and protocol["contrasts"] == CONTRASTS,
        "Frozen budget or contrasts differ",
    )
    need(protocol["population"]["seeds"] == SEEDS, "Seed alias scope differs")
    need(
        protocol["transition"] == ["native_11_0", "native_11_5"]
        and protocol["method"] == "DawnRank",
        "Method or transition differs",
    )
    if (data_class == "PRODUCTION"):
        need(
            {
                key: protocol["population"][key]
                for (key) in (("models", "donors", "normal_reference", "common_genes"))
            }
            == {
                "models": 85,
                "donors": 83,
                "normal_reference": 34,
                "common_genes": 7399,
            },
            "Production population reduced",
        )

def allocate(native, query):
    status = native["status"]
    need(
        status in ("SUCCESS", "SUCCESS_EMPTY", "FAILED", "UNAVAILABLE"),
        "Native state invalid",
    )
    raw = native["scores"]
    need(isinstance(raw, list), "Native scores not a list")
    if (status in ("FAILED", "UNAVAILABLE")):
        need(
            not raw and native.get("reason"),
            "Unavailable state requires empty ranks and reason",
        )
        return {
            "status": status,
            "reason": native["reason"],
            "weights": None,
            "capacity": None,
            "vacancies": None,
            "emitted_count": None,
            "eligible_emitted_count": None,
            "selected_gene_count": None,
            "seed_references": SEEDS[:],
            "seed_semantics": "Five deterministic aliases, not independent replicates",
        }
    need(
        (status == "SUCCESS") == bool(raw),
        "Native status and score availability disagree",
    )
    roots = unique_strings(native["native_root_genes"], "Native roots")
    scores = []
    genes = set()
    for (row) in (raw):
        need(
            isinstance(row, (list, tuple)) and len(row) == 3,
            "Native rank triple invalid",
        )
        gene, score_token, percentile_token = row
        need(
            isinstance(gene, str) and gene and (gene not in genes),
            "Repeated or invalid native gene",
        )
        need(
            isinstance(score_token, str) and isinstance(percentile_token, str),
            "Exact native score tokens required",
        )
        score, percentile = (float(score_token), float(percentile_token))
        need(
            math.isfinite(score)
            and score >= 0
            and math.isfinite(percentile)
            and (0 < percentile < 100),
            "Nonfinite or invalid native score",
        )
        scores.append((gene, score, percentile))
        genes.add(gene)
    need(
        genes == roots and query <= roots,
        "Full native roots or fixed query coverage differs",
    )
    for (left, right) in (zip(scores, scores[1:])):
        need(left[1] >= right[1] and left[2] >= right[2], "Native rank not descending")
        need(
            (left[1] == right[1]) == (left[2] == right[2]),
            "Score/percentile tie partitions differ",
        )
    selected = [row for (row) in (scores) if (row[0] in query)]
    weights, offset = ({}, 0)
    while (offset < len(selected) and offset < 10):
        end = offset + 1
        while (end < len(selected) and selected[end][2] == selected[offset][2]):
            end += 1
        amount = Fraction(min(end - offset, 10 - offset), end - offset)
        weights.update(((row[0], amount) for (row) in (selected[offset:end])))
        offset = end
    capacity = sum(weights.values(), Fraction())
    need(capacity == min(10, len(query)), "Native fixed-query allocation mass differs")
    return {
        "status": status,
        "reason": native.get("reason", ""),
        "weights": [[gene, str(value)] for ((gene, value)) in (sorted(weights.items()))],
        "capacity": str(capacity),
        "vacancies": str(10 - capacity),
        "emitted_count": len(raw),
        "eligible_emitted_count": len(query),
        "selected_gene_count": len(weights),
        "seed_references": SEEDS[:],
        "seed_semantics": "Five deterministic aliases, not independent replicates",
    }

def arm_map(arm, sample, query):
    status = arm["status"]
    need(
        status in ("SUCCESS", "SUCCESS_EMPTY", "FAILED", "UNAVAILABLE"),
        "Arm state invalid",
    )
    if (status in ("FAILED", "UNAVAILABLE")):
        need(
            arm["weights"] is None
            and arm["capacity"] is None
            and (arm["vacancies"] is None),
            "Unavailable arm carries invented mass",
        )
        return None
    result = {}
    for (gene, token) in (arm["weights"]):
        need(
            gene in query and (sample, gene) not in result,
            "Arm gene outside query or repeated",
        )
        value = fraction(token)
        need(0 < value <= 1, "Allocation weight outside unit capacity")
        result[sample, gene] = value
    mass = sum(result.values(), Fraction())
    need(
        mass == fraction(arm["capacity"]) == min(10, len(query)),
        "Common/native query allocation capacity differs",
    )
    need(fraction(arm["vacancies"]) == 10 - mass, "Vacancy accounting differs")
    need(arm["seed_references"] == SEEDS, "Arm seed aliases differ")
    return result

def prepare(protocol, context):
    data_class = context["data_class"]
    validate_protocol(protocol, data_class)
    need(
        context["schema"] == "native_context_discriminator_context_v1",
        "Context schema differs",
    )
    need(
        context["budget"] == 10 and context["master_seeds"] == SEEDS,
        "Context budget or aliases differ",
    )
    genes = unique_strings(context["common_genes"], "Common genes")
    population = protocol["population"]
    need(len(genes) == population["common_genes"], "Common-domain count differs")
    model_ids, sample_ids, donors, requests, output = (set(), set(), set(), set(), [])
    for (model) in (context["models"]):
        mid, sample, donor = (
            model[key] for (key) in (("model_id", "sample_id", "donor_id"))
        )
        need(
            all(
                (isinstance(value, str) and value for (value) in ((mid, sample, donor)))
            ),
            "Model identity invalid",
        )
        need(
            mid not in model_ids and sample not in sample_ids,
            "Repeated model/sample identity",
        )
        model_ids.add(mid)
        sample_ids.add(sample)
        donors.add(donor)
        roots = unique_strings(model["root_genes"], "Common mutation roots")
        query = unique_strings(model["eligible_genes"], "Eligible query")
        need(query <= roots <= genes, "Query/root/common-domain mismatch")
        need(
            set(model["native"]) == {"N0", "N1"}
            and set(model["common"]) == {"C0", "C1"},
            "Four-arm input differs",
        )
        arms = {
            name: allocate(model["native"][name], query) for (name) in (("N0", "N1"))
        }
        arms.update(
            {name: copy.deepcopy(model["common"][name]) for (name) in (("C0", "C1"))}
        )
        maps = {name: arm_map(arms[name], sample, query) for (name) in (ARMS)}
        for (mapping) in (maps.values()):
            if (mapping is not None):
                requests.update(mapping)
        contrasts = {}
        for (name, specification) in (CONTRASTS.items()):
            missing = [arm for (arm) in (specification) if (maps[arm] is None)]
            mapping = (
                {}
                if (missing)
                else combine(
                    [
                        (Fraction(scale), maps[arm])
                        for ((arm, scale)) in (specification.items())
                    ]
                )
            )
            contrasts[name] = coefficient_record(mapping, missing)
        output.append(
            {
                "model_id": mid,
                "sample_id": sample,
                "donor_id": donor,
                "arms": arms,
                "contrasts": contrasts,
            }
        )
    need(
        len(output) == population["models"] and len(donors) == population["donors"],
        "Full model/donor population differs",
    )
    counts = {
        "models": len(output),
        "donors": len(donors),
        "normal_reference": population["normal_reference"],
        "common_genes": len(genes),
        "reused_model_arm_states": 4 * len(output),
        "logical_seed_references": 4 * len(output) * len(SEEDS),
        "reused_native_model_states": 2 * len(output),
        "reused_common_diagonal_states": 2 * len(output),
        "original_native_cohort_jobs": 2,
        "new_native_fit_count": 0,
    }
    return {
        "schema": "native_context_discriminator_allocation_freeze_v1",
        "data_class": data_class,
        "prelabel_freeze": True,
        "counts": counts,
        "budget": 10,
        "master_seeds": SEEDS[:],
        "arms": ARMS[:],
        "contrast_definitions": copy.deepcopy(CONTRASTS),
        "models": output,
        "requested_labels": [
            {"sample_id": sample, "gene": gene} for ((sample, gene)) in (sorted(requests))
        ],
        "label_contract": copy.deepcopy(context["label_contract"]),
        "table5_source": copy.deepcopy(context["table5_source"]),
        "input_bindings": copy.deepcopy(context.get("input_bindings", {})),
        "exposure": "Previously exposed Table5 calls; exploratory local pipeline-context discriminator, not independent confirmation",
        "common_allocation_recomputed": False,
    }

def validate_call(row):
    key = (row["sample_id"], row["gene"])
    need(
        all((isinstance(value, str) and value for (value) in (key))),
        "Label identity invalid",
    )
    state, value = (row["state"], row["value"])
    need(state in ("MEASURED", "EXPLICIT_NULL", "ABSENT"), "Label state invalid")
    need(
        type(value) is int and value in (0, 1)
        if (state == "MEASURED")
        else value is None,
        "Label value semantics differ",
    )
    return (key, (state, value))

def merge_selected(frozen, documents):
    needed = {(row["sample_id"], row["gene"]) for (row) in (frozen["requested_labels"])}
    merged, origins, counts = ({}, {}, [])
    for (index, document) in (enumerate(documents)):
        need(
            document["label_contract"] == frozen["label_contract"]
            and document["data_class"] == frozen["data_class"],
            "Reused label contract differs",
        )
        need(
            document["source_binding"] == frozen["table5_source"],
            "Reused Table5 source differs",
        )
        seen, reused = (set(), 0)
        for (row) in (document["calls"]):
            key, value = validate_call(row)
            need(key not in seen, "Duplicate selected call within source")
            seen.add(key)
            if (key not in needed):
                continue
            if (key in merged):
                need(
                    validate_call(merged[key])[1] == value,
                    "Accepted selected labels conflict",
                )
            else:
                merged[key] = copy.deepcopy(row)
            origins.setdefault(key, []).append(index)
            reused += 1
        counts.append({"source_index": index, "matched_requested_cells": reused})
    return {
        "reused_calls": [merged[key] for (key) in (sorted(merged))],
        "missing_requests": [
            {"sample_id": sample, "gene": gene}
            for ((sample, gene)) in (sorted(needed - set(merged)))
        ],
        "reuse_counts": counts,
        "source_indices_by_cell": [
            {"sample_id": sample, "gene": gene, "source_indices": origins[sample, gene]}
            for ((sample, gene)) in (sorted(merged))
        ],
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
    known, low, high, masks = (Fraction(), Fraction(), Fraction(), Counter())
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
            low += min(coefficient, 0)
            high += max(coefficient, 0)
    low += known
    high += known
    return {
        "known_contribution": str(known),
        "lower": str(low),
        "upper": str(high),
        "per_offered_slot": {
            "known_contribution": str(known / 10),
            "lower": str(low / 10),
            "upper": str(high / 10),
        },
        "sign": sign(low, high),
        "masks": dict(sorted(masks.items())),
    }

def evaluated(record, calls):
    return {
        "available": record["available"],
        "value": bounds(unterms(record["coefficients"]), calls)
        if (record["available"])
        else None,
        **{
            key: copy.deepcopy(record[key])
            for (key) in (
                (
                    "positive_mass",
                    "negative_mass",
                    "half_l1",
                    "capacity_contrast",
                    "missing",
                )
            )
        },
    }

def aggregate(records, names):
    need(records, "Empty aggregation not permitted")
    output = {}
    for (name) in (names):
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
        output[name] = coefficient_record(mapping, missing)
    return output

def mass(arm, sample, calls):
    if (arm["weights"] is None):
        return None
    result = {
        name: Fraction()
        for (name) in (
            (
                "offered",
                "returned",
                "assessed",
                "known_positive",
                "known_zero",
                "explicit_null",
                "absent",
                "unknown",
                "vacant",
            )
        )
    }
    result["offered"] = 10
    for (gene, token) in (arm["weights"]):
        weight = fraction(token)
        state, value = calls[sample, gene]
        result["returned"] += weight
        if (state == "MEASURED"):
            result["assessed"] += weight
            result["known_positive" if (value) else "known_zero"] += weight
        else:
            result["unknown"] += weight
            result["explicit_null" if (state == "EXPLICIT_NULL") else "absent"] += (
                weight
            )
    result["vacant"] = 10 - result["returned"]
    return {key: str(value) for ((key, value)) in (result.items())}

def mean_mass(records):
    need(records, "Empty mass mean")
    if (any((row is None for (row) in (records)))):
        return None
    return {
        name: str(
            sum((fraction(row[name]) for (row) in (records)), Fraction()) / len(records)
        )
        for (name) in (records[0])
    }

def common_maps(reference, frozen):
    need(
        isinstance(reference["binding"], dict)
        and isinstance(reference["readout"], dict),
        "Accepted common reference required",
    )
    original = reference["readout"]
    need(
        original["schema"] == "dawn_factorial_readout_v1"
        and original["data_class"] == frozen["data_class"],
        "Common readout schema/class differs",
    )
    need(
        original["budget"] == 10 and original["master_seeds"] == SEEDS,
        "Common readout budget/aliases differ",
    )
    models, donors = ({}, {})
    for (row) in (original["models"]):
        need(row["model_id"] not in models, "Common duplicate model")
        models[row["model_id"]] = row
    for (row) in (original["donors"]):
        need(row["donor_id"] not in donors, "Common duplicate donor")
        donors[row["donor_id"]] = row
    need(
        set(models) == {row["model_id"] for (row) in (frozen["models"])}
        and set(donors) == {row["donor_id"] for (row) in (frozen["models"])},
        "Common reference population differs",
    )
    for (row) in (frozen["models"]):
        need(
            all(
                (
                    models[row["model_id"]][key] == row[key]
                    for (key) in (("sample_id", "donor_id"))
                )
            ),
            "Common model identity differs",
        )
    return (models, donors, original["population"]["contrasts"]["joint"])

def contrast_sign(record):
    return record["value"]["sign"] if (record["available"]) else "NATIVE_UNAVAILABLE"

def evaluate(protocol, frozen, selected, common_reference):
    validate_protocol(protocol, frozen["data_class"])
    need(
        frozen["schema"] == "native_context_discriminator_allocation_freeze_v1"
        and frozen["prelabel_freeze"] is True,
        "Prelabel freeze required",
    )
    need(
        frozen["budget"] == 10
        and frozen["master_seeds"] == SEEDS
        and (frozen["contrast_definitions"] == CONTRASTS),
        "Allocation contract differs",
    )
    need(
        selected["label_contract"] == frozen["label_contract"]
        and selected["data_class"] == frozen["data_class"],
        "Selected label contract differs",
    )
    need(
        selected["source_binding"] == frozen["table5_source"],
        "Selected label source differs",
    )
    calls = {}
    for (row) in (selected["calls"]):
        key, value = validate_call(row)
        need(key not in calls, "Duplicate selected call")
        calls[key] = value
    need(
        set(calls)
        == {(row["sample_id"], row["gene"]) for (row) in (frozen["requested_labels"])},
        "Exact selected support differs",
    )
    old_models, old_donors, old_population = common_maps(common_reference, frozen)
    output = {
        "schema": "native_context_discriminator_readout_v1",
        "data_class": frozen["data_class"],
        "counts": copy.deepcopy(frozen["counts"]),
        "budget": 10,
        "master_seeds": SEEDS[:],
        "models": [],
        "donors": [],
        "population": {},
        "donor_sign_counts": {},
        "leave_one_donor_out": [],
        "leaveout_summary": {},
        "selected_call_count": len(calls),
        "bounds_type": "Sharp shared-binary-unknown identification bounds, not confidence intervals, p-values or assay-error guarantees",
        "exposure": frozen["exposure"],
        "effect_unit": "Assay-positive allocated gene mass per ten offered slots; exact donor-weighted means",
        "scientifically_complete": all(
            (
                arm["status"] in ("SUCCESS", "SUCCESS_EMPTY")
                for (model) in (frozen["models"])
                for (arm) in (model["arms"].values())
            )
        ),
        "allocation_status_counts": dict(
            sorted(
                Counter(
                    (
                        arm["status"]
                        for (model) in (frozen["models"])
                        for (arm) in (model["arms"].values())
                    )
                ).items()
            )
        ),
    }
    donor_models, by_model = ({}, {})
    for (row) in (frozen["models"]):
        donor_models.setdefault(row["donor_id"], []).append(row)
        contrasts = {
            name: evaluated(row["contrasts"][name], calls) for (name) in (NEW_CONTRASTS)
        }
        contrasts["common_update"] = copy.deepcopy(
            old_models[row["model_id"]]["contrasts"]["joint"]
        )
        record = {
            key: row[key] for (key) in (("model_id", "sample_id", "donor_id"))
        } | {
            "arms": {
                name: {
                    "status": row["arms"][name]["status"],
                    "mass": mass(row["arms"][name], row["sample_id"], calls),
                }
                for (name) in (ARMS)
            },
            "contrasts": contrasts,
        }
        output["models"].append(record)
        by_model[row["model_id"]] = record
    donor_freezes = []
    for (donor, models) in (sorted(donor_models.items())):
        new_records = aggregate(models, NEW_CONTRASTS)
        donor_freezes.append({"donor_id": donor, "contrasts": new_records})
        contrasts = {
            name: evaluated(new_records[name], calls) for (name) in (NEW_CONTRASTS)
        }
        contrasts["common_update"] = copy.deepcopy(
            old_donors[donor]["contrasts"]["joint"]
        )
        need(
            set(old_donors[donor]["model_ids"])
            == {row["model_id"] for (row) in (models)},
            "Common donor membership differs",
        )
        output["donors"].append(
            {
                "donor_id": donor,
                "model_ids": [row["model_id"] for (row) in (models)],
                "arms": {
                    name: {
                        "mass": mean_mass(
                            [
                                by_model[row["model_id"]]["arms"][name]["mass"]
                                for (row) in (models)
                            ]
                        )
                    }
                    for (name) in (ARMS)
                },
                "contrasts": contrasts,
            }
        )
    total = aggregate(donor_freezes, NEW_CONTRASTS)
    contrasts = {name: evaluated(total[name], calls) for (name) in (NEW_CONTRASTS)}
    contrasts["common_update"] = copy.deepcopy(old_population)
    output["population"] = {
        "donor_ids": [row["donor_id"] for (row) in (donor_freezes)],
        "arms": {
            name: {
                "mass": mean_mass(
                    [row["arms"][name]["mass"] for (row) in (output["donors"])]
                )
            }
            for (name) in (ARMS)
        },
        "contrasts": contrasts,
    }
    for (name) in (CONTRASTS):
        output["donor_sign_counts"][name] = dict.fromkeys(
            SIGNS + ["NATIVE_UNAVAILABLE"], 0
        )
        for (row) in (output["donors"]):
            output["donor_sign_counts"][name][
                contrast_sign(row["contrasts"][name])
            ] += 1
    matrix = {left: {right: 0 for (right) in (SIGNS)} for (left) in (SIGNS)}
    identities = {left: {right: [] for (right) in (SIGNS)} for (left) in (SIGNS)}
    unavailable = []
    for (row) in (output["donors"]):
        left, right = (
            contrast_sign(row["contrasts"]["native_update"]),
            contrast_sign(row["contrasts"]["common_update"]),
        )
        if ("NATIVE_UNAVAILABLE" in (left, right)):
            unavailable.append(row["donor_id"])
        else:
            matrix[left][right] += 1
            identities[left][right].append(row["donor_id"])
    output["donor_direction_table"] = {
        "row_axis": "native_update",
        "column_axis": "common_update",
        "signs": SIGNS[:],
        "counts": matrix,
        "donor_ids": identities,
        "available_donors": len(output["donors"]) - len(unavailable),
        "unavailable_donors": len(unavailable),
        "unavailable_donor_ids": unavailable,
    }
    if (len(donor_freezes) > 1):
        for (removed) in (donor_freezes):
            remaining = [
                row
                for (row) in (donor_freezes)
                if (row["donor_id"] != removed["donor_id"])
            ]
            local = aggregate(remaining, NEW_CONTRASTS)
            output["leave_one_donor_out"].append(
                {
                    "removed_donor_id": removed["donor_id"],
                    "remaining_donors": len(remaining),
                    "contrasts": {
                        name: evaluated(local[name], calls)
                        for (name) in (NEW_CONTRASTS)
                    },
                }
            )
    for (name) in (NEW_CONTRASTS):
        values = [row["contrasts"][name] for (row) in (output["leave_one_donor_out"])]
        available = [row["value"] for (row) in (values) if (row["available"])]
        output["leaveout_summary"][name] = {
            "available": len(available),
            "unavailable": len(values) - len(available),
            "sign_counts": dict(
                sorted(Counter((row["sign"] for (row) in (available))).items())
            ),
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
    output["common_reuse"] = {
        "source_binding": copy.deepcopy(common_reference["binding"]),
        "acceptance": copy.deepcopy(common_reference.get("acceptance")),
        "contrast_field": "joint",
        "leaveout_field": "leave_one_donor_out[*].contrasts.joint",
        "leaveout_recomputed": False,
        "source_schema": common_reference["readout"]["schema"],
    }
    native_sign, common_sign = (
        contrast_sign(output["population"]["contrasts"][name])
        for (name) in (("native_update", "common_update"))
    )
    opposite = {native_sign, common_sign} == {"STRICTLY_POSITIVE", "STRICTLY_NEGATIVE"}
    state = (
        "UNAVAILABLE"
        if ("NATIVE_UNAVAILABLE" in (native_sign, common_sign))
        else "OPPOSITE_IDENTIFIED_SIGNS"
        if (opposite)
        else "NO_OPPOSITE_IDENTIFIED_SIGNS"
    )
    output["primary_decision"] = {
        "native_sign": native_sign,
        "common_sign": common_sign,
        "opposite_identified_signs": opposite,
        "status": state,
        "qualification": "Fixed-cohort zero-switching-cost directional comparison; not future or clinical adoption guarantee. Same direction or unresolved is not equivalence.",
    }
    return output
