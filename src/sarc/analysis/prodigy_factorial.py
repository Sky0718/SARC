import math
import re
import struct
from collections import Counter
from fractions import Fraction

SEEDS = [104729, 130363, 155921, 196613, 228017]
GRAPHS = ["11.0", "11.5"]
EVIDENCE = ["H0", "H1"]
LABEL_CONTRACT = "final_Table5_binary_original_mapping"
CONTRASTS = {
    "graph_H0": {"11.5/H0": 1, "11.0/H0": -1},
    "graph_H1": {"11.5/H1": 1, "11.0/H1": -1},
    "evidence_G0": {"11.0/H1": 1, "11.0/H0": -1},
    "evidence_G1": {"11.5/H1": 1, "11.5/H0": -1},
    "interaction": {"11.5/H1": 1, "11.0/H1": -1, "11.5/H0": -1, "11.0/H0": 1},
}

def require(condition, message):
    if (not condition):
        raise ValueError(message)

def fraction(value):
    require(isinstance(value, str), "Fraction must be a string")
    result = Fraction(value)
    require(str(result) == value, "Noncanonical fraction")
    return result

def identities(values, name):
    require(
        isinstance(values, list)
        and all((isinstance(x, str) and x for (x) in (values))),
        name,
    )
    require(len(values) == len(set(values)), "Duplicate " + name)

def source_identity(value):
    require(
        isinstance(value, str) and bool(value.strip()),
        "Explicit scientific source identity required",
    )

def binding(value):
    require(isinstance(value, dict), "Source record required")
    source_identity(value.get("source_id"))

def validate_scope(scope):
    require(
        scope.get("schema") == "common_domain_factorial_analysis_scope_v1",
        "Scope schema",
    )
    require(scope.get("data_class") in ("SYNTHETIC", "PRODUCTION"), "Data class")
    require(scope.get("exploratory") is True, "Exploratory designation required")
    require(
        scope.get("transition") == GRAPHS and scope.get("graph_states") == GRAPHS,
        "Transition",
    )
    require(
        scope.get("evidence_states") == EVIDENCE and scope.get("master_seeds") == SEEDS,
        "Fixed states",
    )
    require(
        type(scope.get("budget")) is int and scope["budget"] == 10, "Fixed budget10"
    )
    require(scope.get("label_contract") == LABEL_CONTRACT, "Label contract")
    require(
        isinstance(scope.get("source_contracts"), list) and scope["source_contracts"],
        "Source references",
    )
    models = scope.get("models")
    require(isinstance(models, list) and models, "Model roster")
    identities([m["model_id"] for (m) in (models)], "models")
    identities([m["sample_id"] for (m) in (models)], "samples")
    for (model) in (models):
        require(
            all(
                (
                    isinstance(model.get(k), str) and model[k]
                    for (k) in (("donor_id", "sample_id"))
                )
            ),
            "Model IDs",
        )
        identities(model["root_genes"], "root genes")
        identities(model["eligible_genes"], "eligible genes")
        require(
            set(model["eligible_genes"]) <= set(model["root_genes"]),
            "Eligibility outside roots",
        )
        source_identity(model["root_axis_id"])
        require(set(model["evidence_id"]) == set(EVIDENCE), "Evidence binding states")
        for (value) in (model["evidence_id"].values()):
            source_identity(value)
        if (scope["data_class"] == "SYNTHETIC"):
            require(
                all(
                    (
                        model[k].startswith("SYN_")
                        for (k) in (("model_id", "donor_id", "sample_id"))
                    )
                ),
                "Fixture IDs",
            )
    donors = len({m["donor_id"] for (m) in (models)})
    require(
        scope.get("expected")
        == {
            "models": len(models),
            "donors": donors,
            "cells": len(models) * 4,
            "seed_states": len(models) * 20,
        },
        "Expected roster counts",
    )
    if (scope["data_class"] == "PRODUCTION"):
        require(
            (len(models), donors) == (85, 83), "Production must retain85models/83donors"
        )
    source_identity(scope["common_domain_id"])
    source_identity(scope["pathway_dictionary_id"])
    require(set(scope["graph_id"]) == set(GRAPHS), "Graph binding states")
    for (value) in (scope["graph_id"].values()):
        source_identity(value)

def allocate(state, model):
    status = state.get("status")
    require(
        status in ("SUCCESS", "SUCCESS_EMPTY", "UNAVAILABLE", "FAILED"),
        "Explicit native status",
    )
    binding(state.get("native_binding"))
    genes, codes = (state.get("genes"), state.get("scores_binary64_le"))
    if (status in ("UNAVAILABLE", "FAILED")):
        require(
            genes is None
            and codes is None
            and isinstance(state.get("reason"), str)
            and state["reason"],
            "Unavailable state must preserve reason and null rank",
        )
        return None
    identities(genes, "emitted genes")
    require(isinstance(codes, list) and len(codes) == len(genes), "Score length")
    require(
        set(genes) <= set(model["root_genes"]),
        "Emitted gene outside complete root axis",
    )
    if (status == "SUCCESS_EMPTY"):
        require(
            not genes and isinstance(state.get("reason"), str) and state["reason"],
            "Empty reason/arrays",
        )
    else:
        require(bool(genes), "Nonempty success expected")
    scores = []
    for (code) in (codes):
        require(
            isinstance(code, str) and re.fullmatch("[0-9a-f]{16}", code),
            "Binary64 encoding",
        )
        value = struct.unpack("<d", bytes.fromhex(code))[0]
        require(
            math.isfinite(value) and value > 0, "Native scores must be finite positive"
        )
        scores.append(value)
    require(
        all((a >= b for ((a, b)) in (zip(scores, scores[1:])))), "Ranking not descending"
    )
    eligible = set(model["eligible_genes"])
    rows = [(g, score) for ((g, score)) in (zip(genes, scores)) if (g in eligible)]
    weights = {}
    offset = 0
    while (offset < len(rows)):
        end = offset + 1
        while (end < len(rows) and rows[end][1] == rows[offset][1]):
            end += 1
        amount = min(end - offset, max(0, 10 - offset))
        for (gene, _) in (rows[offset:end]):
            weights[gene] = Fraction(amount, end - offset)
        offset = end
    require(
        sum(weights.values(), Fraction()) == min(10, len(rows)), "Allocation capacity"
    )
    return weights

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
        require(key not in result, "Duplicate coefficient variable")
        result[key] = fraction(row["coefficient"])
        require(result[key] != 0, "Zero coefficient should be removed")
    return result

def combine(parts):
    result = {}
    for (scale, part) in (parts):
        for (key, value) in (part.items()):
            result[key] = result.get(key, Fraction()) + scale * value
    return {key: value for ((key, value)) in (result.items()) if (value)}

def coefficient_record(name, available, coeff, missing = None):
    if (not available):
        return {
            "contrast": name,
            "available": False,
            "coefficients": None,
            "missing": missing,
        }
    positive = sum((c for (c) in (coeff.values()) if (c > 0)), Fraction())
    negative = -sum((c for (c) in (coeff.values()) if (c < 0)), Fraction())
    return {
        "contrast": name,
        "available": True,
        "coefficients": terms(coeff),
        "positive_mass": str(positive),
        "negative_mass": str(negative),
        "l1": str(positive + negative),
        "half_l1": str((positive + negative) / 2),
        "capacity_contrast": str(positive - negative),
        "all_labels_lower": str(-negative),
        "all_labels_upper": str(positive),
        "algebraic_zero": not coeff,
    }

def freeze(scope, carriers):
    validate_scope(scope)
    require(
        carriers.get("schema") == "common_domain_factorial_native_carriers_v1",
        "Carrier schema",
    )
    require(carriers.get("data_class") == scope["data_class"], "Carrier data class")
    require(carriers.get("declared_complete") is True, "Incomplete carrier declaration")
    cells = carriers.get("cells")
    require(
        isinstance(cells, list) and len(cells) == scope["expected"]["cells"],
        "Full cell count required",
    )
    models = {m["model_id"]: m for (m) in (scope["models"])}
    lookup, signatures, allocations = ({}, set(), [])
    for (cell) in (cells):
        mid, graph, evidence = (
            cell.get("model_id"),
            cell.get("graph"),
            cell.get("evidence"),
        )
        require(
            mid in models and graph in GRAPHS and (evidence in EVIDENCE),
            "Unplanned cell",
        )
        key = (mid, graph, evidence)
        require(key not in lookup, "Duplicate cell")
        model = models[mid]
        signature = cell.get("signature_id")
        require(
            isinstance(signature, str) and signature and (signature not in signatures),
            "Signature identity",
        )
        signatures.add(signature)
        for (field) in (("sample_id", "donor_id", "root_axis_id")):
            require(cell.get(field) == model[field], "Cell model binding: " + field)
        for (field) in (("common_domain_id", "pathway_dictionary_id")):
            require(cell.get(field) == scope[field], "Common input binding: " + field)
        require(cell.get("graph_id") == scope["graph_id"][graph], "Graph binding")
        require(
            cell.get("evidence_id") == model["evidence_id"][evidence],
            "Evidence binding",
        )
        states = cell.get("seed_states")
        require(
            isinstance(states, list) and len(states) == 5,
            "Five explicit seed states required",
        )
        require(
            all((type(s.get("master_seed")) is int for (s) in (states))),
            "Integer seed identities",
        )
        require(
            sorted((s["master_seed"] for (s) in (states))) == sorted(SEEDS),
            "Seed identities/duplicates",
        )
        lookup[key] = {}
        for (state) in (states):
            require(
                state.get("evidence_id") == model["evidence_id"][evidence],
                "H_NOT_SEED_INVARIANT_OR_UNBOUND",
            )
            allocation = allocate(state, model)
            lookup[key][state["master_seed"]] = allocation
            allocations.append(
                {
                    "model_id": mid,
                    "graph": graph,
                    "evidence": evidence,
                    "master_seed": state["master_seed"],
                    "status": state["status"],
                    "reason": state.get("reason"),
                    "weights": None
                    if (allocation is None)
                    else [[g, str(w)] for ((g, w)) in (allocation.items())],
                    "emitted_count": None
                    if (state["genes"] is None)
                    else len(state["genes"]),
                    "eligible_emitted_count": None
                    if (allocation is None)
                    else len(allocation),
                    "capacity": None
                    if (allocation is None)
                    else str(sum(allocation.values(), Fraction())),
                    "eligible_unemitted_count": None
                    if (allocation is None)
                    else len(model["eligible_genes"]) - len(allocation),
                }
            )
    seed_records, model_records, summaries = ([], [], [])
    seed_lookup, model_lookup = ({}, {})
    for (model) in (scope["models"]):
        mid, sample = (model["model_id"], model["sample_id"])
        for (seed) in (SEEDS):
            for (name, signs) in (CONTRASTS.items()):
                selected = [
                    (sign, cell, lookup[(mid,) + tuple(cell.split("/"))][seed])
                    for ((cell, sign)) in (signs.items())
                ]
                missing = [cell for ((_, cell, value)) in (selected) if (value is None)]
                coeff = (
                    {}
                    if (missing)
                    else combine(
                        [
                            (
                                Fraction(sign),
                                {(sample, g): w for ((g, w)) in (value.items())},
                            )
                            for ((sign, _, value)) in (selected)
                        ]
                    )
                )
                row = coefficient_record(name, not missing, coeff, missing)
                row.update(model_id = mid, master_seed = seed)
                seed_records.append(row)
                seed_lookup[mid, seed, name] = row
        for (name) in (CONTRASTS):
            rows = [seed_lookup[mid, seed, name] for (seed) in (SEEDS)]
            available = all((row["available"] for (row) in (rows)))
            coeff = (
                combine(
                    [(Fraction(1, 5), unterms(row["coefficients"])) for (row) in (rows)]
                )
                if (available)
                else {}
            )
            row = coefficient_record(
                name,
                available,
                coeff,
                [r["master_seed"] for (r) in (rows) if (not r["available"])],
            )
            row.update(
                model_id = mid,
                donor_id = model["donor_id"],
                sample_id = sample,
                all_seed_contrasts_algebraically_zero = all(
                    (r.get("algebraic_zero", False) for (r) in (rows))
                ),
            )
            model_records.append(row)
            model_lookup[mid, name] = row
    donor_counts = Counter((m["donor_id"] for (m) in (scope["models"])))
    donor_total = len(donor_counts)
    weights = {
        m["model_id"]: Fraction(1, donor_total * donor_counts[m["donor_id"]])
        for (m) in (scope["models"])
    }
    for (name) in (CONTRASTS):
        available = [
            model_lookup[m["model_id"], name]
            for (m) in (scope["models"])
            if (model_lookup[m["model_id"], name]["available"])
        ]
        missing = [
            m["model_id"]
            for (m) in (scope["models"])
            if (not model_lookup[m["model_id"], name]["available"])
        ]
        coeff = combine(
            [
                (weights[r["model_id"]], unterms(r["coefficients"]))
                for (r) in (available)
            ]
        )
        row = coefficient_record(name, not missing, coeff, missing)
        row.update(
            available_original_weight_coefficients = terms(coeff),
            unavailable_models = missing,
            unavailable_weight = str(
                sum((weights[mid] for (mid) in (missing)), Fraction())
            ),
            available_models = len(available),
            total_models = len(models),
            total_donors = donor_total,
        )
        summaries.append(row)
    requests = sorted(
        {
            (r["sample_id"], r["gene"])
            for (row) in (seed_records)
            if (row["available"])
            for (r) in (row["coefficients"])
        }
    )
    return {
        "schema": "common_domain_factorial_coefficient_freeze_v1",
        "data_class": scope["data_class"],
        "label_contract": LABEL_CONTRACT,
        "scope": scope,
        "accounting_complete": True,
        "scientifically_complete": all((x["available"] for (x) in (seed_records))),
        "scientific_completeness_meaning": "Full native availability and arithmetic; unknown-label bounds may remain unresolved",
        "counts": dict(scope["expected"]),
        "allocations": allocations,
        "seed_contrasts": seed_records,
        "model_contrasts": model_records,
        "population_contrasts": summaries,
        "model_weights": {mid: str(w) for ((mid, w)) in (weights.items())},
        "requested_labels": [{"sample_id": s, "gene": g} for ((s, g)) in (requests)],
        "prelabel_freeze": True,
    }

def bounds(coefficients, calls):
    known, lower, upper = (Fraction(), Fraction(), Fraction())
    masks = Counter()
    for (key, coefficient) in (coefficients.items()):
        state, value = calls[key]
        masks[state] += 1
        if (state == "MEASURED"):
            known += coefficient * value
        else:
            lower += min(coefficient, 0)
            upper += max(coefficient, 0)
    low, high = (known + lower, known + upper)
    sign = (
        "STRICTLY_POSITIVE"
        if (low > 0)
        else "STRICTLY_NEGATIVE"
        if (high < 0)
        else "IDENTIFIED_ZERO"
        if (low == high == 0)
        else "SIGN_UNRESOLVED"
    )
    return {
        "known_contribution": str(known),
        "lower": str(low),
        "upper": str(high),
        "per_offered_slot": {
            "known_contribution": str(known / 10),
            "lower": str(low / 10),
            "upper": str(high / 10),
        },
        "masks": dict(sorted(masks.items())),
        "sign": sign,
    }

def evaluate(frozen, labels = None):
    require(
        frozen.get("schema") == "common_domain_factorial_coefficient_freeze_v1",
        "Freeze schema",
    )
    require(
        frozen.get("prelabel_freeze") is True
        and frozen.get("accounting_complete") is True,
        "Freeze incomplete",
    )
    requests = {(r["sample_id"], r["gene"]) for (r) in (frozen["requested_labels"])}
    require(
        len(requests) == len(frozen["requested_labels"]), "Duplicate frozen requests"
    )
    calls = {}
    if (requests):
        require(
            labels is not None
            and labels.get("schema") == "common_domain_factorial_selected_labels_v1",
            "Labels schema",
        )
        require(labels.get("data_class") == frozen["data_class"], "Label data class")
        require(labels.get("label_contract") == LABEL_CONTRACT, "Label semantics")
        binding(labels.get("source_binding"))
        require(
            isinstance(labels.get("exposure"), str) and labels["exposure"],
            "Exposure statement required",
        )
        for (row) in (labels.get("calls", [])):
            key = (row["sample_id"], row["gene"])
            require(key not in calls and key in requests, "Duplicate/unrequested call")
            state, value = (row["state"], row["value"])
            require(state in ("MEASURED", "EXPLICIT_NULL", "ABSENT"), "Call mask")
            if (state == "MEASURED"):
                require(type(value) is int and value in (0, 1), "Binary call required")
            else:
                require(value is None, "Unknown must have null value")
            calls[key] = (state, value)
        require(set(calls) == requests, "Missing selected calls")
    else:
        require(labels is None, "Do not open labels for algebraic zero")
    output = {
        "schema": "common_domain_factorial_readout_v1",
        "data_class": frozen["data_class"],
        "accounting_complete": True,
        "scientifically_complete": frozen["scientifically_complete"],
        "counts": frozen["counts"],
        "selected_call_count": len(calls),
        "unknown_identity": "literal(sample_id,gene), shared across cells and seeds",
        "bounds_type": "sharp marginal free-binary-completion bounds, not confidence intervals",
    }
    for (section) in (("seed_contrasts", "model_contrasts", "population_contrasts")):
        output[section] = []
        for (row) in (frozen[section]):
            identifiers = {
                key: row[key]
                for (key) in (
                    ("model_id", "donor_id", "sample_id", "master_seed", "contrast")
                )
                if (key in row)
            }
            result = dict(identifiers, available = row["available"], value = None)
            if (row["available"]):
                result["value"] = bounds(unterms(row["coefficients"]), calls)
            if (section == "population_contrasts"):
                result.update(
                    available_original_weight_contribution = bounds(
                        unterms(row["available_original_weight_coefficients"]), calls
                    ),
                    unavailable_models = row["unavailable_models"],
                    unavailable_weight = row["unavailable_weight"],
                )
            output[section].append(result)
    return output
