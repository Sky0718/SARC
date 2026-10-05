from fractions import Fraction

SEEDS = [104729, 130363, 155921, 196613, 228017]
PARTITIONS = ["intersection", "old_only", "new_only"]
SIGNS = ["STRICTLY_POSITIVE", "STRICTLY_NEGATIVE", "IDENTIFIED_ZERO", "SIGN_UNRESOLVED"]
CELL_NAMES = {"old": "G1_H1", "new": "G0_H0"}
LABEL_CONTRACT = "final_Table5_binary_original_mapping"

def need(value, message):
    if (not value):
        raise ValueError(message)

def rational(value):
    need(isinstance(value, str), "Canonical fraction string required")
    answer = Fraction(value)
    need(str(answer) == value, "Noncanonical fraction")
    return answer

def unique_strings(values, name):
    need(
        isinstance(values, list)
        and all((isinstance(value, str) and value for (value) in (values))),
        name + " must be literal nonempty strings",
    )
    need(len(values) == len(set(values)), name + " contains duplicates")

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

def partition(gene, domains):
    values = [name for (name) in (PARTITIONS) if (gene in domains[name])]
    need(
        len(values) == 1,
        "Selected gene outside or duplicated between domain partitions",
    )
    return values[0]

def allocation(cell, domain):
    need(
        cell["status"] in ("SUCCESS", "SUCCESS_EMPTY"),
        "Unsuccessful selected cell must remain a reported failure, not omitted",
    )
    need(
        cell["seed_references"] == SEEDS, "Complete deterministic seed labels required"
    )
    weights = {}
    for (row) in (cell["weights"]):
        need(
            isinstance(row, list) and len(row) == 2,
            "Exact gene/fraction allocation row required",
        )
        gene, text = row
        need(
            isinstance(gene, str) and gene in domain and (gene not in weights),
            "Duplicate or out-of-domain allocated gene",
        )
        weight = rational(text)
        need(0 < weight <= 1, "Positive allocation at most one per gene required")
        weights[gene] = weight
    capacity = sum(weights.values(), Fraction())
    need(
        0 <= capacity <= 10
        and capacity == rational(cell["capacity"])
        and (rational(cell["vacancies"]) == 10 - capacity),
        "Preserved capacity or vacancies differ",
    )
    need(
        (cell["status"] == "SUCCESS_EMPTY") == (cell["emitted_count"] == 0),
        "Native empty-state accounting differs",
    )
    need(
        type(cell["emitted_count"]) is int
        and type(cell["eligible_emitted_count"]) is int
        and (0 <= cell["eligible_emitted_count"] <= cell["emitted_count"]),
        "Native eligible output accounting differs",
    )
    need(
        cell["selected_gene_count"] == len(weights)
        and len(weights) <= cell["eligible_emitted_count"]
        and (capacity == min(10, cell["eligible_emitted_count"])),
        "Exact original cutoff allocation count differs",
    )
    return (
        {
            "status": cell["status"],
            "weights": [
                [gene, str(weight)] for ((gene, weight)) in (sorted(weights.items()))
            ],
            "capacity": str(capacity),
            "vacancies": str(10 - capacity),
        },
        weights,
    )

def prepare_bridge(protocol, old_freeze, new_freeze, old_domain, new_domain):
    need(
        protocol["schema"] == "dawn_middle_version_domain_bridge_v1"
        and protocol["cells"] == CELL_NAMES,
        "Frozen middle-version bridge protocol required",
    )
    scope = protocol["scope"]
    need(
        scope["budget"] == 10
        and scope["seeds"] == SEEDS
        and (scope["version"] == "native_11_5"),
        "Frozen query/version/seed contract differs",
    )
    unique_strings(old_domain, "Old domain")
    unique_strings(new_domain, "New domain")
    need(
        len(old_domain) == scope["old_domain_size"]
        and len(new_domain) == scope["new_domain_size"],
        "Complete pair-specific domains required",
    )
    old_set, new_set = (set(old_domain), set(new_domain))
    domains = {
        "old": sorted(old_set),
        "new": sorted(new_set),
        "intersection": sorted(old_set & new_set),
        "old_only": sorted(old_set - new_set),
        "new_only": sorted(new_set - old_set),
    }
    sets = {name: set(values) for ((name, values)) in (domains.items())}
    need(
        old_freeze["data_class"] == new_freeze["data_class"]
        and old_freeze["data_class"] in ("PRODUCTION", "SYNTHETIC"),
        "Matching data class required",
    )
    data_class = old_freeze["data_class"]
    for (document) in ((old_freeze, new_freeze)):
        need(
            document["schema"] == "dawn_factorial_allocation_freeze_v1"
            and document["prelabel_freeze"] is True,
            "Accepted allocation freeze required",
        )
        need(
            document["budget"] == 10
            and document["master_seeds"] == SEEDS
            and (document["label_contract"] == LABEL_CONTRACT),
            "Accepted budget/seed/label contract differs",
        )
        need(
            document["counts"]
            == {
                "models": scope["models"],
                "donors": scope["donors"],
                "physical_model_cells": 4 * scope["models"],
                "logical_seed_references": 20 * scope["models"],
            },
            "Accepted full census differs",
        )
        need(
            len(document["models"]) == scope["models"],
            "Complete model population required",
        )
    need(
        old_freeze["table5_source"] == new_freeze["table5_source"],
        "Same original assay source identity required",
    )
    if (data_class == "PRODUCTION"):
        need(
            (
                scope["models"],
                scope["donors"],
                scope["old_domain_size"],
                scope["new_domain_size"],
            )
            == (85, 83, 7399, 8620),
            "Complete frozen production scope required",
        )
    old_models, new_models = ({}, {})
    for (document, index) in (((old_freeze, old_models), (new_freeze, new_models))):
        samples = set()
        for (model) in (document["models"]):
            identity = tuple(
                (model[key] for (key) in (("model_id", "sample_id", "donor_id")))
            )
            need(
                all((isinstance(value, str) and value for (value) in (identity))),
                "Literal model/sample/donor identities required",
            )
            if (data_class == "SYNTHETIC"):
                need(
                    all((value.startswith("SYN_") for (value) in (identity))),
                    "Invented fixture identities required",
                )
            need(
                identity[0] not in index and identity[1] not in samples,
                "Duplicate model or sample",
            )
            index[identity[0]] = model
            samples.add(identity[1])
        need(
            len({model["donor_id"] for (model) in (document["models"])})
            == scope["donors"],
            "Complete donor population required",
        )
    need(set(old_models) == set(new_models), "Paired model coverage differs")
    requests, models = (set(), [])
    for (old) in (old_freeze["models"]):
        new = new_models[old["model_id"]]
        need(
            all((old[key] == new[key] for (key) in (("sample_id", "donor_id")))),
            "Paired sample or donor identity differs",
        )
        old_cell, old_weights = allocation(old["cells"][CELL_NAMES["old"]], old_set)
        new_cell, new_weights = allocation(new["cells"][CELL_NAMES["new"]], new_set)
        ledger = []
        for (gene) in (sorted(set(old_weights) | set(new_weights))):
            first, second = (
                old_weights.get(gene, Fraction()),
                new_weights.get(gene, Fraction()),
            )
            ledger.append(
                {
                    "gene": gene,
                    "partition": partition(gene, sets),
                    "old_weight": str(first),
                    "new_weight": str(second),
                    "coefficient": str(second - first),
                }
            )
            requests.add((old["sample_id"], gene))
        models.append(
            {key: old[key] for (key) in (("model_id", "sample_id", "donor_id"))}
            | {"old": old_cell, "new": new_cell, "ledger": ledger}
        )
    return {
        "schema": "dawn_middle_version_allocation_freeze_v1",
        "data_class": data_class,
        "prelabel_freeze": True,
        "budget": 10,
        "master_seeds": SEEDS,
        "label_contract": LABEL_CONTRACT,
        "table5_source": old_freeze["table5_source"],
        "counts": {
            "models": len(models),
            "donors": scope["donors"],
            "selected_model_cells": 2 * len(models),
            "logical_seed_references": 10 * len(models),
            "requested_sample_genes": len(requests),
        },
        "domains": domains,
        "models": models,
        "requested_labels": [
            {"sample_id": sample, "gene": gene} for ((sample, gene)) in (sorted(requests))
        ],
    }

def cache_calls(document, frozen):
    need(
        document["schema"] == "common_domain_factorial_selected_labels_v1"
        and document["data_class"] == frozen["data_class"],
        "Accepted selected-label cache required",
    )
    need(
        document["label_contract"] == LABEL_CONTRACT
        and document["source_binding"] == frozen["table5_source"],
        "Original assay source and contract differ",
    )
    output = {}
    for (row) in (document["calls"]):
        sample, gene = (row["sample_id"], row["gene"])
        need(
            isinstance(sample, str) and sample and isinstance(gene, str) and gene,
            "Literal label identity required",
        )
        key = (sample, gene)
        need(key not in output, "Duplicate selected call identity")
        state, value = (row["state"], row["value"])
        need(
            state in ("MEASURED", "EXPLICIT_NULL", "ABSENT"),
            "Unknown selected call state",
        )
        need(
            type(value) is int and value in (0, 1)
            if (state == "MEASURED")
            else value is None,
            "Preserve exact binary versus unknown semantics",
        )
        output[key] = (state, value)
    return output

def bound_stats(coefficients, calls):
    coefficients = {key: value for ((key, value)) in (coefficients.items()) if (value)}
    known, lower, upper = (Fraction(), Fraction(), Fraction())
    for (key, coefficient) in (coefficients.items()):
        state, value = calls[key]
        if (state == "MEASURED"):
            known += coefficient * value
        else:
            lower += min(Fraction(), coefficient)
            upper += max(Fraction(), coefficient)
    lower, upper = (known + lower, known + upper)
    positive = sum(
        (value for (value) in (coefficients.values()) if (value > 0)), Fraction()
    )
    negative = -sum(
        (value for (value) in (coefficients.values()) if (value < 0)), Fraction()
    )
    return {
        "known_contribution": str(known),
        "lower": str(lower),
        "upper": str(upper),
        "sign": sign(lower, upper),
        "positive_mass": str(positive),
        "negative_mass": str(negative),
        "half_l1": str((positive + negative) / 2),
        "capacity_contrast": str(positive - negative),
        "per_offered_slot": {
            "known_contribution": str(known / 10),
            "lower": str(lower / 10),
            "upper": str(upper / 10),
        },
    }

def cell_mass(weights, calls):
    mass = {
        name: Fraction()
        for (name) in (
            (
                "returned",
                "measured_positive",
                "measured_negative",
                "explicit_null",
                "absent",
            )
        )
    }
    for (key, weight) in (weights.items()):
        state, value = calls[key]
        mass["returned"] += weight
        field = (
            ("measured_positive" if (value) else "measured_negative")
            if (state == "MEASURED")
            else "explicit_null"
            if (state == "EXPLICIT_NULL")
            else "absent"
        )
        mass[field] += weight
    mass["assessed"] = mass["measured_positive"] + mass["measured_negative"]
    mass["unknown"] = mass["explicit_null"] + mass["absent"]
    mass["vacancies"] = 10 - mass["returned"]
    lower, upper = (
        mass["measured_positive"],
        mass["measured_positive"] + mass["unknown"],
    )
    output = {key: str(value) for ((key, value)) in (mass.items())}
    output["capture"] = {
        "known_contribution": str(lower),
        "lower": str(lower),
        "upper": str(upper),
        "per_offered_slot": {
            "known_contribution": str(lower / 10),
            "lower": str(lower / 10),
            "upper": str(upper / 10),
        },
    }
    return output

def allocation_metrics(ledger):
    pairs = [
        (rational(row["old_weight"]), rational(row["new_weight"])) for (row) in (ledger)
    ]
    old = sum((a for ((a, b)) in (pairs)), Fraction())
    new = sum((b for ((a, b)) in (pairs)), Fraction())
    shared = sum((min(a, b) for ((a, b)) in (pairs)), Fraction())
    lost, gained = (old - shared, new - shared)
    union = old + new - shared
    return {
        "old_mass": str(old),
        "new_mass": str(new),
        "shared_mass": str(shared),
        "lost_mass": str(lost),
        "gained_mass": str(gained),
        "union_mass": str(union),
        "half_l1": str((lost + gained) / 2),
        "weighted_jaccard": str(shared / union) if (union) else None,
        "undefined_jaccard_models": int(union == 0),
    }

def aggregate_metrics(records, scales):
    fields = (
        "old_mass",
        "new_mass",
        "shared_mass",
        "lost_mass",
        "gained_mass",
        "union_mass",
        "half_l1",
    )
    answer = {
        name: str(
            sum(
                (
                    scale * rational(record[name])
                    for ((record, scale)) in (zip(records, scales))
                ),
                Fraction(),
            )
        )
        for (name) in (fields)
    }
    answer["weighted_jaccard"] = (
        None
        if (any((record["weighted_jaccard"] is None for (record) in (records))))
        else str(
            sum(
                (
                    scale * rational(record["weighted_jaccard"])
                    for ((record, scale)) in (zip(records, scales))
                ),
                Fraction(),
            )
        )
    )
    answer["undefined_jaccard_models"] = sum(
        (record["undefined_jaccard_models"] for (record) in (records))
    )
    return answer

def summarise(models, scales, calls):
    need(
        len(models) == len(scales) and sum(scales, Fraction()) == 1,
        "Complete normalised hierarchical model weights required",
    )
    old, new, coefficients = ({}, {}, {})
    parts = {name: {} for (name) in (PARTITIONS)}
    for (model, scale) in (zip(models, scales)):
        for (row) in (model["ledger"]):
            key = (model["sample_id"], row["gene"])
            first, second, coefficient = (
                scale * rational(row["old_weight"]),
                scale * rational(row["new_weight"]),
                scale * rational(row["coefficient"]),
            )
            old[key] = old.get(key, Fraction()) + first
            new[key] = new.get(key, Fraction()) + second
            coefficients[key] = coefficients.get(key, Fraction()) + coefficient
            target = parts[row["partition"]]
            target[key] = target.get(key, Fraction()) + coefficient
    return {
        "cells": {"old": cell_mass(old, calls), "new": cell_mass(new, calls)},
        "bridge": bound_stats(coefficients, calls),
        "partitions": {
            name: bound_stats(mapping, calls) for ((name, mapping)) in (parts.items())
        },
        "allocation": aggregate_metrics(
            [allocation_metrics(model["ledger"]) for (model) in (models)], scales
        ),
    }

def finish_from_prepared(protocol, frozen, old_labels, new_labels):
    need(
        frozen["schema"] == "dawn_middle_version_allocation_freeze_v1"
        and frozen["prelabel_freeze"] is True,
        "Bridge allocation freeze must precede linking labels",
    )
    need(
        protocol["schema"] == "dawn_middle_version_domain_bridge_v1"
        and protocol["cells"] == CELL_NAMES
        and (protocol["scope"]["seeds"] == SEEDS)
        and (protocol["scope"]["budget"] == 10),
        "Protocol contract differs",
    )
    need(
        frozen["budget"] == 10
        and frozen["master_seeds"] == SEEDS
        and (frozen["counts"]["models"] == protocol["scope"]["models"])
        and (frozen["counts"]["donors"] == protocol["scope"]["donors"]),
        "Frozen complete scope differs",
    )
    old_calls, new_calls = (
        cache_calls(old_labels, frozen),
        cache_calls(new_labels, frozen),
    )
    overlap = set(old_calls) & set(new_calls)
    need(
        all((old_calls[key] == new_calls[key] for (key) in (overlap))),
        "Overlapping accepted call states disagree",
    )
    requests = [
        (row["sample_id"], row["gene"]) for (row) in (frozen["requested_labels"])
    ]
    need(
        len(requests)
        == len(set(requests))
        == frozen["counts"]["requested_sample_genes"],
        "Exact selected request census differs",
    )
    combined = old_calls | new_calls
    need(
        set(requests) <= set(combined),
        "Required accepted selected call missing; no imputation or full-assay fallback",
    )
    calls = {key: combined[key] for (key) in (requests)}
    donor_models = {}
    for (model) in (frozen["models"]):
        donor_models.setdefault(model["donor_id"], []).append(model)
    need(len(donor_models) == protocol["scope"]["donors"], "Full donor support changed")
    output = {
        "schema": "dawn_middle_version_bridge_readout_v1",
        "data_class": frozen["data_class"],
        "counts": frozen["counts"],
        "domains": frozen["domains"],
        "budget": 10,
        "master_seeds": SEEDS,
        "comparison": "new-domain minus old-domain at nominal native_11_5",
        "uncertainty": protocol["uncertainty"],
        "interpretation": protocol["interpretation"],
        "seed_semantics": "Five deterministic aliases, not independent replicates",
        "models": [],
        "donors": [],
        "population": {},
        "donor_sign_counts": dict.fromkeys(SIGNS, 0),
        "changed_allocation_models": 0,
        "changed_allocation_model_ids": [],
        "jaccard_aggregation": "Shared/union allocation mass ratio per model, equal model mean within donor, then equal donor mean; not a ratio of pooled masses. Both-empty is undefined and is never dropped from the population.",
        "donor_ids_by_sign": {name: [] for (name) in (SIGNS)},
        "ledger": [],
        "label_accounting": {
            "old_cache_calls": len(old_calls),
            "new_cache_calls": len(new_calls),
            "overlapping_cache_calls": len(overlap),
            "linked_bridge_calls": len(calls),
            "new_full_assay_access": False,
        },
    }
    population_scales = {}
    for (donor, models) in (sorted(donor_models.items())):
        result = summarise(models, [Fraction(1, len(models))] * len(models), calls)
        output["donors"].append(
            {
                "donor_id": donor,
                "model_ids": [model["model_id"] for (model) in (models)],
                **result,
            }
        )
        output["donor_sign_counts"][result["bridge"]["sign"]] += 1
        output["donor_ids_by_sign"][result["bridge"]["sign"]].append(donor)
        for (model) in (models):
            population_scales[model["model_id"]] = Fraction(
                1, len(donor_models) * len(models)
            )
    for (model) in (frozen["models"]):
        output["models"].append(
            {key: model[key] for (key) in (("model_id", "sample_id", "donor_id"))}
            | summarise([model], [Fraction(1)], calls)
        )
        for (row) in (model["ledger"]):
            state, value = calls[model["sample_id"], row["gene"]]
            output["ledger"].append(
                {key: model[key] for (key) in (("model_id", "sample_id", "donor_id"))}
                | row
                | {
                    "state": state,
                    "value": value,
                    "model_population_weight": str(
                        population_scales[model["model_id"]]
                    ),
                }
            )
    output["ledger"].sort(key = lambda row: (row["sample_id"], row["gene"]))
    output["population"] = summarise(
        frozen["models"],
        [population_scales[model["model_id"]] for (model) in (frozen["models"])],
        calls,
    )
    output["changed_allocation_model_ids"] = [
        model["model_id"]
        for (model) in (output["models"])
        if (rational(model["allocation"]["half_l1"]) > 0)
    ]
    output["changed_allocation_models"] = len(output["changed_allocation_model_ids"])
    need(
        len(output["ledger"]) == len(calls),
        "Complete sample-gene ledger census differs",
    )
    return output

def analyse(
    protocol, old_freeze, new_freeze, old_labels, new_labels, old_domain, new_domain
):
    frozen = prepare_bridge(protocol, old_freeze, new_freeze, old_domain, new_domain)
    return finish_from_prepared(protocol, frozen, old_labels, new_labels)
