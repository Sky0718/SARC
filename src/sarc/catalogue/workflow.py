from collections import Counter

from . import (
    comparisons,
    decisions,
    development_details,
    features,
    groups,
    matrix,
    optimism,
    references,
    review_adapter,
    review_rules,
    selection,
    selection_adapter,
)
from .bootstrap import PairedBootstrap
from .constants import (
    ENDPOINTS,
    FACTORS,
    METHODS,
    PARAMETERS,
    REFERENCES,
    SEEDS,
    TRANSITIONS,
    network_states,
)
from .encoding import decode_rationals

def registry_index(registry, development):
    records = (
        registry["cells"]
        if (isinstance(registry, dict) and "cells" in registry)
        else registry
    )
    if (isinstance(records, dict)):
        records = list(records.values())
    index = {row["cell_id"]: row for (row) in (records)}
    if (len(index) != len(records) or len(index) != (162 if (development) else 69)):
        raise ValueError("The complete distinct catalogue cell registry is required")
    if (development):
        expected = set()
        for (context) in (selection.CONTEXTS):
            for (method) in (METHODS):
                for (network) in (network_states(method)):
                    for (seed) in (SEEDS if (method == "PRODIGY") else (None,)):
                        expected.add(
                            (context, method, network, PARAMETERS[method], seed)
                        )
            for (network) in (("native_11_0", "native_11_5", "native_12_0")):
                for (parameter) in (("mu_1", "mu_10")):
                    expected.add((context, "DawnRank", network, parameter, None))
                for (parameter) in (("alpha_0.01", "alpha_0.1")):
                    expected.add((context, "PRODIGY", network, parameter, 104729))
        observed = [
            (
                row["context"],
                row["method"],
                row["network_id"],
                row["parameter_id"],
                row["master_seed"],
            )
            for (row) in (records)
        ]
        if (set(observed) != expected or len(observed) != len(expected)):
            raise ValueError(
                "The complete original development factor grid is required"
            )
    return index

def complete_scores(registry, raw_cells, reference_rows, development):
    registry = registry_index(registry, development)
    contexts = selection.CONTEXTS if (development) else ("COAD_TCGA",)
    if (set(reference_rows) != set(contexts)):
        raise ValueError(
            "References must retain the exact development or protected cohorts"
        )
    for (context) in (contexts):
        population = list(reference_rows[context])
        if (len(population) != (36 if (development) else 396)):
            raise ValueError("The complete source sample population is required")
        matrix.validate_references(reference_rows[context], population)
    if (not development):
        matrix.validate_registry(registry, list(reference_rows["COAD_TCGA"]))
    records = {row["cell_id"]: row["samples"] for (row) in (raw_cells)}
    if (len(records) != len(raw_cells) or set(records) != set(registry)):
        raise ValueError(
            "Exactly one complete raw result per registered cell is required"
        )
    completed = {}
    for (cell_id, cell) in (registry.items()):
        if (cell["sample_order"] != list(reference_rows[cell["context"]])):
            raise ValueError("Cell and reference sample order differs")
        raw = records[cell_id]
        if (isinstance(raw, list)):
            identities = [row["sample_id"] for (row) in (raw)]
            if (identities != cell["sample_order"]):
                raise ValueError(
                    "Raw output does not preserve the complete ordered population"
                )
            raw = {row["sample_id"]: row for (row) in (raw)}
        completed[cell_id] = matrix.score_cell(
            cell, raw, reference_rows[cell["context"]]
        )
    return registry, completed

def sample_metrics(completed):
    result = []
    for (value) in (completed.values()):
        cell = value["cell"]
        for (sample, row) in (value["samples"].items()):
            for ((support, reference, k), measurement) in (row["measurements"].items()):
                result.append(
                    {
                        "cell_id": cell["cell_id"],
                        **{name: cell[name] for (name) in (FACTORS)},
                        "master_seed": cell["master_seed"],
                        "sample_id": sample,
                        "method_status": row["method_status"],
                        "failure_reason": row.get("reason"),
                        "support": support,
                        "reference": reference,
                        **measurement,
                    }
                )
    return result

def cohort_metrics(reduced):
    return [
        {
            **dict(zip(FACTORS, key)),
            "seeds": group["expected_seeds"],
            "support": endpoint[0],
            "reference": endpoint[1],
            "k": endpoint[2],
            "retrieval": summary,
            "recall": groups.recall_summary(group["measurements"][endpoint]),
        }
        for ((key, group)) in (reduced.items())
        for ((endpoint, summary)) in (group["cohort_summaries"].items())
    ]

def coverage(completed):
    states = Counter(
        row["method_status"]
        for (value) in (completed.values())
        for (row) in (value["samples"].values())
    )
    return {
        "logical_cells": len(completed),
        "sample_cell_states": sum(states.values()),
        "native_status_counts": dict(states),
        "seeds_reduced_before_samples": True,
        "catalogue_recovery_is_not_functional_validation": True,
    }

def development(registry, raw_cells, reference_rows, feature_records):
    feature_records = decode_rationals(feature_records)
    registry, completed = complete_scores(registry, raw_cells, reference_rows, True)
    values = list(completed.values())
    reduced = groups.group_cells(values)
    contrasts, reversals = groups.comparisons(reduced)
    variability, seeded = development_details.variability_records(values, reduced)
    dataset = selection_adapter.prepare_selection(completed, reference_rows, registry)
    selectors = selection.evaluate_development(dataset, feature_records)
    strata = review_adapter.prepare_review(completed, reference_rows, registry)[
        "strata"
    ]
    return {
        "coverage": coverage(completed),
        "sample_metrics": sample_metrics(completed),
        "cohort_summaries": cohort_metrics(reduced),
        "paired_contrasts": contrasts,
        "reversal_margins": reversals,
        "secondary_retuning": groups.tuning_trials(values),
        "leading_membership_stability": development_details.membership_records(
            values, reference_rows
        ),
        "stochastic_sample_variability": variability,
        "paired_seed_release_effects": seeded,
        "selection": selectors,
        "review": review_rules.evaluate_review(strata),
        "leave_one_context_out_review": review_rules.loco_review(strata),
    }

def protected(registry, raw_cells, reference_rows, frozen_choices, frozen_review):
    frozen_choices = decode_rationals(frozen_choices)
    frozen_review = decode_rationals(frozen_review)
    registry, completed = complete_scores(registry, raw_cells, reference_rows, False)
    source = reference_rows["COAD_TCGA"]
    reduced = matrix.reduce_cells(completed, registry, source)
    bootstrap = PairedBootstrap(list(source))
    contrasts, reversals = comparisons.comparisons(reduced, bootstrap)
    strata = decisions.review_strata(completed, source)
    return {
        "coverage": coverage(completed),
        "sample_metrics": sample_metrics(completed),
        "cohort_summaries": cohort_metrics(reduced),
        "paired_contrasts": contrasts,
        "reversal_margins": reversals,
        "leading_membership_stability": list(
            comparisons.membership_records(completed, source)
        ),
        "stochastic_sample_variability": list(
            comparisons.variability_records(completed, reduced, bootstrap)
        ),
        "paired_seed_release_effects": list(
            comparisons.seeded_contrasts(completed, bootstrap)
        ),
        "selection": decisions.selection_records(reduced, frozen_choices, bootstrap),
        "review": [
            decisions.review_records(stratum, frozen_review, bootstrap)
            for (stratum) in (strata)
        ],
    }

def cohort_features(samples, old_network, new_network, persistent):
    return features.cohort_features(samples, old_network, new_network, persistent)

def frozen_choices(selection_model, target_features):
    model = decode_rationals(selection_model)
    target = decode_rationals(target_features)
    if (set(target) != set(TRANSITIONS)):
        raise ValueError("Both protected input-only transitions are required")
    result = []
    for (transition) in (TRANSITIONS):
        for (strategy) in (decisions.STRATEGIES):
            prediction = (
                selection.predict_frozen(
                    model, target[transition], transition, strategy
                )
                if (strategy in selection.FAMILIES)
                else model["transitions"][transition]["simple_rules"][strategy]
            )
            result.append(
                {
                    "context": "COAD_TCGA",
                    "transition": transition,
                    "strategy": strategy,
                    "prediction": prediction,
                }
            )
    decisions.validate_choices(result)
    return result

def reference_join(primary_rows, secondary_rows, cohort_inputs):
    known = {
        row["symbol"] for (row) in (primary_rows) if (row["type"] == "Known Cancer")
    }
    if (len(known) != 711):
        raise ValueError("The original NCG6 global known-gene class contains 711 genes")
    catalogues = references.build_catalogues(primary_rows, secondary_rows, known)
    output = {}
    for (cohort) in (cohort_inputs):
        context = cohort["context"]
        if (context not in catalogues or context in output):
            raise ValueError("Unknown or repeated catalogue context")
        samples = cohort["sample_order"]
        if (len(samples) != (396 if (context == "COAD_TCGA") else 36) or len(
            set(samples)
        ) != len(samples)):
            raise ValueError("Complete original sample order is required")
        mutation_rows = cohort["mutation_rows"]
        labels = [row["gene"] for (row) in (mutation_rows)]
        states = references.literal_states(labels)
        if (any(
            len(row["values"]) != len(samples)
            or any(value not in (0, 1) for (value) in (row["values"]))
            for (row) in (mutation_rows)
        )):
            raise ValueError(
                "Binary mutation rows must cover the complete ordered samples"
            )
        common_support = set(cohort["common_nonrelease_identity_candidate_genes"])
        result = {}
        for (index, sample) in (enumerate(samples)):
            original = {
                row["gene"] for (row) in (mutation_rows) if (row["values"][index] != 0)
            }
            valid = {gene for (gene) in (original) if (not states[gene])}
            common = original & common_support
            if (not common <= valid):
                raise ValueError(
                    "Common eligibility contains an unresolved literal mutation label"
                )
            result[sample] = {
                "sample": sample,
                "original_nonzero_mutation_labels": sorted(original),
                "native_eligible_literal_gene_labels": sorted(valid),
                "native_unresolved_mutation_labels": [
                    {"literal_label": gene, "states": states[gene]}
                    for (gene) in (sorted(original - valid))
                ],
                "native_valid_but_common_excluded_labels": sorted(valid - common),
                "common_eligible_literal_gene_labels": sorted(common),
                "common_reference_positive_labels": {
                    name: sorted(common & genes)
                    for ((name, genes)) in (catalogues[context].items())
                },
                "native_reference_positive_labels": {
                    name: sorted(valid & genes)
                    for ((name, genes)) in (catalogues[context].items())
                },
            }
        output[context] = result
    return output

def mean_bank(records, populations, biological_units):
    contexts = {}
    for (row) in (records):
        context = row["context"]
        endpoint = tuple(row["endpoint"])
        candidate = tuple(row["candidate"])
        target = contexts.setdefault(context, {}).setdefault(endpoint, {})
        if (candidate in target):
            raise ValueError("Repeated context-endpoint-candidate vector")
        target[candidate] = decode_rationals(row["values"])
    return optimism.MeanBank(contexts, populations, biological_units)

def selection_optimism(
    records, populations, biological_units, designs, ordinary_records = None
):
    bank = mean_bank(records, populations, biological_units)
    ordinary = None
    if (ordinary_records is not None):
        contexts = {row["context"] for (row) in (ordinary_records)}
        ordinary = mean_bank(
            ordinary_records,
            {name: populations[name] for (name) in (contexts)},
            {name: biological_units[name] for (name) in (contexts)},
        )
    return optimism.evaluate(bank, designs, 2000, 20260927, ordinary)
