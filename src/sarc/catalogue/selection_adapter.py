from fractions import Fraction

from .selection import (
    CONTEXTS,
    METHODS,
    NATIVE,
    PRIMARY,
    TRANSITIONS,
    mean,
    state_networks,
)

REFERENCES = (
    "NCG6_primary_all",
    "NCG6_known_subgroup",
    "NCG6_candidate_subgroup",
    "CancerMine2019_secondary",
)
SCALES = (10, 1, 5, 20)
SEEDS = (104729, 130363, 155921, 196613, 228017)
PARAMETERS = {"DawnRank": "mu_3", "PersonaDrive": "original", "PRODIGY": "alpha_0.05"}
ENDPOINTS = tuple(
    (support, reference, k)
    for (support) in (("common", "native"))
    for (reference) in (REFERENCES)
    for (k) in (SCALES)
)

def project_input_samples(references):
    if (set(references) != set(CONTEXTS)):
        raise ValueError("Input projection is development-only")
    result = {}
    for (context) in (CONTEXTS):
        if (len(references[context]) != 36):
            raise ValueError("The complete 36-sample input population is required")
        result[context] = [
            {
                "sample_id": sample,
                "eligible_ids": list(row["common_eligible_literal_gene_labels"]),
            }
            for ((sample, row)) in (references[context].items())
        ]
    return result

def prepare_selection(cells, references, registry):
    if (
        (len(registry) != 162)
        or (set(cells) != set(registry))
        or (set(references) != set(CONTEXTS))
    ):
        raise ValueError(
            "The accepted complete 162-cell development matrix is required"
        )
    populations = {context: tuple(references[context]) for (context) in (CONTEXTS)}
    if (any(
        (len(population) != 36) or (len(set(population)) != 36)
        for (population) in (populations.values())
    )):
        raise ValueError("Frozen sample population is incomplete")
    grouped = {}
    for (cell_id, result) in (cells.items()):
        cell = result["cell"]
        if (
            (cell != registry[cell_id])
            or (cell["cell_id"] != cell_id)
            or (cell["context"] not in CONTEXTS)
            or (cell["method"] not in METHODS)
        ):
            raise ValueError("Misrouted or protected cell")
        population = populations[cell["context"]]
        if ((tuple(cell["sample_order"]) != population) or (
            tuple(result["samples"]) != population
        )):
            raise ValueError("Cell and reference populations differ")
        for (sample) in (population):
            current = result["samples"][sample]
            status = current["method_status"]
            if (status not in ("SUCCESS", "SUCCESS_EMPTY", "FAILED", "UNAVAILABLE")):
                raise ValueError("Unknown method status")
            valid = status in ("SUCCESS", "SUCCESS_EMPTY")
            if ((
                valid and ((current["genes"] is None) or (current["scores"] is None))
            ) or (
                (not valid)
                and ((current["genes"] is not None) or (current["scores"] is not None))
            )):
                raise ValueError("A failed or successful raw state is misrepresented")
            if (valid and (
                (len(current["genes"]) != len(current["scores"]))
                or (bool(current["genes"]) != (status == "SUCCESS"))
            )):
                raise ValueError("Successful empty output is misrepresented")
            if (set(current["measurements"]) != set(ENDPOINTS)):
                raise ValueError("Required primary or secondary endpoints are missing")
            if (any(
                (measurement["hits"] is not None) != valid
                for (measurement) in (current["measurements"].values())
            )):
                raise ValueError("Unavailable performance must not be filled with zero")
        key = (
            cell["context"],
            cell["method"],
            cell["network_id"],
            cell["parameter_id"],
        )
        by_seed = grouped.setdefault(key, {})
        if (cell["master_seed"] in by_seed):
            raise ValueError("Repeated seed within a cohort cell group")
        by_seed[cell["master_seed"]] = result
    for (key, by_seed) in (grouped.items()):
        expected = (
            SEEDS
            if (key[1] == "PRODIGY" and key[3] == "alpha_0.05")
            else (104729,)
            if (key[1] == "PRODIGY")
            else (None,)
        )
        if (set(by_seed) != set(expected)):
            raise ValueError("Incomplete registered repeat group")
    rows = {}
    means = {}
    networks = tuple(
        dict.fromkeys(
            NATIVE
            + tuple(
                network
                for (transition) in (TRANSITIONS)
                for (network) in (state_networks(transition).values())
            )
        )
    )
    for (context) in (CONTEXTS):
        for (method) in (METHODS):
            for (network) in (networks):
                key = (context, method, network, PARAMETERS[method])
                if (key not in grouped):
                    raise ValueError(
                        "A required native or actual input-intervention group is absent"
                    )
                expected = SEEDS if (method == "PRODIGY") else (None,)
                for (endpoint) in (ENDPOINTS):
                    support, reference, k = endpoint
                    reduced = {}
                    for (sample) in (populations[context]):
                        reference_row = references[context][sample]
                        eligible = frozenset(
                            reference_row[support + "_eligible_literal_gene_labels"]
                        )
                        positives = frozenset(
                            reference_row[support + "_reference_positive_labels"][
                                reference
                            ]
                        )
                        if (not positives.issubset(eligible)):
                            raise ValueError(
                                "Reference positives exceed frozen eligibility"
                            )
                        values = [
                            grouped[key][seed]["samples"][sample]["measurements"][
                                endpoint
                            ]
                            for (seed) in (expected)
                        ]
                        if (any(
                            (value["eligible_count"] != len(eligible))
                            or (value["positive_count"] != len(positives))
                            or (value["k"] != k)
                            for (value) in (values)
                        )):
                            raise ValueError(
                                "Per-seed support or reference identity changed"
                            )
                        hits = mean(value["hits"] for (value) in (values))
                        reduced[sample] = {
                            "hits": hits,
                            "lower": mean(value["lower"] for (value) in (values)),
                            "upper": mean(value["upper"] for (value) in (values)),
                            "eligible_ids": eligible,
                            "positive_ids": positives,
                            "k": k,
                            "recall": mean(value["recall"] for (value) in (values)),
                            "repeat_count": len(expected),
                            "complete_repeat_count": sum(
                                value["hits"] is not None for (value) in (values)
                            ),
                            "seed_values": {
                                str(seed): value["hits"]
                                for ((seed, value)) in (zip(expected, values))
                            },
                            "status": (
                                "SUCCESS"
                                if (hits is not None)
                                else "UNAVAILABLE_REQUIRED_REPEAT"
                            ),
                            "native_output_count": mean(
                                value["native_output_count"] for (value) in (values)
                            ),
                            "common_output_count": mean(
                                value["common_output_count"] for (value) in (values)
                            ),
                            "positive_output_count": mean(
                                value["positive_output_count"] for (value) in (values)
                            ),
                        }
                    rows[(context, method, network, endpoint)] = reduced
                    means[(context, method, network, endpoint)] = mean(
                        record["hits"] for (record) in (reduced.values())
                    )
    return {
        "contexts": CONTEXTS,
        "populations": populations,
        "endpoints": ENDPOINTS,
        "rows": rows,
        "cohort_means": means,
        "complete_input_cells": 162,
        "consumed_primary_groups": (len(CONTEXTS) * len(METHODS) * len(networks)),
        "retuning_results_not_used_for_rule_choice": True,
    }
