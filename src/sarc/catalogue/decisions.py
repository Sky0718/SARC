import itertools

from . import review_rules as rules
from . import selection
from .constants import (
    CONTEXT,
    ENDPOINTS,
    METHODS,
    NATIVE,
    REFERENCES,
    SCALES,
    SEEDS,
    TRANSITIONS,
    group_key,
)

STRATEGIES = ("strong", "candidate", "development_single_best", "previous_release_best")

def validate_choices(choices):
    expected = set(itertools.product(TRANSITIONS, STRATEGIES))
    observed = [(row["transition"], row["strategy"]) for (row) in (choices)]
    if (
        len(observed) != 8
        or set(observed) != expected
        or any(row["context"] != CONTEXT for (row) in (choices))
    ):
        raise ValueError(
            "Frozen input-only choices do not cover both complete decisions"
        )
    for (row) in (choices):
        prediction = row["prediction"]
        chosen = prediction["choice"]
        if (chosen not in METHODS and chosen is not None):
            raise ValueError("Frozen choice is outside the fixed portfolio")
        if ((chosen is None) != (
            prediction["status"] == "UNAVAILABLE_TRAINING_PORTFOLIO"
        )):
            raise ValueError("Abstention and prediction status differ")
    return {
        (row["transition"], row["strategy"]): row["prediction"]["choice"]
        for (row) in (choices)
    }

def selection_records(groups, choices, bootstrap):
    fixed = validate_choices(choices)
    population = groups[group_key(METHODS[0], NATIVE[0])]["population"]
    result = []
    for (transition) in (TRANSITIONS):
        network = selection.state_networks(transition)["b"]
        selected = {
            strategy: fixed[(transition, strategy)] for (strategy) in (STRATEGIES)
        } | {"uniform_random": None, "posthoc_oracle": None}
        for (endpoint) in (ENDPOINTS):
            rows = {
                method: groups[group_key(method, network)]["measurements"][endpoint]
                for (method) in (METHODS)
            }
            intervals = bootstrap.selection_intervals(rows, selected)
            coverage = selection.portfolio_summary(population, rows)
            for (strategy, choice) in (selected.items()):
                loss = (
                    selection.expected_random_loss(population, rows)
                    if (strategy == "uniform_random")
                    else selection.oracle_loss(population, rows)
                    if (strategy == "posthoc_oracle")
                    else selection.portfolio_loss(population, rows, choice)
                )
                result.append(
                    {
                        "held_context": CONTEXT,
                        "transition": transition,
                        "strategy": strategy,
                        "endpoint": list(endpoint),
                        "choice": choice,
                        "loss": loss,
                        "outcome_coverage": coverage,
                        "conditional_interval": intervals[strategy],
                        "frozen_input_only_choice": strategy in STRATEGIES,
                    }
                )
    if (len(result) != 384):
        raise ValueError("Two decisions by six strategies by 32 endpoints are required")
    return {
        "status": "PROTECTED_FROZEN_CHOICE_EVALUATION_NO_REFIT",
        "decision_count": 2,
        "decisions_are_independent_replicates": False,
        "training_performed": False,
        "endpoints": result,
    }

def review_strata(completed, references):
    indexed = {
        (
            value["cell"]["method"],
            value["cell"]["network_id"],
            value["cell"]["master_seed"],
        ): value
        for (value) in (completed.values())
        if (value["cell"]["network_id"] in NATIVE)
    }
    expected = {
        (method, network, seed)
        for (method) in (METHODS)
        for (network) in (NATIVE)
        for (seed) in (SEEDS if (method == "PRODIGY") else (None,))
    }
    if (set(indexed) != expected or len(indexed) != 21):
        raise ValueError(
            "Complete 21 native seeded states are required for protected review"
        )
    for (method, transition, support, k) in (itertools.product(
        METHODS, TRANSITIONS, ("common", "native"), SCALES
    )):
        old, new = ("native_" + value for (value) in (transition.split("_to_")))
        seeds = SEEDS if (method == "PRODIGY") else (None,)
        prepared = {}
        for (sample, source) in (references.items()):
            eligible = rules.identities(
                source[support + "_eligible_literal_gene_labels"]
            )
            pairs = []
            for (seed) in (seeds):
                before, after = (
                    indexed[(method, old, seed)],
                    indexed[(method, new, seed)],
                )
                pairs.append(
                    {
                        "master_seed": seed,
                        "old_cell_id": before["cell"]["cell_id"],
                        "new_cell_id": after["cell"]["cell_id"],
                        "features": rules.pair_features(
                            before["samples"][sample],
                            after["samples"][sample],
                            eligible,
                            k,
                        ),
                    }
                )
            prepared[sample] = pairs
        for (reference) in (REFERENCES):
            rows = []
            for (sample, pairs) in (prepared.items()):
                positive = rules.identities(
                    references[sample][support + "_reference_positive_labels"][
                        reference
                    ]
                )
                joined = [
                    pair | {"outcome": rules.pair_outcome(pair["features"], positive)}
                    for (pair) in (pairs)
                ]
                row = rules.reduce_seed_pairs(sample, method, joined)
                endpoint = (support, reference, k)
                previous = [
                    indexed[(method, old, seed)]["samples"][sample]["measurements"][
                        endpoint
                    ]
                    for (seed) in (seeds)
                ]
                following = [
                    indexed[(method, new, seed)]["samples"][sample]["measurements"][
                        endpoint
                    ]
                    for (seed) in (seeds)
                ]
                before_hits = selection.mean(item["hits"] for (item) in (previous))
                after_hits = selection.mean(item["hits"] for (item) in (following))
                direct = {
                    "old_hits": before_hits,
                    "new_hits": after_hits,
                    "delta_hits": None
                    if (before_hits is None or after_hits is None)
                    else after_hits - before_hits,
                    "delta_lower": selection.mean(
                        after["lower"] - before["upper"]
                        for ((before, after)) in (zip(previous, following))
                    ),
                    "delta_upper": selection.mean(
                        after["upper"] - before["lower"]
                        for ((before, after)) in (zip(previous, following))
                    ),
                }
                if (any(
                    row["outcome"][name] != value for ((name, value)) in (direct.items())
                ) or row["actionable"] != (
                    before_hits is not None and after_hits is not None
                )):
                    raise ValueError(
                        "Independent direct-score and structural identities disagree"
                    )
                rows.append(row)
            identity = "::".join(
                (CONTEXT, method, transition, support, reference, str(k))
            )
            yield (
                {
                    "stratum_id": identity,
                    "context": CONTEXT,
                    "method": method,
                    "transition": transition,
                    "support": support,
                    "reference": reference,
                    "k": k,
                    "rows": rows,
                }
            )

def review_records(stratum, frozen_model, bootstrap):
    if (set(frozen_model) != {"strong", "candidate"}):
        raise ValueError("Both frozen development review families are required")
    curves = rules.curves(stratum, frozen_model)
    allocations = {}
    for (curve) in (curves["curves"]):
        key = curve["rule"] + "::" + str(curve["allocation"]["budget_fraction"])
        existing = allocations.setdefault(key, curve["allocation"])
        if (existing != curve["allocation"]):
            raise ValueError("Review allocation was changed with the event threshold")
    if (len(allocations) != 21 or len(curves["curves"]) != 63):
        raise ValueError("Full frozen review budget and event grid is required")
    intervals = {
        delta: bootstrap.review_intervals(stratum["rows"], allocations, delta)
        for (delta) in (rules.DELTAS)
    }
    for (curve) in (curves["curves"]):
        key = curve["rule"] + "::" + str(curve["allocation"]["budget_fraction"])
        curve["conditional_interval"] = intervals[curve["endpoints"]["delta"]][key]
    return curves
