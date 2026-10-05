import itertools
from fractions import Fraction
from ..transfer import scoring as core
from . import review_rules as rules

SUPPORTS = ("common", "native")
REFERENCES = (
    "NCG6_primary_all",
    "NCG6_known_subgroup",
    "NCG6_candidate_subgroup",
    "CancerMine2019_secondary",
)
SCALES = (10, 1, 5, 20)
RELEASES = ("11_0", "11_5", "12_0")
PARAMETERS = {"DawnRank": "mu_3", "PRODIGY": "alpha_0.05", "PersonaDrive": "original"}

def accepted_core():
    return core

def native_cells(cells, references, registry):
    if (
        (not isinstance(cells, dict))
        or (not isinstance(registry, dict))
        or (len(registry) != 162)
        or (set(cells) != set(registry))
    ):
        raise ValueError("Complete accepted 162-cell registry is required")
    if (set(references) != set(rules.CONTEXTS)):
        raise ValueError("Only complete development reference contexts are allowed")
    populations = {}
    selected = {}
    for (cell_id, cell) in (registry.items()):
        if (
            (cell_id != cell["cell_id"])
            or (cell["context"] not in rules.CONTEXTS)
            or (cell["method"] not in rules.METHODS)
            or (cell["node"] not in ("NCS-S03", "NCS-S04"))
        ):
            raise ValueError("Foreign logical cell")
        order = cell["sample_order"]
        rules.identities(order)
        if ((len(order) != 36) or (list(references[cell["context"]]) != order)):
            raise ValueError("Frozen full sample order differs")
        populations.setdefault(cell["context"], order)
        if (populations[cell["context"]] != order):
            raise ValueError("Cells do not share the frozen population")
        record = cells[cell_id]
        if ((record["cell"] != cell) or (list(record["samples"]) != order)):
            raise ValueError("Accepted cell metadata or full sample order differs")
        if (
            (cell["node"] != "NCS-S03")
            or (cell["parameter_id"] != PARAMETERS[cell["method"]])
            or (
                cell["network_id"]
                not in tuple("native_" + release for (release) in (RELEASES))
            )
        ):
            continue
        seed = cell["master_seed"]
        expected = rules.SEEDS if (cell["method"] == "PRODIGY") else (None,)
        if (seed not in expected):
            raise ValueError("Unexpected primary native seed")
        key = (cell["context"], cell["method"], cell["network_id"][7:], seed)
        if (key in selected):
            raise ValueError("Duplicate primary native state")
        selected[key] = record
    expected = {
        (context, method, release, seed)
        for ((context, method, release)) in (
            itertools.product(rules.CONTEXTS, rules.METHODS, RELEASES)
        )
        for (seed) in (rules.SEEDS if (method == "PRODIGY") else (None,))
    }
    if ((set(selected) != expected) or (set(populations) != set(rules.CONTEXTS))):
        raise ValueError("Incomplete 42-cell primary native population")
    return selected, populations

def eligible_labels(source, support):
    return rules.identities(source[support + "_eligible_literal_gene_labels"])

def positive_labels(source, support, reference):
    catalogues = source[support + "_reference_positive_labels"]
    if (set(catalogues) != set(REFERENCES)):
        raise ValueError("Frozen four-catalogue population differs")
    positive = rules.identities(catalogues[reference])
    if (not positive.issubset(eligible_labels(source, support))):
        raise ValueError("Positive reference is outside frozen eligibility")
    return positive

def stratum_identity(context, method, transition, support, reference, k):
    return "::".join((context, method, transition, support, reference, str(k)))

def prepare_features(cells, references, registry):
    selected, populations = native_cells(cells, references, registry)
    strata = []
    for (context, method, transition, support, k) in (itertools.product(
        rules.CONTEXTS, rules.METHODS, rules.TRANSITIONS, SUPPORTS, SCALES
    )):
        old, new = transition.split("_to_")
        seeds = rules.SEEDS if (method == "PRODIGY") else (None,)
        rows = []
        for (sample) in (populations[context]):
            eligible = eligible_labels(references[context][sample], support)
            pairs = []
            for (seed) in (seeds):
                previous = selected[(context, method, old, seed)]
                following = selected[(context, method, new, seed)]
                features = rules.pair_features(
                    previous["samples"][sample],
                    following["samples"][sample],
                    eligible,
                    k,
                )
                pairs.append(
                    {
                        "master_seed": seed,
                        "old_cell_id": previous["cell"]["cell_id"],
                        "new_cell_id": following["cell"]["cell_id"],
                        "features": features,
                    }
                )
            rows.append({"sample_id": sample, "pairs": pairs})
        strata.append(
            {
                "context": context,
                "method": method,
                "transition": transition,
                "support": support,
                "k": k,
                "rows": rows,
            }
        )
    return {
        "schema": "s05_review_label_free_features_v1",
        "status": "OUTPUT_DERIVED_FEATURES_NO_REFERENCE_LABEL_JOIN",
        "consumed_cell_ids": [
            record["cell"]["cell_id"] for (record) in (selected.values())
        ],
        "sample_orders": populations,
        "strata": strata,
    }

def attach_outcomes(prepared, references):
    if ((prepared["schema"] != "s05_review_label_free_features_v1") or (
        set(references) != set(rules.CONTEXTS)
    )):
        raise ValueError("Invalid prepared review feature record")
    output = []
    for (source) in (prepared["strata"]):
        context, method, transition, support, k = (
            source[key]
            for (key) in (("context", "method", "transition", "support", "k"))
        )
        if (
            (context not in rules.CONTEXTS)
            or (method not in rules.METHODS)
            or (transition not in rules.TRANSITIONS)
            or (support not in SUPPORTS)
            or (k not in SCALES)
        ):
            raise ValueError("Foreign prepared review stratum")
        if ([row["sample_id"] for (row) in (source["rows"])] != prepared[
            "sample_orders"
        ][context]):
            raise ValueError("Prepared sample order changed")
        for (reference) in (REFERENCES):
            rows = []
            for (row) in (source["rows"]):
                sample = row["sample_id"]
                ref = references[context][sample]
                positive = positive_labels(ref, support, reference)
                pairs = []
                for (pair) in (row["pairs"]):
                    features = pair["features"]
                    if (features["_eligible"] != eligible_labels(ref, support)):
                        raise ValueError(
                            "Reference join changed non-release eligibility"
                        )
                    pairs.append(
                        pair | {"outcome": rules.pair_outcome(features, positive)}
                    )
                rows.append(rules.reduce_seed_pairs(sample, method, pairs))
            output.append(
                {
                    "stratum_id": stratum_identity(
                        context, method, transition, support, reference, k
                    ),
                    "context": context,
                    "method": method,
                    "transition": transition,
                    "support": support,
                    "reference": reference,
                    "k": k,
                    "rows": rows,
                }
            )
    expected = set(
        itertools.product(
            rules.CONTEXTS, rules.METHODS, rules.TRANSITIONS, SUPPORTS, SCALES
        )
    )
    observed = [
        (row["context"], row["method"], row["transition"], row["support"], row["k"])
        for (row) in (prepared["strata"])
    ]
    if ((len(observed) != len(expected)) or (set(observed) != expected)):
        raise ValueError("Incomplete prepared review feature population")
    return {
        "schema": "s05_review_strata_v1",
        "status": "COMPLETE_DEVELOPMENT_REVIEW_INPUTS_NOT_FITTED",
        "stratum_count": len(output),
        "sample_orders": prepared["sample_orders"],
        "consumed_cell_ids": prepared["consumed_cell_ids"],
        "strata": output,
    }

def prepare_review(cells, references, registry):
    return attach_outcomes(prepare_features(cells, references, registry), references)

def direct_reference(cells, references, registry):
    selected, populations = native_cells(cells, references, registry)
    score_list = accepted_core().score_list
    scores = {}
    for (key, cell) in (selected.items()):
        context = key[0]
        for (sample) in (populations[context]):
            output = cell["samples"][sample]
            state = output["method_status"]
            if (state not in ("SUCCESS", "SUCCESS_EMPTY", "FAILED", "UNAVAILABLE")):
                raise ValueError("Unknown method state")
            if (
                (
                    (state in ("SUCCESS", "SUCCESS_EMPTY"))
                    != (output["genes"] is not None)
                )
                or ((state == "SUCCESS") and (not output["genes"]))
                or ((state == "SUCCESS_EMPTY") and bool(output["genes"]))
            ):
                raise ValueError("Method state and native roster disagree")
            for (support, reference, k) in (itertools.product(
                SUPPORTS, REFERENCES, SCALES
            )):
                source = references[context][sample]
                scores[(key, sample, support, reference, k)] = score_list(
                    output["genes"],
                    output["scores"],
                    eligible_labels(source, support),
                    positive_labels(source, support, reference),
                    k,
                )
    strata = []
    for (context, method, transition, support, reference, k) in (itertools.product(
        rules.CONTEXTS, rules.METHODS, rules.TRANSITIONS, SUPPORTS, REFERENCES, SCALES
    )):
        old, new = transition.split("_to_")
        seeds = rules.SEEDS if (method == "PRODIGY") else (None,)
        rows = []
        for (sample) in (populations[context]):
            previous = [
                scores[((context, method, old, seed), sample, support, reference, k)]
                for (seed) in (seeds)
            ]
            following = [
                scores[((context, method, new, seed), sample, support, reference, k)]
                for (seed) in (seeds)
            ]
            old_hits = (
                None
                if (any(row["hits"] is None for (row) in (previous)))
                else sum((row["hits"] for (row) in (previous)), Fraction(0))
                / len(seeds)
            )
            new_hits = (
                None
                if (any(row["hits"] is None for (row) in (following)))
                else sum((row["hits"] for (row) in (following)), Fraction(0))
                / len(seeds)
            )
            outcome = {
                "old_hits": old_hits,
                "new_hits": new_hits,
                "delta_hits": None
                if ((old_hits is None) or (new_hits is None))
                else new_hits - old_hits,
                "delta_lower": sum(
                    (
                        after["lower"] - before["upper"]
                        for ((before, after)) in (zip(previous, following))
                    ),
                    Fraction(0),
                )
                / len(seeds),
                "delta_upper": sum(
                    (
                        after["upper"] - before["lower"]
                        for ((before, after)) in (zip(previous, following))
                    ),
                    Fraction(0),
                )
                / len(seeds),
                "positive_count": previous[0]["positive_count"],
            }
            rows.append(
                {
                    "sample_id": sample,
                    "actionable": (old_hits is not None) and (new_hits is not None),
                    "master_seeds": list(seeds),
                    "expected_seed_count": len(seeds),
                    "complete_seed_count": sum(
                        (before["hits"] is not None) and (after["hits"] is not None)
                        for ((before, after)) in (zip(previous, following))
                    ),
                    "features": None,
                    "outcome": outcome,
                }
            )
        strata.append(
            {
                "stratum_id": stratum_identity(
                    context, method, transition, support, reference, k
                ),
                "context": context,
                "method": method,
                "transition": transition,
                "support": support,
                "reference": reference,
                "k": k,
                "rows": rows,
            }
        )
    return {
        "schema": "s05_direct_reference_strata_v1",
        "status": "EXACT_RAW_LIST_CATALOGUE_PROCESSING_NO_STRUCTURAL_FEATURES",
        "stratum_count": len(strata),
        "sample_orders": populations,
        "strata": strata,
    }

def conditional_review_summary(rows, allocation, indices, delta = 1):
    core = accepted_core()
    if (indices != core.bootstrap_indices(len(rows), 2000, 20260926)):
        raise ValueError("Conditional intervals require exact frozen F2 paired indices")
    point = rules.budget_endpoints(rows, allocation, delta)
    names = (
        "captured_yield",
        "captured_yield_lower",
        "captured_yield_upper",
        "precision",
        "known_event_capture_rate",
        "full_event_capture_rate",
        "full_event_capture_lower",
        "full_event_capture_upper",
        "full_event_rate",
    )
    values = {name: [] for (name) in (names)}
    counts = {name: [] for (name) in (point["counts"])}
    for (draw) in (indices):
        result = rules.reweight_allocation(rows, allocation, draw, delta)
        for (name) in (names):
            values[name].append(result[name])
        for (name) in (counts):
            counts[name].append(result["counts"][name])
    intervals = {}
    for (name, vector) in ((values | counts).items()):
        available = [float(value) for (value) in (vector) if (value is not None)]
        intervals[name] = {
            "available_draws": len(available),
            "draws": len(indices),
            "lower_95": core.percentile(sorted(available), 0.025)
            if (len(available) == len(indices))
            else None,
            "upper_95": core.percentile(sorted(available), 0.975)
            if (len(available) == len(indices))
            else None,
            "status": "CONDITIONAL_FIXED_ALLOCATION_PERCENTILE"
            if (len(available) == len(indices))
            else "UNAVAILABLE_ZERO_DENOMINATOR_OR_MISSING_DRAW",
        }
    return {
        "status": "CONDITIONAL_SAMPLE_COMPOSITION_NOT_TRAINED_RULE_GENERALISATION",
        "point": point,
        "draw_count": 2000,
        "seed": 20260926,
        "allocation_recomputed": False,
        "parameters_refitted": False,
        "allocated_sample_ids": list(allocation["allocated"]),
        "intervals": intervals,
    }
