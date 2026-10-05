import collections
import itertools
import json
from fractions import Fraction
from .continuous import (
    require,
    number,
    mean,
    key,
    metric,
    paired_seed,
    aggregate_seeds,
    seed_list,
    METHODS,
    POOLS,
    BUDGETS,
    POLICIES,
    PRIMARY_PAIRS,
    DEGREE_PAIRS,
    KEY_FIELDS,
)

MODELS = tuple(
    (
        "ACH-" + suffix
        for (suffix) in (
            (
                "003332",
                "003334",
                "003337",
                "003346",
                "003387",
                "003400",
                "003401",
                "003402",
            )
        )
    )
)
MASTERS = (104729, 130363, 155921, 196613, 228017)

def validate_endpoints(endpoints):
    require(
        isinstance(endpoints, dict) and set(endpoints).issubset(MODELS),
        "Endpoints restricted to exact eight models",
    )
    support = {}
    for (model) in (MODELS):
        scores = endpoints.get(model, {})
        require(isinstance(scores, dict), "Model endpoint must be gene to score/null")
        require(
            all((isinstance(g, str) and g and (g == g.strip()) for (g) in (scores))),
            "Literal canonical gene keys required",
        )
        require(
            all((v is None or number(v) for (v) in (scores.values()))),
            "Nonfinite/string/Boolean scores rejected",
        )
        finite = sorted((g for ((g, v)) in (scores.items()) if (v is not None)))
        support[model] = dict(
            sample_endpoint_record_present = model in endpoints,
            assessed_gene_count = len(finite),
            finite_zero_count = sum(
                (v == 0 for (v) in (scores.values()) if (v is not None))
            ),
            assessed_genes = finite,
            score_threshold_used = False,
        )
    return support

def full_model_summary(rows, value_field):
    require(
        [r["sample_id"] for (r) in (rows)] == list(MODELS),
        "Exact eight-model order required",
    )
    available = [r for (r) in (rows) if (r[value_field] is not None)]
    values = [r[value_field] for (r) in (available)]
    full = len(available) == len(MODELS)
    return dict(
        status = "IDENTIFIED_DESCRIPTIVE_N8"
        if (full)
        else "UNIDENTIFIED_INCOMPLETE_EIGHT_MODEL_SUPPORT",
        equal_model_mean = mean(values) if (full) else None,
        model_range = [min(values), max(values)] if (full) else None,
        expected_model_count = 8,
        valid_complete_seed_model_count = len(available),
        descriptive_available_model_mean = mean(values),
        descriptive_available_model_ids = [r["sample_id"] for (r) in (available)],
        unavailable_model_ids = [
            r["sample_id"] for (r) in (rows) if (r[value_field] is None)
        ],
        model_results = rows,
        equal_weight_per_model = True,
        leave_one_model_out_means = [
            dict(
                omitted_sample_id = omitted["sample_id"],
                remaining_sample_ids = [
                    r["sample_id"] for (r) in (rows) if (r is not omitted)
                ],
                mean = mean([r[value_field] for (r) in (rows) if (r is not omitted)]),
            )
            for (omitted) in (rows)
        ]
        if (full)
        else None,
        inferential_p_value = None,
        confidence_interval = None,
        independent_confirmation_claim = False,
    )

def _bounds(scored, selection):
    scored["unassayed_returned_mass"] = (
        scored["R"] - scored["A"]
        if (scored["R"] is not None and scored["A"] is not None)
        else None
    )
    scored["unreturned_offered_mass"] = (
        scored["k"] - scored["R"] if (scored["R"] is not None) else None
    )
    all_returned_assayed = (
        scored["selection_valid"]
        and (not scored.get("unassessed_genes"))
        and (scored["M"] is not None)
    )
    exact_returned = (
        sum(
            (
                Fraction(g["weight_numerator"], g["weight_denominator"])
                for (g) in (selection.get("genes", []))
            ),
            Fraction(0),
        )
        if (scored["selection_valid"])
        else None
    )
    offered_identified = all_returned_assayed and exact_returned == scored["k"]
    scored["full_returned_mean"] = scored["M"] if (all_returned_assayed) else None
    scored["full_returned_mean_bound"] = (
        "POINT_IDENTIFIED"
        if (all_returned_assayed)
        else "UNDEFINED_EMPTY_RETURN"
        if (exact_returned == 0)
        else "NOT_IDENTIFIED_UNBOUNDED_WITHOUT_SCORE_BOUNDS"
        if (scored["selection_valid"])
        else "UNAVAILABLE_SOURCE_SELECTION"
    )
    scored["unconditional_offered_budget_utility"] = None
    if (offered_identified):
        scored["unconditional_offered_budget_utility"] = scored["M"]
        scored["unconditional_offered_budget_bound"] = "POINT_IDENTIFIED"
    elif (scored["selection_valid"] and exact_returned < scored["k"]):
        scored["unconditional_offered_budget_bound"] = (
            "NOT_DEFINED_WITHOUT_VACANCY_UTILITY_ASSUMPTION"
        )
    else:
        scored["unconditional_offered_budget_bound"] = (
            "NOT_IDENTIFIED_UNBOUNDED_WITHOUT_SCORE_BOUNDS"
            if (scored["selection_valid"])
            else "UNAVAILABLE_SOURCE_SELECTION"
        )
    return scored

def paired_offered_utility(left, right):
    left_value, right_value = (
        arm["unconditional_offered_budget_utility"] for (arm) in ((left, right))
    )
    point = left_value is not None and right_value is not None
    return dict(
        unconditional_delta_M = left_value - right_value if (point) else None,
        unconditional_delta_bound = "POINT_IDENTIFIED"
        if (point)
        else "NOT_IDENTIFIED_OR_NOT_ANALYSED",
        shared_unknown_coefficient_cancellation_analysed = False,
    )

def analyze(selections, endpoints, master_seeds = MASTERS):
    require(
        list(master_seeds) == list(MASTERS)
        and all((type(s) is int for (s) in (master_seeds))),
        "Exact five frozen PRODIGY seeds required",
    )
    support = validate_endpoints(endpoints)
    expected = [
        (model, method, pool, seed, policy, k)
        for ((model, method, pool, policy, k)) in (
            itertools.product(MODELS, METHODS, POOLS, POLICIES, BUDGETS)
        )
        for (seed) in (seed_list(method, master_seeds))
    ]
    allowed, provided = (set(expected), {})
    for (row) in (selections):
        require(
            isinstance(row, dict) and type(row.get("k")) is int,
            "Selection object and integer budget required",
        )
        require(
            type(row.get("master_seed")) is int
            if (row.get("method") == "PRODIGY")
            else row.get("master_seed") is None,
            "Exact seed type required",
        )
        actual = key(row)
        require(actual in allowed, "Selection outside frozen full grid")
        require(actual not in provided, "Duplicate selection group")
        provided[actual] = row
    metrics, indexed = ([], {})
    for (group) in (expected):
        row = provided.get(
            group,
            dict(
                zip(KEY_FIELDS, group),
                native_status = "UNAVAILABLE",
                status = "MISSING_GROUP",
                reason = "Expected selection row absent",
            ),
        )
        scored = _bounds(metric(row, endpoints.get(group[0], {})), row)
        scored.update(
            model_id = group[0],
            group_present_in_input = group in provided,
            assessed_genes = support[group[0]]["assessed_genes"],
            persona_pool_alias_not_independent = group[1] == "PersonaDrive",
        )
        metrics.append(scored)
        indexed[group] = scored
    seed_contrasts, model_contrasts, group_contrasts = ([], [], [])
    for (method, pool, k, pair) in (itertools.product(
        METHODS, POOLS, BUDGETS, PRIMARY_PAIRS + DEGREE_PAIRS
    )):
        left_policy, right_policy = pair
        primary = pair in PRIMARY_PAIRS and pool == "all34" and (k == 10)
        family = (
            "prespecified_primary_pair"
            if (pair in PRIMARY_PAIRS)
            else "degree_control_same_paired_support"
        )
        models = []
        for (model) in (MODELS):
            rows = []
            for (seed) in (seed_list(method, master_seeds)):
                left, right = (
                    indexed[model, method, pool, seed, p, k] for (p) in (pair)
                )
                r = paired_seed(left, right, *pair)
                r.update(
                    family = family,
                    primary = primary,
                    **paired_offered_utility(left, right),
                )
                rows.append(r)
                seed_contrasts.append(r)
            reduced = aggregate_seeds(rows, seed_list(method, master_seeds))
            reduced.update(family = family, primary = primary, model_id = model)
            offered = [r["unconditional_delta_M"] for (r) in (rows)]
            reduced.update(
                unconditional_delta_M = mean(offered)
                if (all((v is not None for (v) in (offered))))
                else None,
                unconditional_delta_bound = "POINT_IDENTIFIED"
                if (all((v is not None for (v) in (offered))))
                else "NOT_IDENTIFIED_OR_NOT_ANALYSED",
            )
            models.append(reduced)
            model_contrasts.append(reduced)
        grouped = full_model_summary(models, "delta_M")
        grouped.update(
            {
                field: models[0][field]
                for (field) in (
                    ("method", "normal_pool", "k", "left_policy", "right_policy")
                )
            }
        )
        grouped["equal_model_delta_mean"] = grouped.pop("equal_model_mean")
        grouped["model_delta_range"] = grouped.pop("model_range")
        offered = [r["unconditional_delta_M"] for (r) in (models)]
        grouped.update(
            family = family,
            primary = primary,
            unconditional_delta_M = mean(offered)
            if (all((v is not None for (v) in (offered))))
            else None,
            unconditional_delta_bound = "POINT_IDENTIFIED"
            if (all((v is not None for (v) in (offered))))
            else "NOT_IDENTIFIED_OR_NOT_ANALYSED",
        )
        group_contrasts.append(grouped)
    model_metrics, group_metrics = ([], [])
    for (method, pool, policy, k) in (itertools.product(METHODS, POOLS, POLICIES, BUDGETS)):
        model_rows = []
        for (model) in (MODELS):
            rows = [
                indexed[model, method, pool, seed, policy, k]
                for (seed) in (seed_list(method, master_seeds))
            ]
            full = all((r["M"] is not None for (r) in (rows)))
            reduced = dict(
                sample_id = model,
                model_id = model,
                method = method,
                normal_pool = pool,
                policy = policy,
                k = k,
                expected_seeds = seed_list(method, master_seeds),
                valid_seed_count = sum((r["M"] is not None for (r) in (rows))),
                status = "IDENTIFIED_COMPLETE_SEED_SUPPORT"
                if (full)
                else "UNIDENTIFIED_INCOMPLETE_SEED_SUPPORT",
                missing_seed_ids = [
                    r["master_seed"] for (r) in (rows) if (r["M"] is None)
                ],
                M = mean([r["M"] for (r) in (rows)]) if (full) else None,
                seed_results = rows,
            )
            for (field) in (("R", "A", "returned_coverage", "assessed_coverage")):
                reduced[field] = (
                    mean([r[field] for (r) in (rows)])
                    if (all((r[field] is not None for (r) in (rows))))
                    else None
                )
            model_rows.append(reduced)
            model_metrics.append(reduced)
        grouped = full_model_summary(model_rows, "M")
        grouped.update(
            method = method,
            normal_pool = pool,
            policy = policy,
            k = k,
            persona_pool_alias_not_independent = method == "PersonaDrive",
        )
        for (field) in (("R", "A", "returned_coverage", "assessed_coverage")):
            grouped[field + "_equal_model_mean"] = (
                mean([r[field] for (r) in (model_rows)])
                if (all((r[field] is not None for (r) in (model_rows))))
                else None
            )
        group_metrics.append(grouped)
    return dict(
        status = "DESCRIPTIVE_EIGHT_MODEL_EXTERNAL_VALIDATION_NOT_INDEPENDENT_DONOR_CONFIRMATION",
        model_order = list(MODELS),
        master_seeds = list(MASTERS),
        expected_selection_groups = len(expected),
        provided_selection_groups = len(provided),
        missing_selection_groups = len(expected) - len(provided),
        metric_status_counts = dict(
            collections.Counter((r["status"] for (r) in (metrics)))
        ),
        endpoint_support = support,
        metrics = metrics,
        paired_seed_contrasts = seed_contrasts,
        paired_model_contrasts = model_contrasts,
        group_contrasts = group_contrasts,
        primary_group_contrasts = [r for (r) in (group_contrasts) if (r["primary"])],
        model_policy_metrics = model_metrics,
        group_policy_metrics = group_metrics,
        interpretation = dict(
            M = "Negative weighted Chronos conditional mean among assayed returned genes; A=0 is unidentified",
            coverage = ["R/k returned", "A/k assayed"],
            unknown_scores = "Absent/null stay unknown; finite zero is assayed",
            unconditional_estimand = "Point identified for full assayed offered support; returned mean identified if all returns assayed; missing score bounds absent and vacancy utility unspecified; shared-unknown cancellation not analysed",
            seed_aggregation = "All exact five PRODIGY seeds must be identified; deterministic methods one null seed",
            model_aggregation = "Equal-eight-model mean only if every model identified; available-only descriptions separate",
            inference = "No p-values, bootstrap intervals, independent-donor or population-generalisation claim",
            adverse_results = "All negative/null/unavailable results retained",
            prior_three_models = "Closed and separate; never pooled into eleven",
            persona_pools = "Logical aliases, not independent observations",
            source_limitations = [
                "Availability selected eight models",
                "Mixed GDC STAR tumour and deposited normal count processing",
                "ACH-003346 RNA85A versus MAF85B",
                "Same model name does not prove identical aliquot, passage, or donor independence",
            ],
        ),
    )

def normalize(rows):
    output = []
    status_map = {
        "available_nonempty": "AVAILABLE",
        "available_empty": "AVAILABLE",
        "unavailable": "UNAVAILABLE",
    }
    for (raw) in (rows):
        require(isinstance(raw, dict), "Source selection must be object")
        require(
            raw.get("sample_id") == raw.get("model_id")
            and raw.get("sample_id") in MODELS,
            "Exact model/sample identity required",
        )
        require(raw.get("support") == "native", "Native support required")
        require(raw.get("status") in status_map, "Unrecognised source selection status")
        require(isinstance(raw.get("genes"), list), "Explicit genes list required")
        row = dict(
            raw,
            raw_policy = raw["policy"],
            raw_status = raw["status"],
            raw_native_status = raw["native_status"],
        )
        row["policy"] = "V0DEGREE" if (raw["policy"] == "V0_DEGREE") else raw["policy"]
        row["status"] = status_map[raw["status"]]
        if ("returned_numerator" in raw or "returned_denominator" in raw):
            require(
                "returned_numerator" in raw
                and "returned_denominator" in raw
                and (raw["returned_numerator"] == raw.get("returned_mass_numerator"))
                and (
                    raw["returned_denominator"] == raw.get("returned_mass_denominator")
                ),
                "Legacy returned fraction aliases disagree",
            )
        if (raw["status"] == "unavailable"):
            require(
                raw["genes"] == []
                and all(
                    (
                        field in raw and raw[field] is None
                        for (field) in (
                            (
                                "returned_mass",
                                "returned_mass_numerator",
                                "returned_mass_denominator",
                            )
                        )
                    )
                ),
                "Unavailable selection mass must remain explicit null, not empty zero",
            )
            output.append(row)
            continue
        mass = Fraction(0)
        for (gene) in (raw["genes"]):
            n, d = (gene.get("weight_numerator"), gene.get("weight_denominator"))
            require(
                type(n) is int and type(d) is int and (d > 0),
                "Exact integer weight fraction required",
            )
            exact = Fraction(n, d)
            require(
                0 < exact <= 1
                and number(gene.get("weight"))
                and (float(exact) == gene["weight"]),
                "Selection weight/fraction changed",
            )
            mass += exact
        require(
            type(raw.get("returned_mass_numerator")) is int
            and type(raw.get("returned_mass_denominator")) is int
            and (raw["returned_mass_denominator"] > 0),
            "Returned exact fraction required",
        )
        require(
            mass
            == Fraction(
                raw["returned_mass_numerator"], raw["returned_mass_denominator"]
            )
            and number(raw.get("returned_mass"))
            and (float(mass) == raw["returned_mass"]),
            "Returned mass disagrees with full genes",
        )
        require(
            raw["status"] != "available_empty" or not raw["genes"],
            "Empty status has genes",
        )
        require(
            raw["status"] != "available_nonempty" or bool(raw["genes"]),
            "Nonempty status lacks genes",
        )
        output.append(row)
    return output
