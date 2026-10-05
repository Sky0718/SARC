import collections
import itertools
import json
import math

SAMPLE_MODELS = {
    "HCM-SANG-0265-C18": "ACH-003388",
    "HCM-SANG-0270-C20": "ACH-003390",
    "HCM-SANG-0286-C20": "ACH-003391",
}
METHODS = ("PRODIGY", "DawnRank", "PersonaDrive")
POOLS = ("all34", "distal16", "proximal18")
BUDGETS = (1, 5, 10, 20)
POLICIES = ("FIXED_V0", "FIXED_V1", "FIXED_V2", "RC", "STABLE_ONLY_3", "V0DEGREE")
PRIMARY_PAIRS = (
    ("FIXED_V2", "FIXED_V0"),
    ("RC", "FIXED_V2"),
    ("STABLE_ONLY_3", "FIXED_V2"),
)
DEGREE_PAIRS = tuple(((p, "V0DEGREE") for (p) in (POLICIES) if (p != "V0DEGREE")))
KEY_FIELDS = ("sample_id", "method", "normal_pool", "master_seed", "policy", "k")

def require(value, message):
    if (not value):
        raise ValueError(message)

def number(value):
    return type(value) in (int, float) and math.isfinite(value)

def mean(values):
    return math.fsum(values) / len(values) if (values) else None

def key(row):
    return tuple((row.get(field) for (field) in (KEY_FIELDS)))

def seed_list(method, master_seeds):
    return list(master_seeds) if (method == "PRODIGY") else [None]

def validate_endpoints(endpoints):
    require(
        isinstance(endpoints, dict) and set(endpoints).issubset(SAMPLE_MODELS),
        "Endpoints restricted to exact3 sample identifiers",
    )
    support = {}
    for (sample) in (SAMPLE_MODELS):
        scores = endpoints.get(sample, {})
        require(
            isinstance(scores, dict),
            "Each sample endpoint must be gene->finite number ornull",
        )
        require(
            all((isinstance(g, str) and g and (g == g.strip()) for (g) in (scores))),
            "Exact canonical gene keys required",
        )
        require(
            all((value is None or number(value) for (value) in (scores.values()))),
            "Nonfinite/string/Boolean Chronos value rejected; no silent replacement",
        )
        finite = sorted((g for ((g, value)) in (scores.items()) if (value is not None)))
        support[sample] = {
            "sample_endpoint_record_present": sample in endpoints,
            "assessed_gene_count": len(finite),
            "finite_zero_count": sum(
                (value == 0 for (value) in (scores.values()) if (value is not None))
            ),
            "assessed_genes": finite,
            "score_threshold_used": False,
        }
    return support

def metric(selection, scores):
    source_native = selection.get("native_status")
    source_status = selection.get("status")
    result = {field: selection.get(field) for (field) in (KEY_FIELDS)}
    result.update(
        source_native_status = source_native,
        source_status = source_status,
        source_reason = selection.get("reason"),
        R = None,
        A = None,
        M = None,
        returned_coverage = None,
        assessed_coverage = None,
        selection_valid = False,
    )
    admitted_source = (
        source_native == "SUCCESS"
        or (source_native == "SUCCESS_EMPTY" and selection.get("genes") == [])
        or (
            source_native == "SOURCE_BOUND_UNTRAINED_V0_DEGREE"
            and selection.get("policy") == "V0DEGREE"
        )
    )
    if (not admitted_source or source_status != "AVAILABLE"):
        result.update(
            status = "UNAVAILABLE_SOURCE_SELECTION",
            reason = selection.get("reason")
            or f"native_status={source_native};status={source_status}",
        )
        return result
    genes = selection.get("genes")
    if (not isinstance(genes, list)):
        return dict(
            result, status = "INVALID_SELECTION", reason = "genes must be explicit list"
        )
    seen = set()
    weights = []
    assessed = []
    weighted = []
    missing = []
    for (item) in (genes):
        if (not isinstance(item, dict)):
            return dict(
                result, status = "INVALID_SELECTION", reason = "non-object gene entry"
            )
        gene = item.get("gene")
        weight = item.get("weight")
        if (
            not isinstance(gene, str)
            or not gene
            or gene != gene.strip()
            or (gene in seen)
        ):
            return dict(
                result,
                status = "INVALID_SELECTION",
                reason = "empty/noncanonical/duplicate gene identifier",
            )
        if (not number(weight) or not 0 < weight <= 1):
            return dict(
                result,
                status = "INVALID_SELECTION",
                reason = "weight must be finite in(0,1]",
            )
        seen.add(gene)
        weights.append(weight)
        if (gene not in scores):
            missing.append(
                {"gene": gene, "weight": weight, "availability": "ABSENT_ENDPOINT_GENE"}
            )
        elif (scores[gene] is None):
            missing.append(
                {
                    "gene": gene,
                    "weight": weight,
                    "availability": "EXPLICIT_NULL_ENDPOINT",
                }
            )
        else:
            assessed.append(weight)
            weighted.append(weight * scores[gene])
    returned = math.fsum(weights)
    offered = selection["k"]
    if (returned > offered + 1e-09):
        return dict(
            result,
            status = "INVALID_SELECTION",
            reason = "returned fractional mass exceeds fixed offered budget",
        )
    assayed = math.fsum(assessed)
    effect = (
        -math.fsum((value / assayed for (value) in (weighted)))
        if (assayed > 0)
        else None
    )
    result.update(
        R = returned,
        A = assayed,
        M = effect,
        returned_coverage = returned / offered,
        assessed_coverage = assayed / offered,
        selection_valid = True,
        selected_gene_count = len(genes),
        unassessed_genes = missing,
        status = "AVAILABLE_CONDITIONAL_ASSESSED_MEAN"
        if (assayed > 0)
        else "UNIDENTIFIED_ZERO_ASSESSED_MASS",
        reason = None if (assayed > 0) else "A=0; unavailable mean, notzero",
    )
    return result

def paired_seed(left, right, left_policy, right_policy):
    common = {field: left[field] for (field) in (KEY_FIELDS) if (field != "policy")}
    valid = left["M"] is not None and right["M"] is not None
    return dict(
        common,
        left_policy = left_policy,
        right_policy = right_policy,
        status = "PAIRED_AVAILABLE" if (valid) else "PAIRED_UNIDENTIFIED",
        delta_M = left["M"] - right["M"] if (valid) else None,
        left_M = left["M"],
        right_M = right["M"],
        left_R = left["R"],
        right_R = right["R"],
        left_A = left["A"],
        right_A = right["A"],
        left_returned_coverage = left["returned_coverage"],
        right_returned_coverage = right["returned_coverage"],
        left_assessed_coverage = left["assessed_coverage"],
        right_assessed_coverage = right["assessed_coverage"],
        left_status = left["status"],
        right_status = right["status"],
        left_reason = left.get("reason"),
        right_reason = right.get("reason"),
        same_sample_seed_endpoint_universe = True,
        complete_case_replacement_of_estimand = False,
    )

def aggregate_seeds(rows, expected_seeds):
    require(
        len(rows) == len(expected_seeds)
        and [r["master_seed"] for (r) in (rows)] == list(expected_seeds),
        "Exact full paired seed list/order required",
    )
    available = [r["delta_M"] for (r) in (rows) if (r["delta_M"] is not None)]
    full = len(available) == len(expected_seeds)
    first = rows[0]
    result = {
        field: first[field]
        for (field) in (
            ("sample_id", "method", "normal_pool", "k", "left_policy", "right_policy")
        )
    }
    result.update(
        expected_seeds = list(expected_seeds),
        expected_seed_count = len(expected_seeds),
        valid_paired_seed_count = len(available),
        status = "IDENTIFIED_COMPLETE_SEED_SUPPORT"
        if (full)
        else "UNIDENTIFIED_INCOMPLETE_SEED_SUPPORT",
        delta_M = mean(available) if (full) else None,
        descriptive_available_seed_mean = mean(available),
        descriptive_available_seed_range = [min(available), max(available)]
        if (available)
        else None,
        equal_weight_per_seed = True,
        missing_seed_ids = [r["master_seed"] for (r) in (rows) if (r["delta_M"] is None)],
        unavailable_states = [
            {
                "master_seed": r["master_seed"],
                "left_status": r["left_status"],
                "right_status": r["right_status"],
                "left_reason": r["left_reason"],
                "right_reason": r["right_reason"],
            }
            for (r) in (rows)
            if (r["delta_M"] is None)
        ],
    )
    for (field) in ((
        "left_M",
        "right_M",
        "left_R",
        "right_R",
        "left_A",
        "right_A",
        "left_returned_coverage",
        "right_returned_coverage",
        "left_assessed_coverage",
        "right_assessed_coverage",
    )):
        values = [r[field] for (r) in (rows)]
        result[field + "_equal_seed_mean"] = (
            mean(values) if (all((v is not None for (v) in (values)))) else None
        )
    return result

def aggregate_models(rows):
    require(
        [r["sample_id"] for (r) in (rows)] == list(SAMPLE_MODELS),
        "Exact3-model order/support required",
    )
    complete = [r for (r) in (rows) if (r["delta_M"] is not None)]
    identified = len(complete) == 3
    values = [r["delta_M"] for (r) in (complete)]
    first = rows[0]
    result = {
        field: first[field]
        for (field) in (("method", "normal_pool", "k", "left_policy", "right_policy"))
    }
    result.update(
        status = "IDENTIFIED_DESCRIPTIVE_N3"
        if (identified)
        else "UNIDENTIFIED_INCOMPLETE_THREE_MODEL_SUPPORT",
        equal_model_delta_mean = mean(values) if (identified) else None,
        model_delta_range = [min(values), max(values)] if (identified) else None,
        expected_model_count = 3,
        valid_complete_seed_model_count = len(complete),
        descriptive_available_model_mean = mean(values),
        descriptive_available_model_ids = [r["sample_id"] for (r) in (complete)],
        unavailable_model_ids = [
            r["sample_id"] for (r) in (rows) if (r["delta_M"] is None)
        ],
        model_results = rows,
        equal_weight_per_model = True,
        inferential_p_value = None,
        confidence_interval = None,
        independent_confirmation_claim = False,
    )
    result["leave_one_model_out_means"] = (
        [
            {
                "omitted_sample_id": omitted["sample_id"],
                "remaining_sample_ids": [
                    r["sample_id"] for (r) in (rows) if (r is not omitted)
                ],
                "mean": mean([r["delta_M"] for (r) in (rows) if (r is not omitted)]),
            }
            for (omitted) in (rows)
        ]
        if (identified)
        else None
    )
    return result

def analyze(selections, endpoints, master_seeds):
    require(
        isinstance(master_seeds, (list, tuple))
        and len(master_seeds) == 5
        and (len(set(master_seeds)) == 5)
        and all((type(x) is int and x > 0 for (x) in (master_seeds))),
        "Exactly five frozen PRODIGY master seeds required",
    )
    endpoint_support = validate_endpoints(endpoints)
    expected = []
    for (sample, method, pool, policy, k) in (itertools.product(
        SAMPLE_MODELS, METHODS, POOLS, POLICIES, BUDGETS
    )):
        for (seed) in (seed_list(method, master_seeds)):
            expected.append((sample, method, pool, seed, policy, k))
    allowed = set(expected)
    provided = {}
    for (row) in (selections):
        require(isinstance(row, dict), "Selection rows must be objects")
        require(type(row.get("k")) is int, "Offered budget must be integer, notBoolean")
        require(
            type(row.get("master_seed")) is int
            if (row.get("method") == "PRODIGY")
            else row.get("master_seed") is None,
            "ExplicitPRODIGY integer seed or deterministic null required",
        )
        actual = key(row)
        require(
            actual in allowed,
            "Selection outside exact sample/method/pool/seed/policy/budget grid",
        )
        require(actual not in provided, "Duplicate selection group; no last-write-wins")
        provided[actual] = row
    metrics = []
    indexed = {}
    for (group) in (expected):
        row = provided.get(group)
        if (row is None):
            row = dict(
                zip(KEY_FIELDS, group),
                native_status = "UNAVAILABLE",
                status = "MISSING_GROUP",
                reason = "Expected selection row absent",
            )
        scored = metric(row, endpoints.get(group[0], {}))
        scored["group_present_in_input"] = group in provided
        scored["model_id"] = SAMPLE_MODELS[group[0]]
        scored["assessed_genes"] = endpoint_support[group[0]]["assessed_genes"]
        metrics.append(scored)
        indexed[group] = scored
    seed_contrasts = []
    model_contrasts = []
    group_contrasts = []
    for (method, pool, k, pair) in (itertools.product(
        METHODS, POOLS, BUDGETS, PRIMARY_PAIRS + DEGREE_PAIRS
    )):
        left_policy, right_policy = pair
        models = []
        family = (
            "prespecified_primary_pair"
            if (pair in PRIMARY_PAIRS)
            else "degree_control_same_paired_support"
        )
        primary = pair in PRIMARY_PAIRS and pool == "all34" and (k == 10)
        for (sample) in (SAMPLE_MODELS):
            rows = []
            for (seed) in (seed_list(method, master_seeds)):
                left = indexed[sample, method, pool, seed, left_policy, k]
                right = indexed[sample, method, pool, seed, right_policy, k]
                r = paired_seed(left, right, left_policy, right_policy)
                r.update(family = family, primary = primary)
                rows.append(r)
                seed_contrasts.append(r)
            model = aggregate_seeds(rows, seed_list(method, master_seeds))
            model.update(family = family, primary = primary, model_id = SAMPLE_MODELS[sample])
            models.append(model)
            model_contrasts.append(model)
        group = aggregate_models(models)
        group.update(family = family, primary = primary)
        group_contrasts.append(group)
    return {
        "status": "DESCRIPTIVE_THREE_MODEL_CROSS_EXPERIMENT_ANALYSIS_NOT_INDEPENDENT_CONFIRMATION",
        "sample_model_mapping": SAMPLE_MODELS.copy(),
        "master_seeds": list(master_seeds),
        "expected_selection_groups": len(expected),
        "provided_selection_groups": len(provided),
        "missing_selection_groups": len(expected) - len(provided),
        "metric_status_counts": dict(
            collections.Counter((r["status"] for (r) in (metrics)))
        ),
        "endpoint_support": endpoint_support,
        "metrics": metrics,
        "paired_seed_contrasts": seed_contrasts,
        "paired_model_contrasts": model_contrasts,
        "group_contrasts": group_contrasts,
        "primary_group_contrasts": [r for (r) in (group_contrasts) if (r["primary"])],
        "interpretation": {
            "M": "Minus weighted Chronos mean conditional on assessed selected gene mass; M isunidentified whenA=0",
            "coverage": ["R/k returned mass", "A/k assessed mass"],
            "unknown_scores": "Absent/null scores stayunavailable, finitezero isassessed; no imputation",
            "unassessed_selection_mean_bound": "Unbounded without additional score bounds; no full-selection effect claim",
            "paired_support": "Same fixed model,seed,pool,budget and endpoint gene universe; no endpoint-sign/essentiality filter",
            "seed_aggregation": "Equal mean ofpaired deltas onlyifall5PRODIGY seeds valid; deterministicothers exactlyone",
            "model_aggregation": "Equal mean onlyifall3models identified; available-only descriptions clearlyseparate",
            "inference": "Descriptive n=3 only; no p-values,confidence intervals or donor-independence claim",
            "adverse_results": "Negative andzeroeffects are preserved; no tuning,thresholding or baseline weakening",
        },
    }
