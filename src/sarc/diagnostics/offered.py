import math

from statistics import mean

SEEDS = [104729, 130363, 155921, 196613, 228017]

METHODS = ["PRODIGY", "DawnRank", "PersonaDrive"]

POLICIES = ["FIXED_V0", "FIXED_V1", "FIXED_V2", "RC", "STABLE_ONLY_3", "V0_DEGREE"]

MODELS = {
    "eight": [
        "ACH-003332",
        "ACH-003334",
        "ACH-003337",
        "ACH-003346",
        "ACH-003387",
        "ACH-003400",
        "ACH-003401",
        "ACH-003402",
    ],
    "three": ["HCM-SANG-0265-C18", "HCM-SANG-0270-C20", "HCM-SANG-0286-C20"],
}

def score_total(record):
    r, a, m = record.get("R"), record.get("A"), record.get("M")
    if (not record.get("selection_valid") or record.get("source_status") != "AVAILABLE"):
        return None, "unavailable_prediction"
    if (r is None or a is None or not math.isfinite(r) or not math.isfinite(a)):
        return None, "unavailable_coverage"
    if (record.get("unassessed_genes") or r != a):
        return None, "unassessed_returned_mass"
    if (r == 0):
        if (a != 0 or m is not None):
            raise ValueError("invalid_empty_state")
        return 0.0, "valid_empty_sum"
    if (m is None or not math.isfinite(m)):
        return None, "unavailable_conditional_score"
    return a * m / record["k"], "complete_assessment"

def required_seeds(method):
    return SEEDS if (method == "PRODIGY") else [None]

def complete_mean(values):
    return mean(values) if (all(v is not None for (v) in (values))) else None

def calculate(name, data):
    selected = [x for (x) in (data["metrics"]) if (x["normal_pool"] == "all34" and x["k"] == 10)]
    expected = {
        (m, p, s, d)
        for (m) in (METHODS)
        for (p) in (POLICIES)
        for (s) in (required_seeds(m))
        for (d) in (MODELS[name])
    }
    units, lookup = [], {}
    for (record) in (selected):
        policy = "V0_DEGREE" if (record["policy"] == "V0DEGREE") else record["policy"]
        key = record["method"], policy, record["master_seed"], record["sample_id"]
        if (key in lookup or key not in expected):
            raise ValueError(("duplicate_or_unexpected_unit", key))
        q, state = score_total(record)
        unit = dict(
            method = key[0],
            policy = key[1],
            seed = key[2],
            model = key[3],
            R = record["R"],
            A = record["A"],
            M = record["M"],
            Q = q,
            S = None if (q is None) else q * 10,
            state = state,
            original_metric_status = record["status"],
            unassessed_genes = record.get("unassessed_genes", []),
        )
        if (q is not None and record["R"] > 0):
            assert math.isclose(q, record["R"] * record["M"] / 10, rel_tol = 1e-14, abs_tol = 1e-14)
        lookup[key] = unit
        units.append(unit)
    assert set(lookup) == expected
    model_rows, group_rows, contrast_rows = [], [], []
    for (method) in (METHODS):
        for (policy) in (POLICIES):
            rows = []
            for (model) in (MODELS[name]):
                items = [lookup[method, policy, seed, model] for (seed) in (required_seeds(method))]
                row = dict(
                    method = method,
                    policy = policy,
                    model = model,
                    Q = complete_mean([x["Q"] for (x) in (items)]),
                    M = complete_mean([x["M"] for (x) in (items)]),
                    R = mean(x["R"] for (x) in (items)),
                    A = mean(x["A"] for (x) in (items)),
                    empty_seeds = [x["seed"] for (x) in (items) if (x["state"] == "valid_empty_sum")],
                    unavailable_seeds = [x["seed"] for (x) in (items) if (x["Q"] is None)],
                )
                rows.append(row)
                model_rows.append(row)
            qmean = complete_mean([x["Q"] for (x) in (rows)])
            group_rows.append(
                dict(
                    method = method,
                    policy = policy,
                    Q = qmean,
                    B = None if (qmean is None) else 10 * qmean,
                    M = complete_mean([x["M"] for (x) in (rows)]),
                    R = mean(x["R"] for (x) in (rows)),
                    A = mean(x["A"] for (x) in (rows)),
                    complete_Q_models = sum(x["Q"] is not None for (x) in (rows)),
                    expected_models = len(rows),
                    unavailable_Q_models = [x["model"] for (x) in (rows) if (x["Q"] is None)],
                )
            )
        pairs = [(p, "V0_DEGREE") for (p) in (POLICIES) if (p != "V0_DEGREE")] + [
            ("FIXED_V2", "FIXED_V0"),
            ("RC", "FIXED_V2"),
            ("STABLE_ONLY_3", "FIXED_V2"),
        ]
        for (left, right) in (pairs):
            values = []
            for (model) in (MODELS[name]):
                differences = []
                for (seed) in (required_seeds(method)):
                    l, r = lookup[method, left, seed, model], lookup[method, right, seed, model]
                    differences.append(
                        None if (l["Q"] is None or r["Q"] is None) else l["Q"] - r["Q"]
                    )
                values.append(dict(model = model, delta_Q = complete_mean(differences)))
            delta = complete_mean([x["delta_Q"] for (x) in (values)])
            if (delta is not None):
                gl = next(
                    x for (x) in (group_rows) if (x["method"] == method and x["policy"] == left)
                )
                gr = next(
                    x for (x) in (group_rows) if (x["method"] == method and x["policy"] == right)
                )
                assert math.isclose(delta, gl["Q"] - gr["Q"], rel_tol = 1e-12, abs_tol = 1e-12)
            contrast_rows.append(
                dict(
                    method = method,
                    left_policy = left,
                    right_policy = right,
                    delta_Q = delta,
                    complete_models = sum(x["delta_Q"] is not None for (x) in (values)),
                    expected_models = len(values),
                    model_results = values,
                )
            )
    return dict(
        model_order = MODELS[name],
        units = units,
        model_policy = model_rows,
        group_policy = group_rows,
        contrasts = contrast_rows,
        interpretation = "Exploratory measured-score total per offered slot; distinct from conditional M and clinical utility",
    )
