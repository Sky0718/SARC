import math
from fractions import Fraction
import numpy as np

DEVELOPMENT = ("COAD_CCLE", "LUAD_CCLE")
METHODS = ("DawnRank", "PersonaDrive", "PRODIGY")
MASTERS = (104729, 130363, 155921, 196613, 228017)
FIT_CONTROLS = (271828, 314159)
HOLDOUT = 173205
FEATURES = (
    "log_edit_count",
    "mutated_edit_fraction",
    "deg_edit_fraction",
    "between_module_edit_fraction",
    "mean_endpoint_log_degree",
    "mean_consumed_cost_change",
)
MODELS = {
    "strict_input_only": FEATURES,
    "edit_volume_only": FEATURES[:1],
    "degree_incidence_only": (FEATURES[1], FEATURES[4]),
    "additive_volume_degree": (FEATURES[0], FEATURES[1], FEATURES[4]),
    "old_output_margin": ("old_native_rank_margin",),
    "input_plus_margin": FEATURES + ("old_native_rank_margin",),
}
PENALTIES = (0.1, 1.0, 10.0)

def expit(value):
    value = np.asarray(value, dtype = float)
    result = np.empty_like(value)
    nonnegative = value >= 0
    result[nonnegative] = 1 / (1 + np.exp(-value[nonnegative]))
    exponential = np.exp(value[~nonnegative])
    result[~nonnegative] = exponential / (1 + exponential)
    return result

def exact(value):
    if (isinstance(value, Fraction)):
        return value
    if (isinstance(value, dict) and set(value) >= {"numerator", "denominator"}):
        return Fraction(value["numerator"], value["denominator"])
    if (type(value) is int):
        return Fraction(value)
    raise ValueError(
        "Endpoint event requires exact inherited rational values, not rounded floats"
    )

def primary_event(method, old_hits, new_hits):
    expected = list(MASTERS) if (method == "PRODIGY") else [None]
    if (list(old_hits) != expected or list(new_hits) != expected):
        raise ValueError("Complete ordered paired algorithm seed states required")
    if (any((old_hits[seed] is None or new_hits[seed] is None for (seed) in (expected)))):
        return {
            "status": "UNAVAILABLE",
            "event": None,
            "signed_mean_delta": None,
            "per_seed_event_frequency_secondary": None,
            "expected_seed_count": len(expected),
        }
    differences = [
        exact(new_hits[seed]) - exact(old_hits[seed]) for (seed) in (expected)
    ]
    mean_delta = sum(differences, Fraction(0)) / len(expected)
    return {
        "status": "AVAILABLE",
        "event": int(abs(mean_delta) >= 1),
        "signed_mean_delta": mean_delta,
        "absolute_mean_delta": abs(mean_delta),
        "per_seed_event_frequency_secondary": Fraction(
            sum((abs(value) >= 1 for (value) in (differences))), len(expected)
        ),
        "seed_differences": differences,
        "expected_seed_count": len(expected),
        "seeds_are_biological_replicates": False,
    }

def old_native_margin(genes, scores, eligible, native_score_meaningful):
    if (not native_score_meaningful):
        return None
    if (
        len(genes) != len(scores)
        or len(genes) != len(set(genes))
        or any((not math.isfinite(float(value)) for (value) in (scores)))
        or any((float(a) < float(b) for ((a, b)) in (zip(scores, scores[1:]))))
    ):
        raise ValueError("Native ranked scores are invalid")
    allowed = set(eligible)
    native = [
        float(score) for ((gene, score)) in (zip(genes, scores)) if (gene in allowed)
    ]
    return None if (len(native) < 11) else native[9] - native[10]

def allowed_development(rows):
    identities = []
    for (row) in (rows):
        identities.append(row["row_id"])
        if (
            row["cohort"] not in DEVELOPMENT
            or row["method"] not in METHODS
            or row["control_seed"] not in (None, *FIT_CONTROLS)
        ):
            raise ValueError(
                "Protected, graph-holdout or out-of-scope row in development fitting"
            )
        expected = list(MASTERS) if (row["method"] == "PRODIGY") else [None]
        if (row["master_seeds"] != expected):
            raise ValueError(
                "Algorithm seeds were reduced or expanded into pseudo-people"
            )
        if (row["target"] not in (None, 0, 1)):
            raise ValueError(
                "Primary method/sample/intervention target must be binary or unavailable"
            )
        if (row["availability"] not in ("AVAILABLE", "UNAVAILABLE") or (
            row["availability"] == "AVAILABLE"
        ) != (row["target"] is not None)):
            raise ValueError("Availability must not be encoded as a negative endpoint")
    if (len(identities) != len(set(identities))):
        raise ValueError("Duplicated logical development row")

def columns(model):
    if (model not in MODELS):
        raise ValueError("Unfrozen predictor family")
    return list(MODELS[model]) + ["method_DawnRank", "method_PersonaDrive"]

def vector(row, model):
    result = []
    for (feature) in (MODELS[model]):
        value = row["features"].get(feature)
        if (feature == "mean_consumed_cost_change" and row["method"] in (
            "DawnRank",
            "PersonaDrive",
        )):
            if (
                value is not None
                or row["features"].get("cost_structurally_not_applicable") is not True
            ):
                raise ValueError(
                    "Binary method must explicitly retain structurally inapplicable costs"
                )
            value = 0.0
        if (value is None):
            return None
        value = (
            float(Fraction(value["numerator"], value["denominator"]))
            if (isinstance(value, dict))
            else float(value)
        )
        if (not math.isfinite(value)):
            raise ValueError("Nonfinite predictor input")
        result.append(value)
    return result + [
        float(row["method"] == "DawnRank"),
        float(row["method"] == "PersonaDrive"),
    ]

def available(rows, model):
    (retained, matrix, omitted) = ([], [], [])
    for (row) in (rows):
        values = vector(row, model)
        if (row["availability"] == "AVAILABLE" and values is not None):
            retained.append(row)
            matrix.append(values)
        else:
            omitted.append(
                {
                    "row_id": row["row_id"],
                    "reason": "OUTCOME_UNAVAILABLE"
                    if (row["availability"] != "AVAILABLE")
                    else "MODEL_INPUT_UNAVAILABLE",
                }
            )
    return (
        retained,
        np.asarray(matrix, dtype = float).reshape((len(matrix), len(columns(model)))),
        omitted,
    )

def ridge_fit(x, y, penalty, continuous_count):
    x = np.asarray(x, dtype = float)
    y = np.asarray(y, dtype = float)
    if (
        x.ndim != 2
        or len(x) != len(y)
        or (not len(y))
        or (not np.isin(y, [0, 1]).all())
        or (not np.isfinite(x).all())
        or (penalty not in PENALTIES)
    ):
        raise ValueError("Malformed complete development fitting matrix")
    centre = np.zeros(x.shape[1])
    scale = np.ones(x.shape[1])
    centre[:continuous_count] = x[:, :continuous_count].mean(axis = 0)
    scale[:continuous_count] = x[:, :continuous_count].std(axis = 0, ddof = 0)
    scale[scale == 0] = 1.0
    matrix = np.column_stack((np.ones(len(x)), (x - centre) / scale))
    if (np.all(y == y[0])):
        return {
            "status": "CONSTANT_CLASS",
            "constant_probability": float(y[0]),
            "penalty": penalty,
            "centre": centre.tolist(),
            "scale": scale.tolist(),
            "intercept": None,
            "coefficients": None,
            "rows": len(y),
        }

    def objective(theta):
        linear = matrix @ theta
        value = np.mean(np.logaddexp(0, linear) - y * linear) + 0.5 * penalty * np.dot(
            theta[1:], theta[1:]
        )
        gradient = matrix.T @ (expit(linear) - y) / len(y)
        gradient[1:] += penalty * theta[1:]
        return (value, gradient)

    initial = np.zeros(matrix.shape[1])
    initial[0] = math.log(float(y.mean()) / (1 - float(y.mean())))
    theta = initial
    converged = False
    regularisation = np.diag([0.0] + [penalty] * (matrix.shape[1] - 1))
    for (iteration) in (range(100)):
        (value, gradient) = objective(theta)
        if (np.max(np.abs(gradient)) <= 1e-09):
            converged = True
            break
        probabilities = expit(matrix @ theta)
        hessian = (
            matrix.T
            @ (matrix * (probabilities * (1 - probabilities))[:, None])
            / len(y)
            + regularisation
        )
        direction = np.linalg.solve(hessian, gradient)
        step = 1.0
        for (backtrack) in (range(60)):
            proposed = theta - step * direction
            if (objective(proposed)[0] <= value - 0.0001 * step * np.dot(
                gradient, direction
            )):
                theta = proposed
                break
            step *= 0.5
        else:
            raise ValueError("Ridge Newton line search failed")
    (value, gradient) = objective(theta)
    if (not np.isfinite(theta).all() or np.max(np.abs(gradient)) > 1e-07):
        raise ValueError(
            "Ridge numerical convergence failed; never accept failed coefficients"
        )
    return {
        "status": "FITTED",
        "penalty": penalty,
        "centre": centre.tolist(),
        "scale": scale.tolist(),
        "intercept": float(theta[0]),
        "coefficients": theta[1:].tolist(),
        "rows": len(y),
        "objective": float(value),
        "gradient_max_abs": float(np.max(np.abs(gradient))),
        "solver_reported_success": converged,
        "solver": "Newton exact gradient/Hessian with Armijo backtracking",
        "iterations": iteration + 1,
    }

def predict(model, x):
    x = np.asarray(x, dtype = float)
    if (x.ndim != 2 or x.shape[1] != len(model["centre"])):
        raise ValueError("Predictor feature axes differ")
    if (model["status"] == "CONSTANT_CLASS"):
        return np.full(len(x), model["constant_probability"])
    return expit(
        model["intercept"]
        + (x - model["centre"]) / model["scale"] @ np.asarray(model["coefficients"])
    )

def select_penalty(trials):
    if ([trial["penalty"] for (trial) in (trials)] != list(PENALTIES) or any(
        (not math.isfinite(trial["mean_heldout_brier"]) for (trial) in (trials))
    )):
        raise ValueError("Complete prespecified penalty grid required")
    return min(
        trials, key = lambda trial: (trial["mean_heldout_brier"], -trial["penalty"])
    )["penalty"]

def diagnostics(rows, probabilities):
    if (len(rows) != len(probabilities) or not rows):
        raise ValueError("Complete nonempty paired prediction carrier required")
    y = np.asarray([row["target"] for (row) in (rows)], dtype = float)
    p = np.asarray(probabilities, dtype = float)
    if (not np.isfinite(p).all() or np.any(p < 0) or np.any(p > 1)):
        raise ValueError("Invalid prediction probabilities")
    order = sorted(
        range(len(rows)),
        key = lambda i: (-float(p[i]), rows[i]["row_id"].encode("utf-8")),
    )
    count = math.ceil(0.25 * len(rows))
    event_count = int(y.sum())
    captured = int(y[order[:count]].sum())
    calibration = []
    for (index) in (range(10)):
        selected = (p >= index / 10) & (p < (index + 1) / 10 if (index < 9) else p <= 1)
        calibration.append(
            {
                "lower": index / 10,
                "upper": (index + 1) / 10,
                "n": int(selected.sum()),
                "mean_probability": float(p[selected].mean())
                if (selected.any())
                else None,
                "event_fraction": float(y[selected].mean())
                if (selected.any())
                else None,
            }
        )
    return {
        "brier": float(np.mean((p - y) ** 2)),
        "calibration": calibration,
        "budget_fraction": 0.25,
        "allocated_count": count,
        "capture": None if (event_count == 0) else captured / event_count,
        "captured_events": captured,
        "total_events": event_count,
        "allocated_row_ids": [rows[i]["row_id"] for (i) in (order[:count])],
        "full_order": [rows[i]["row_id"] for (i) in (order)],
        "tie_rule": "descending risk then frozen UTF-8 row identity",
    }

def fit_family(rows, name):
    allowed_development(rows)
    (retained, x, omitted) = available(rows, name)
    if (set((row["cohort"] for (row) in (retained))) != set(DEVELOPMENT)):
        return {
            "status": "UNAVAILABLE_BOTH_COMPLETE_DEVELOPMENT_CONTEXTS_REQUIRED",
            "model": name,
            "omitted": omitted,
        }
    y = np.asarray([row["target"] for (row) in (retained)], dtype = float)
    trials = []
    for (penalty) in (PENALTIES):
        folds = []
        for (heldout) in (DEVELOPMENT):
            train = np.asarray([row["cohort"] != heldout for (row) in (retained)])
            test = ~train
            model = ridge_fit(x[train], y[train], penalty, len(MODELS[name]))
            probabilities = predict(model, x[test])
            folds.append(
                {
                    "heldout": heldout,
                    "training": [
                        cohort for (cohort) in (DEVELOPMENT) if (cohort != heldout)
                    ],
                    "training_rows": int(train.sum()),
                    "heldout_rows": int(test.sum()),
                    "heldout_brier": float(np.mean((probabilities - y[test]) ** 2)),
                    "training_only_transform": {
                        "centre": model["centre"],
                        "scale": model["scale"],
                    },
                }
            )
        trials.append(
            {
                "penalty": penalty,
                "folds": folds,
                "mean_heldout_brier": sum((fold["heldout_brier"] for (fold) in (folds)))
                / len(folds),
            }
        )
    selected = select_penalty(trials)
    fitted = ridge_fit(x, y, selected, len(MODELS[name]))
    return {
        "status": "FITTED_DEVELOPMENT_ONLY",
        "model": name,
        "columns": columns(name),
        "selected_penalty": selected,
        "trials": trials,
        "fitted": fitted,
        "retained_row_ids": [row["row_id"] for (row) in (retained)],
        "omitted": omitted,
        "training_contexts": list(DEVELOPMENT),
        "excluded_graph_control_seed": HOLDOUT,
        "transfer_evaluated": False,
        "in_sample_metrics_not_transfer": diagnostics(retained, predict(fitted, x)),
    }
