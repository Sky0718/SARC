import json
import math
import numpy as np
from scipy.optimize import brentq, minimize
from scipy.special import expit, logit
from sklearn.preprocessing import StandardScaler

SEED = 20260930
TEMPERATURE = 0.5
PROBABILITY_CLIP = 1e-06
COEFFICIENT_BOUND = 10.0
LOCALISATION = (
    "added_mutant_contact_edges",
    "removed_mutant_contact_edges",
    "mutant_candidates_touched_fraction",
    "mean_mutant_relative_edit_burden",
    "max_mutant_relative_edit_burden",
    "mean_mutant_deg_onehop_edit_burden",
)
INPUT_FEATURES = (
    "log_edit_count",
    "mutated_edit_fraction",
    "deg_edit_fraction",
    "between_module_edit_fraction",
    "mean_endpoint_log_degree",
    "mean_consumed_cost_change",
)
INVARIANT = (
    "affine_slope_signed_log1p",
    "affine_intercept_over_old_mean_scale",
    "residual_mean_absolute_over_old_mass",
    "residual_rms_over_old_mass",
    "residual_max_absolute_over_old_mass",
    "residual_concentration",
    "positive_residual_mass_over_old_mass",
    "negative_residual_mass_over_old_mass",
    "sigmoid_oldscore_weighted_residual_over_old_mass",
    "squared_standardised_oldscore_weighted_residual_over_old_mass",
    "candidate_entry_fraction",
    "candidate_exit_fraction",
    "shared_candidate_fraction",
    "degenerate_old_variance",
    "no_shared_candidates",
    "affine_reversal",
)
TRACKS = {
    "old_inputs": [*INPUT_FEATURES, *LOCALISATION],
    "invariant_inputs": [*INPUT_FEATURES, *LOCALISATION, *INVARIANT],
}
METHODS = ("DawnRank", "PersonaDrive", "PRODIGY")

def require(condition, message):
    if (not condition):
        raise ValueError(message)

def patient(row):
    return (row["cohort"], row["sample_id"])

def row_ids(rows):
    identities = [row["row_id"] for (row) in (rows)]
    require(
        identities and len(identities) == len(set(identities)),
        "Nonempty unique row identities required",
    )
    require(
        all(
            (
                isinstance(value, str) and value
                for (row) in (rows)
                for (value) in ((row["row_id"], row["cohort"], row["sample_id"]))
            )
        ),
        "Explicit nonempty identities required",
    )
    return identities

def make_crossfit_plan(rows, partition_order, n_splits = 3, seed = SEED):
    row_ids(rows)
    require(
        n_splits == 3 and seed == SEED, "Frozen threefold patient crossfit required"
    )
    assignment = {}
    for (cohort) in (sorted({row["cohort"] for (row) in (rows)})):
        groups = {patient(row) for (row) in (rows) if (row["cohort"] == cohort)}
        require(
            len(groups) >= n_splits,
            "At least three training patients per cohort required",
        )
        priority = partition_order[cohort]
        require(
            len(priority) == len(set(priority))
            and {sample for (_, sample) in (groups)} <= set(priority),
            "Complete unique archived R7 sample priority required",
        )
        ordered = [
            (cohort, sample) for (sample) in (priority) if ((cohort, sample) in groups)
        ]
        assignment.update(
            {group: index % n_splits for (index, group) in (enumerate(ordered))}
        )
    return [
        {
            "fold_id": fold,
            "fit_row_ids": [
                row["row_id"] for (row) in (rows) if (assignment[patient(row)] != fold)
            ],
            "predict_row_ids": [
                row["row_id"] for (row) in (rows) if (assignment[patient(row)] == fold)
            ],
        }
        for (fold) in (range(n_splits))
    ]

def validate_identity(identity):
    require(
        isinstance(identity, dict)
        and set(identity) == {"baseline_spec", "feature_names", "source_id", "data_id"},
        "Source-bound anchor identity required",
    )
    require(
        identity["feature_names"]
        == [*LOCALISATION, "method_DawnRank", "method_PersonaDrive"],
        "Exact R2 localisation-only feature order required",
    )
    spec = identity["baseline_spec"]
    require(
        isinstance(spec, dict) and spec.get("family") == "extra_trees",
        "R2 ExtraTrees anchor required",
    )
    require(
        spec.get("params")
        in (
            {"depth": 6, "min_leaf": 32, "trees": 96},
            {"depth": 10, "min_leaf": 64, "trees": 96},
        ),
        "Exact R2 ExtraTrees parameters required",
    )
    require(
        spec.get("spec_id") == "extra_trees_" + str(spec["params"]["depth"]),
        "Anchor specification identity differs",
    )
    require(
        all(
            (
                isinstance(identity[key], str) and bool(identity[key])
                for (key) in (("source_id", "data_id"))
            )
        ),
        "Named source and carrier identities required",
    )

def validate_anchor(anchor, rows, trainrows, crossfit, partition_order):
    require(
        isinstance(anchor, dict), "Anchor provenance package required, not bare scores"
    )
    require(
        anchor.get("row_ids") == row_ids(rows), "Anchor prediction row order differs"
    )
    require(
        anchor.get("training_row_ids") == row_ids(trainrows),
        "Anchor training rows must equal current outer training rows",
    )
    validate_identity(anchor.get("identity"))
    probabilities = np.asarray(anchor.get("scores"), dtype = float)
    require(
        probabilities.shape == (len(rows),)
        and np.isfinite(probabilities).all()
        and ((probabilities >= 0) & (probabilities <= 1)).all(),
        "Finite anchor probabilities required",
    )
    if (crossfit):
        require(
            anchor.get("mode") == "patient_crossfit",
            "Training anchors must be patient crossfit within current training, never saved outer OOF",
        )
        require(
            anchor.get("crossfit_folds")
            == make_crossfit_plan(trainrows, partition_order),
            "Anchor fold populations differ from deterministic training-only patient crossfit",
        )
    else:
        require(
            anchor.get("mode") == "outer_training_fit",
            "Test anchor must come from exact current training fit",
        )
        require(
            not set(row_ids(rows)).intersection(row_ids(trainrows)),
            "Test rows overlap current training rows",
        )
    return probabilities

def build_crossfit_anchors(rows, baseline_identity, fit_predict_base, partition_order):
    validate_identity(baseline_identity)
    plan = make_crossfit_plan(rows, partition_order)
    index = {row["row_id"]: position for (position, row) in (enumerate(rows))}
    scores = np.full(len(rows), np.nan)
    receipts = []
    for (fold) in (plan):
        fit_rows = [rows[index[identity]] for (identity) in (fold["fit_row_ids"])]
        predict_rows = [rows[index[identity]] for (identity) in (fold["predict_row_ids"])]
        result = fit_predict_base(fit_rows, predict_rows, baseline_identity)
        predicted = validate_anchor(
            result, predict_rows, fit_rows, False, partition_order
        )
        require(
            result["identity"] == baseline_identity,
            "Crossfit callback changed anchor specification or sources",
        )
        scores[[index[identity] for (identity) in (fold["predict_row_ids"])]] = predicted
        receipts.append(
            {key: value for (key, value) in (result.items()) if (key != "scores")}
        )
    package = {
        "mode": "patient_crossfit",
        "row_ids": row_ids(rows),
        "training_row_ids": row_ids(rows),
        "scores": scores.tolist(),
        "identity": baseline_identity,
        "crossfit_folds": plan,
        "fit_receipts": receipts,
    }
    validate_anchor(package, rows, rows, True, partition_order)
    return package

def specifications(feature_names):
    names = list(feature_names)
    require(
        names in TRACKS.values(),
        "Exact frozen old-input or invariant-input residual feature order required",
    )
    result = [
        {
            "spec_id": "r7_anchor_baseline",
            "ridge": 0.0,
            "budget_weight": 0.0,
            "alpha": 0,
            "feature_names": names,
        }
    ]
    result.extend(
        (
            {
                "spec_id": "r7_residual_ridge_" + str(ridge) + "_budget_" + str(weight),
                "ridge": ridge,
                "budget_weight": weight,
                "alpha": 1,
                "feature_names": names,
            }
            for (ridge) in ((0.03, 0.3))
            for (weight) in ((0.0, 0.5))
        )
    )
    return result

def design_matrices(trainrows, testrows, trainX, testX, feature_names):
    training, testing = (
        np.asarray(trainX, dtype = float),
        np.asarray(testX, dtype = float),
    )
    count = len(feature_names)
    require(
        training.shape == (len(trainrows), count + 2)
        and testing.shape == (len(testrows), count + 2),
        "Exact declared feature and method-bit dimensions required",
    )
    require(
        np.isfinite(training).all() and np.isfinite(testing).all(),
        "Finite full-support input features required",
    )
    for (rows, matrix) in (((trainrows, training), (testrows, testing))):
        require(all((row["method"] in METHODS for (row) in (rows))), "Unknown method")
        expected = np.asarray(
            [
                [
                    float(row["method"] == "DawnRank"),
                    float(row["method"] == "PersonaDrive"),
                ]
                for (row) in (rows)
            ]
        )
        require(
            np.array_equal(matrix[:, count:], expected),
            "Method bits differ from row methods",
        )
    scaler = StandardScaler().fit(training[:, :count])
    transformed = []
    for (matrix) in ((training, testing)):
        continuous, bits = (scaler.transform(matrix[:, :count]), matrix[:, count:])
        transformed.append(
            np.column_stack(
                (
                    np.ones(len(matrix)),
                    continuous,
                    bits,
                    continuous * bits[:, [0]],
                    continuous * bits[:, [1]],
                )
            )
        )
    names = [
        "intercept",
        *feature_names,
        "method_DawnRank",
        "method_PersonaDrive",
        *[name + ":DawnRank" for (name) in (feature_names)],
        *[name + ":PersonaDrive" for (name) in (feature_names)],
    ]
    return (
        *transformed,
        {
            "mean": scaler.mean_.tolist(),
            "scale": scaler.scale_.tolist(),
            "variance": scaler.var_.tolist(),
            "n_samples_seen": int(scaler.n_samples_seen_),
            "expanded_feature_names": names,
        },
    )

def soft_budget_capture(scores, labels, temperature = TEMPERATURE):
    scores, labels = (np.asarray(scores, dtype = float), np.asarray(labels, dtype = float))
    require(
        scores.ndim == 1
        and scores.size
        and (labels.shape == scores.shape)
        and np.isfinite(scores).all(),
        "Finite nonempty budget group required",
    )
    require(
        np.isfinite(labels).all()
        and np.isin(labels, [0.0, 1.0]).all()
        and (temperature == TEMPERATURE),
        "Binary targets and frozen temperature required",
    )
    selected = math.ceil(0.25 * len(scores))
    events = float(labels.sum())
    if (selected == len(scores)):
        return (
            1.0 if (events) else 0.0,
            np.zeros_like(scores),
            {
                "n_rows": len(scores),
                "selected_budget": selected,
                "soft_selected": float(selected),
                "events": events,
                "threshold": None,
            },
        )
    radius = temperature * (50.0 + math.log(len(scores)))
    threshold = brentq(
        lambda value: float(expit((scores - value) / temperature).sum()) - selected,
        float(scores.min() - radius),
        float(scores.max() + radius),
        xtol = 1e-12,
        rtol = 1e-14,
    )
    selected_probabilities = expit((scores - threshold) / temperature)
    sensitivity = selected_probabilities * (1.0 - selected_probabilities) / temperature
    require(
        abs(float(selected_probabilities.sum()) - selected) < 1e-08
        and float(sensitivity.sum()) > 0,
        "Budget threshold failed numerical admission",
    )
    if (events):
        capture = float(labels @ selected_probabilities / events)
        gradient = (
            sensitivity
            * (labels - float(labels @ sensitivity) / float(sensitivity.sum()))
            / events
        )
    else:
        capture, gradient = (0.0, np.zeros_like(scores))
    return (
        capture,
        gradient,
        {
            "n_rows": len(scores),
            "selected_budget": selected,
            "soft_selected": float(selected_probabilities.sum()),
            "events": events,
            "threshold": float(threshold),
        },
    )

def objective(weights, design, anchor_logits, labels, groups, ridge, budget_weight):
    scores = anchor_logits + design @ weights
    probabilities = expit(scores)
    loss = float(
        np.mean(np.logaddexp(0.0, scores) - labels * scores)
        + 0.5 * ridge * (weights @ weights)
    )
    score_gradient = (probabilities - labels) / len(labels)
    details = []
    for (group) in (sorted(set(groups))):
        indices = np.flatnonzero(np.asarray(groups) == group)
        capture, capture_gradient, detail = soft_budget_capture(
            scores[indices], labels[indices]
        )
        loss -= budget_weight * capture / len(set(groups))
        score_gradient[indices] -= budget_weight * capture_gradient / len(set(groups))
        details.append({"cohort": group, "soft_capture": capture, **detail})
    gradient = design.T @ score_gradient + ridge * weights
    require(
        math.isfinite(loss) and np.isfinite(gradient).all(),
        "Nonfinite residual objective or derivative",
    )
    return (loss, gradient, details)

def fit_predict(
    trainrows,
    testrows,
    trainX,
    testX,
    y,
    trainanchor,
    testanchor,
    spec,
    partition_order,
):
    require(
        isinstance(spec, dict)
        and spec in specifications(spec.get("feature_names", [])),
        "Exact finite residual specification required",
    )
    row_ids(trainrows)
    row_ids(testrows)
    test_base = validate_anchor(testanchor, testrows, trainrows, False, partition_order)
    labels = np.asarray(y, dtype = float)
    require(
        labels.shape == (len(trainrows),)
        and np.isfinite(labels).all()
        and np.isin(labels, [0.0, 1.0]).all(),
        "Complete binary training labels required",
    )
    diagnostics = {
        "specification": spec,
        "anchor_identity": testanchor["identity"],
        "seed": SEED,
        "actual_anchor_fit_performed_here": False,
        "training_rows": len(trainrows),
        "training_patients": len({patient(row) for (row) in (trainrows)}),
        "test_rows": len(testrows),
        "temperature": TEMPERATURE,
        "probability_clip": PROBABILITY_CLIP,
    }
    if (spec["alpha"] == 0):
        return {
            "score": test_base.tolist(),
            "baseline_score": test_base.tolist(),
            "diagnostics": {
                **diagnostics,
                "status": "EXACT_BASELINE_IDENTITY",
                "optimisation_performed": False,
            },
        }
    training_base = validate_anchor(
        trainanchor, trainrows, trainrows, True, partition_order
    )
    require(
        trainanchor["identity"] == testanchor["identity"],
        "Training and test anchor provenance identities differ",
    )
    training, testing, scale = design_matrices(
        trainrows, testrows, trainX, testX, spec["feature_names"]
    )
    anchor_logits = logit(
        np.clip(training_base, PROBABILITY_CLIP, 1.0 - PROBABILITY_CLIP)
    )
    groups = [row["cohort"] for (row) in (trainrows)]
    evaluate = lambda weights: objective(
        weights,
        training,
        anchor_logits,
        labels,
        groups,
        spec["ridge"],
        spec["budget_weight"],
    )
    result = minimize(
        lambda weights: evaluate(weights)[:2],
        np.zeros(training.shape[1]),
        method = "L-BFGS-B",
        jac = True,
        bounds = [(-COEFFICIENT_BOUND, COEFFICIENT_BOUND)] * training.shape[1],
        options = {"maxiter": 300, "ftol": 1e-12, "gtol": 1e-07, "maxls": 50},
    )
    final_loss, final_gradient, budget_details = evaluate(result.x)
    optimiser = {
        "success": bool(result.success),
        "status": int(result.status),
        "message": str(result.message),
        "iterations": int(result.nit),
        "evaluations": int(result.nfev),
        "objective": final_loss,
        "maximum_absolute_gradient": float(np.abs(final_gradient).max()),
        "coefficient_bound": COEFFICIENT_BOUND,
        "coefficients_at_bound": int(
            np.count_nonzero(np.abs(result.x) >= COEFFICIENT_BOUND - 1e-08)
        ),
    }
    require(
        result.success and np.isfinite(result.x).all(),
        "Residual optimisation failed: " + json.dumps(optimiser, sort_keys = True),
    )
    predictions = expit(
        logit(np.clip(test_base, PROBABILITY_CLIP, 1.0 - PROBABILITY_CLIP))
        + testing @ result.x
    )
    require(np.isfinite(predictions).all(), "Nonfinite residual predictions")
    diagnostics.update(
        status = "FITTED_RESIDUAL",
        optimisation_performed = True,
        optimiser = optimiser,
        coefficients = result.x.tolist(),
        scaler = scale,
        crossfit_folds = trainanchor["crossfit_folds"],
        budget_groups = budget_details,
        ridge_includes_intercept = True,
        unseen_cohort_robustness_guaranteed = False,
    )
    return {
        "score": predictions.tolist(),
        "baseline_score": test_base.tolist(),
        "diagnostics": diagnostics,
    }
