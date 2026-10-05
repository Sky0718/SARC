import json
import math
import time
from dataclasses import dataclass
from fractions import Fraction
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

SEED = 20260930
FEATURES = (
    "log_edit_count",
    "mutated_edit_fraction",
    "deg_edit_fraction",
    "between_module_edit_fraction",
    "mean_endpoint_log_degree",
    "mean_consumed_cost_change",
)
LOCALISATION = (
    "added_mutant_contact_edges",
    "removed_mutant_contact_edges",
    "mutant_candidates_touched_fraction",
    "mean_mutant_relative_edit_burden",
    "max_mutant_relative_edit_burden",
    "mean_mutant_deg_onehop_edit_burden",
)
MECHANISM = (
    "mean_normalised_absolute_direct_delta",
    "max_normalised_absolute_direct_delta",
    "mean_normalised_absolute_reweighting_delta",
    "max_normalised_absolute_reweighting_delta",
    "mean_normalised_absolute_interaction_delta",
    "max_normalised_absolute_interaction_delta",
    "mean_normalised_absolute_total_delta",
    "max_normalised_absolute_total_delta",
    "similarity_change_fraction",
    "mean_absolute_similarity_change",
    "candidate_entry_fraction",
    "candidate_exit_fraction",
)
TRACKS = {
    "edit_volume_only": [FEATURES[0]],
    "degree_incidence_only": [FEATURES[1], FEATURES[4]],
    "additive_volume_degree": [FEATURES[0], FEATURES[1], FEATURES[4]],
    "strict_input_only": list(FEATURES),
    "localisation_only": list(LOCALISATION),
    "strict_plus_localisation": [*FEATURES, *LOCALISATION],
    "strict_localisation_mechanism": [*FEATURES, *LOCALISATION, *MECHANISM],
}
METHODS = ("DawnRank", "PersonaDrive", "PRODIGY")
MASTERS = (104729, 130363, 155921, 196613, 228017)
COHORTS = ("COAD_CCLE", "LUAD_CCLE")
TRANSITIONS = ("11_0_to_11_5", "11_5_to_12_0")
EXPECTED_RUNTIME = {
    "numpy": "2.2.6",
    "scipy": "1.15.3",
    "sklearn": "1.7.2",
    "threadpoolctl": "3.6.0",
    "xgboost": "3.2.0",
    "catboost": "1.2.10",
}
EXPECTED_PYTHON = (3, 10, 6)
SELECTION_RULE = "MAX_EQUAL_INNER_CONTEXT_CAPTURE25_THEN_MIN_BRIER_THEN_SPEC_ID"
SPLIT_PREFIX = "R6_CALIBRATION_20260930"
FAMILIES = ("xgb_ndcg_topk", "xgb_pairwise_mean", "xgb_binary_logistic")
CALIBRATION_PENALTY = 0.005
COMMON_PARAMS = dict(
    n_estimators = 120,
    learning_rate = 0.05,
    min_child_weight = 16,
    reg_lambda = 10,
    reg_alpha = 0,
    subsample = 1,
    colsample_bytree = 1,
    max_bin = 128,
    tree_method = "hist",
    n_jobs = 1,
    random_state = SEED,
)

def require(condition, message):
    if (not condition):
        raise ValueError(message)

def number(value):
    if (isinstance(value, dict) and {"numerator", "denominator"} <= set(value)):
        value = Fraction(value["numerator"], value["denominator"])
    value = float(value)
    require(math.isfinite(value), "Nonfinite feature")
    return value

def vector(row, track):
    require(track in TRACKS and row["method"] in METHODS, "Undeclared track or method")
    result = []
    for (name) in (TRACKS[track]):
        value = row["features"].get(name)
        if (name == "mean_consumed_cost_change" and row["method"] != "PRODIGY"):
            require(
                value is None
                and row["features"].get("cost_structurally_not_applicable") is True,
                "Binary method structural cost contract differs",
            )
            value = 0.0
        if (value is None):
            return None
        result.append(number(value))
    return result + [
        float(row["method"] == "DawnRank"),
        float(row["method"] == "PersonaDrive"),
    ]

def specifications():
    return [
        dict(
            spec_id = f"r6_{family}_depth_{depth}",
            family = family,
            params = {"max_depth": depth},
            complexity_order = index,
        )
        for (index, (family, depth)) in (enumerate(
            ((family, depth) for (family) in (FAMILIES) for (depth) in ((2, 4)))
        ))
    ]

specs = specifications

def biological_halves(rows, partition_order):
    require(rows, "Nonempty usable training rows required")
    ids = [row["row_id"] for (row) in (rows)]
    require(
        all((isinstance(identity, str) and identity for (identity) in (ids)))
        and len(ids) == len(set(ids)),
        "Unique nonempty training row IDs required",
    )
    cohort_samples = {}
    for (row) in (rows):
        cohort, sample = (row["cohort"], row["sample_id"])
        require(
            cohort in COHORTS and isinstance(sample, str) and sample,
            "Development cohort and biological sample IDs required",
        )
        cohort_samples.setdefault(cohort, set()).add(sample)
    maps = [{}, {}]
    for (cohort, samples) in (sorted(cohort_samples.items())):
        priority = partition_order[cohort]
        require(
            len(priority) == len(set(priority)) and samples <= set(priority),
            "Complete unique archived R6 sample priority required",
        )
        ordered = [sample for (sample) in (priority) if (sample in samples)]
        require(
            len(ordered) >= 2,
            "Every required cohort needs biological samples in both halves",
        )
        for (half) in (range(2)):
            maps[half][cohort] = ordered[half::2]
    halves = []
    for (half) in (range(2)):
        indexes = [
            index
            for (index, row) in (enumerate(rows))
            if (row["sample_id"] in maps[half][row["cohort"]])
        ]
        halves.append(
            {
                "half": half,
                "cohort_samples": maps[half],
                "row_ids": [ids[index] for (index) in (indexes)],
                "indexes": indexes,
            }
        )
    require(
        set(halves[0]["row_ids"]).isdisjoint(halves[1]["row_ids"])
        and set(halves[0]["row_ids"]) | set(halves[1]["row_ids"]) == set(ids),
        "Calibration halves must be disjoint and exhaustive",
    )
    return halves

def calibration_objective(parameters, z, labels):
    slope, intercept = parameters
    logits = slope * z + intercept
    loss = (
        np.mean(np.logaddexp(0.0, logits) - labels * logits)
        + CALIBRATION_PENALTY * slope * slope
    )
    residual = expit(logits) - labels
    gradient = np.array(
        [np.mean(residual * z) + 2 * CALIBRATION_PENALTY * slope, np.mean(residual)]
    )
    return (float(loss), gradient)

@dataclass
class SigmoidCalibrator:
    mean: float
    sd: float
    slope: float | None
    intercept: float | None
    constant: float | None

    def predict(self, scores):
        scores = np.asarray(scores, dtype = np.float64)
        require(
            scores.ndim == 1 and np.isfinite(scores).all(),
            "Finite one-dimensional prediction margins required",
        )
        if (self.constant is not None):
            return np.full(scores.shape, self.constant, dtype = np.float64)
        values = expit(self.slope * ((scores - self.mean) / self.sd) + self.intercept)
        require(np.isfinite(values).all(), "Nonfinite calibrated probabilities")
        return values

def fit_calibrator(scores, labels, row_ids):
    scores, labels = (
        np.asarray(scores, dtype = np.float64),
        np.asarray(labels, dtype = np.float64),
    )
    require(
        scores.ndim == labels.ndim == 1
        and len(scores) == len(labels) == len(row_ids)
        and (len(scores) > 0),
        "Aligned nonempty calibration rows required",
    )
    require(
        np.isfinite(scores).all() and np.isin(labels, (0, 1)).all(),
        "Finite margins and binary calibration labels required",
    )
    mean, sd, prevalence = (
        float(scores.mean()),
        float(scores.std(ddof = 0)),
        float(labels.mean()),
    )
    require(math.isfinite(mean) and math.isfinite(sd), "Nonfinite calibration moments")
    receipt = dict(
        row_ids = list(row_ids),
        raw_scores = scores.tolist(),
        labels = labels.astype(int).tolist(),
        n = len(scores),
        positive = int(labels.sum()),
        prevalence = prevalence,
        mean = mean,
        sd = sd,
    )
    if (prevalence in (0.0, 1.0) or sd == 0.0):
        mode = (
            "constant_class_prevalence"
            if (prevalence in (0.0, 1.0))
            else "zero_score_sd_prevalence"
        )
        receipt.update(
            mode = mode,
            slope = None,
            intercept = None,
            constant_probability = prevalence,
            convergence = {"success": True, "optimizer_used": False},
        )
        return (SigmoidCalibrator(mean, sd, None, None, prevalence), receipt)
    z = (scores - mean) / sd
    clipped = np.clip(prevalence, 1e-06, 1 - 1e-06)
    initial = np.array([1.0, math.log(clipped / (1 - clipped))])
    result = minimize(
        calibration_objective,
        initial,
        args = (z, labels),
        method = "L-BFGS-B",
        jac = True,
        bounds = ((0.0, None), (None, None)),
        options = {"maxiter": 500, "ftol": 1e-12, "gtol": 1e-08},
    )
    require(
        result.success
        and np.isfinite(result.x).all()
        and math.isfinite(float(result.fun))
        and np.isfinite(result.jac).all(),
        "Calibration optimizer failed or returned nonfinite values: "
        + str(result.message),
    )
    slope, intercept = map(float, result.x)
    require(slope >= 0, "Negative calibration slope forbidden")
    receipt.update(
        mode = "monotone_sigmoid",
        slope = slope,
        intercept = intercept,
        constant_probability = None,
        convergence = {
            "success": True,
            "optimizer_used": True,
            "status": int(result.status),
            "message": str(result.message),
            "iterations": int(result.nit),
            "function_evaluations": int(result.nfev),
            "objective": float(result.fun),
            "gradient": result.jac.tolist(),
            "initial": initial.tolist(),
            "method": "L-BFGS-B",
            "bounds": [[0.0, None], [None, None]],
            "options": {"maxiter": 500, "ftol": 1e-12, "gtol": 1e-08},
            "slope_squared_penalty": CALIBRATION_PENALTY,
        },
    )
    return (SigmoidCalibrator(mean, sd, slope, intercept, None), receipt)

def learner_parameters(spec, query_sizes):
    require(spec in specifications(), "Exact finite R6 specification required")
    require(
        query_sizes
        and all((isinstance(size, int) and size > 0 for (size) in (query_sizes.values()))),
        "Positive cohort query lengths required",
    )
    result = {**COMMON_PARAMS, **spec["params"]}
    if (spec["family"] != "xgb_binary_logistic"):
        require(
            len(set(query_sizes.values())) == 1,
            "Ranking requires equal cohort query lengths; no budget adaptation",
        )
    if (spec["family"] == "xgb_ndcg_topk"):
        result.update(
            objective = "rank:ndcg",
            lambdarank_pair_method = "topk",
            lambdarank_num_pair_per_sample = math.ceil(
                0.25 * next(iter(query_sizes.values()))
            ),
            ndcg_exp_gain = False,
        )
    elif (spec["family"] == "xgb_pairwise_mean"):
        result.update(
            objective = "rank:pairwise",
            lambdarank_pair_method = "mean",
            lambdarank_num_pair_per_sample = 8,
        )
    else:
        result.update(objective = "binary:logistic")
    return result

@dataclass
class Component:
    learner: object
    constant_base_score: float | None
    calibrator: SigmoidCalibrator

    def margin(self, matrix):
        if (self.constant_base_score is not None):
            return np.full(len(matrix), self.constant_base_score, dtype = np.float64)
        return np.asarray(
            self.learner.predict(matrix, output_margin = True), dtype = np.float64
        )

    def predict(self, matrix):
        return self.calibrator.predict(self.margin(matrix))

class CalibratedEnsemble:
    classes_ = np.array([0, 1])

    def __init__(self, components, feature_count):
        self.components, self.n_features_in_ = (components, feature_count)

    def predict_components(self, matrix):
        matrix = np.asarray(matrix, dtype = np.float64)
        require(
            matrix.ndim == 2
            and matrix.shape[1] == self.n_features_in_
            and np.isfinite(matrix).all(),
            "Finite declared feature matrix required",
        )
        result = []
        for (component) in (self.components):
            scores = component.margin(matrix) if (len(matrix)) else np.empty(0)
            result.append(
                {
                    "raw_scores": scores.tolist(),
                    "calibrated_probabilities": component.calibrator.predict(
                        scores
                    ).tolist(),
                }
            )
        return result

    def predict_proba(self, matrix):
        matrix = np.asarray(matrix, dtype = np.float64)
        require(
            matrix.ndim == 2
            and matrix.shape[1] == self.n_features_in_
            and np.isfinite(matrix).all(),
            "Finite declared feature matrix required",
        )
        if (not len(matrix)):
            return np.empty((0, 2), dtype = np.float64)
        probabilities = np.mean(
            [component.predict(matrix) for (component) in (self.components)], axis = 0
        )
        require(
            np.isfinite(probabilities).all()
            and ((probabilities >= 0) & (probabilities <= 1)).all(),
            "Invalid ensemble probability",
        )
        return np.column_stack((1 - probabilities, probabilities))

def fit_calibrated(spec, track, usable_rows, vectors, partition_order):
    from xgboost import XGBClassifier, XGBRanker

    started, cpu_started = (time.monotonic(), time.process_time())
    require(spec in specifications() and track in TRACKS, "Undeclared R6 fit")
    rows = list(usable_rows)
    halves = biological_halves(rows, partition_order)
    if (isinstance(vectors, dict)):
        vectors = [vectors[row["row_id"]] for (row) in (rows)]
    matrix = np.asarray(vectors, dtype = np.float64)
    require(
        matrix.shape == (len(rows), len(TRACKS[track]) + 2)
        and np.isfinite(matrix).all(),
        "Aligned finite R6 feature matrix required",
    )
    expected = [vector(row, track) for (row) in (rows)]
    require(
        all((value is not None for (value) in (expected)))
        and np.array_equal(matrix, np.asarray(expected, dtype = np.float64)),
        "Features must exactly match input-only declared row vectors",
    )
    require(
        all((row["target"] in (0, 1) for (row) in (rows))),
        "Every usable row needs its unchanged binary target",
    )
    labels = np.asarray([row["target"] for (row) in (rows)], dtype = np.int32)
    components, receipts = ([], [])
    for (half) in (range(2)):
        component_started, component_cpu = (time.monotonic(), time.process_time())
        fit_indexes = sorted(
            halves[half]["indexes"],
            key = lambda index: (rows[index]["cohort"], rows[index]["row_id"]),
        )
        calibration_indexes = halves[1 - half]["indexes"]
        query_sizes = {
            cohort: sum((rows[index]["cohort"] == cohort for (index) in (fit_indexes)))
            for (cohort) in (sorted(halves[half]["cohort_samples"]))
        }
        params = learner_parameters(spec, query_sizes)
        fit_labels = labels[fit_indexes]
        constant = float(fit_labels.mean()) if (len(np.unique(fit_labels)) == 1) else None
        learner, config = (None, None)
        learner_started, learner_cpu = (time.monotonic(), time.process_time())
        if (constant is None):
            if (spec["family"] == "xgb_binary_logistic"):
                learner = XGBClassifier(**params)
                learner.fit(matrix[fit_indexes], fit_labels)
            else:
                learner = XGBRanker(**params)
                learner.fit(
                    matrix[fit_indexes], fit_labels, group = list(query_sizes.values())
                )
            config = json.loads(learner.get_booster().save_config())
        fit_cost = dict(
            elapsed_seconds = time.monotonic() - learner_started,
            process_cpu_seconds = time.process_time() - learner_cpu,
        )
        temporary = Component(learner, constant, None)
        scores = temporary.margin(matrix[calibration_indexes])
        calibration_started, calibration_cpu = (time.monotonic(), time.process_time())
        calibrator, calibration = fit_calibrator(
            scores,
            labels[calibration_indexes],
            [rows[index]["row_id"] for (index) in (calibration_indexes)],
        )
        calibration_cost = dict(
            elapsed_seconds = time.monotonic() - calibration_started,
            process_cpu_seconds = time.process_time() - calibration_cpu,
        )
        temporary.calibrator = calibrator
        components.append(temporary)
        receipts.append(
            dict(
                fit_half = half,
                calibration_half = 1 - half,
                fit_row_ids = [rows[index]["row_id"] for (index) in (fit_indexes)],
                calibration_row_ids = calibration["row_ids"],
                cohort_query_sizes = query_sizes,
                learner_params = params,
                booster_config = config,
                constant_base_score = constant,
                fit_positive = int(fit_labels.sum()),
                fit_n = len(fit_labels),
                calibration = calibration,
                learner_cost = fit_cost,
                calibration_cost = calibration_cost,
                elapsed_seconds = time.monotonic() - component_started,
                process_cpu_seconds = time.process_time() - component_cpu,
            )
        )
    model = CalibratedEnsemble(components, matrix.shape[1])
    diagnostics = dict(
        schema = "R6_GROUP_CALIBRATION_01",
        specification = spec,
        track = track,
        train_row_ids = [row["row_id"] for (row) in (rows)],
        biological_halves = [
            {key: value for (key, value) in (half.items()) if (key != "indexes")}
            for (half) in (halves)
        ],
        components = receipts,
        no_all_data_refit = True,
        ensemble_rule = "MEAN_OF_TWO_COMPONENT_CALIBRATED_PROBABILITIES",
        elapsed_seconds = time.monotonic() - started,
        process_cpu_seconds = time.process_time() - cpu_started,
    )
    return (model, diagnostics)
