import json
import math

SEED = 20260930
PAIR_CAP = 8192
METHODS = ("DawnRank", "PersonaDrive", "PRODIGY")

def require(value, message):
    if (not value):
        raise ValueError(message)

def specifications():
    values = []
    for (alpha) in ((0.5, 2.0)):
        for (penalty) in ((0.003, 0.03)):
            values.append(
                (
                    "pairwise_a" + str(alpha) + "_l" + str(penalty),
                    "composite_pairwise_logistic",
                    {
                        "alpha": alpha,
                        "penalty": penalty,
                        "pair_cap": PAIR_CAP,
                        "max_iter": 500,
                    },
                )
            )
    for (penalty) in ((0.003, 0.03)):
        values.append(
            (
                "expert_logistic_l" + str(penalty),
                "method_expert_logistic",
                {"penalty": penalty, "max_iter": 1000},
            )
        )
    for (leaves, minimum) in (((7, 64), (15, 128))):
        values.append(
            (
                "expert_hist_" + str(leaves),
                "method_expert_hist",
                {"leaves": leaves, "iterations": 150, "min_leaf": minimum},
            )
        )
    return [
        {
            "spec_id": identity,
            "family": family,
            "params": parameters,
            "complexity_order": order,
        }
        for (order, (identity, family, parameters)) in (enumerate(values))
    ]

def specs():
    return specifications()

def arrays(x, continuous_count, y = None):
    import numpy as np

    values = np.asarray(x, dtype = float)
    require(
        type(continuous_count) is int and continuous_count > 0,
        "Positive declared continuous feature count required",
    )
    require(
        values.ndim == 2 and values.shape[1] == continuous_count + 2,
        "Only declared continuous features and two method indicators allowed",
    )
    require(np.isfinite(values).all(), "Nonfinite predictor is not available input")
    bits = values[:, -2:]
    require(
        np.isin(bits, (0.0, 1.0)).all() and (bits.sum(axis = 1) <= 1).all(),
        "Invalid method indicator encoding",
    )
    if (y is None):
        return values
    targets = np.asarray(y)
    require(
        targets.ndim == 1
        and len(targets) == len(values)
        and (len(targets) > 0)
        and np.isin(targets, (0, 1)).all(),
        "Complete nonempty binary training targets required",
    )
    return (values, targets.astype(float))

def method_indices(x):
    import numpy as np

    return np.where(x[:, -2] == 1, 0, np.where(x[:, -1] == 1, 1, 2))

def training_pairs(y, cap = PAIR_CAP, seed = SEED):
    import numpy as np

    require(type(cap) is int and cap > 0, "Positive fixed training-pair cap required")
    targets = np.asarray(y)
    require(
        targets.ndim == 1 and np.isin(targets, (0, 1)).all(),
        "Binary training pairs required",
    )
    positive = np.flatnonzero(targets == 1)
    negative = np.flatnonzero(targets == 0)
    total = len(positive) * len(negative)
    if (total == 0):
        return np.empty((0, 2), dtype = np.int64)
    slots = (
        np.arange(total, dtype = np.int64)
        if (total <= cap)
        else np.sort(np.random.default_rng(seed).choice(total, size = cap, replace = False))
    )
    return np.column_stack(
        (positive[slots // len(negative)], negative[slots % len(negative)])
    ).astype(np.int64)

def composite_loss_gradient(theta, design, y, pair_differences, alpha, penalty):
    import numpy as np
    from scipy.special import expit

    coefficients, intercept = (theta[:-1], theta[-1])
    scores = design @ coefficients + intercept
    residual = expit(scores) - y
    loss = np.mean(np.logaddexp(0.0, scores) - y * scores) + 0.5 * penalty * (
        coefficients @ coefficients
    )
    gradient = np.concatenate(
        (design.T @ residual / len(y) + penalty * coefficients, [np.mean(residual)])
    )
    if (len(pair_differences)):
        margins = pair_differences @ coefficients
        loss += alpha * np.mean(np.logaddexp(0.0, -margins))
        gradient[:-1] -= alpha * pair_differences.T @ expit(-margins) / len(margins)
    return (float(loss), gradient)

class CompositePairwiseLogistic:
    def __init__(
        self,
        continuous_count,
        training_rows,
        alpha,
        penalty,
        pair_cap = PAIR_CAP,
        max_iter = 500,
    ):
        self.continuous_count = continuous_count
        self.training_rows = training_rows
        self.alpha = alpha
        self.penalty = penalty
        self.pair_cap = pair_cap
        self.max_iter = max_iter

    def fit(self, x, y):
        import numpy as np
        from scipy.optimize import minimize
        from sklearn.preprocessing import StandardScaler

        values, targets = arrays(x, self.continuous_count, y)
        require(
            len(values) == self.training_rows,
            "Declared usable training population differs",
        )
        require(
            self.alpha > 0
            and self.penalty > 0
            and (self.max_iter == 500)
            and (self.pair_cap == PAIR_CAP),
            "Unfrozen composite objective",
        )
        self.classes_ = np.asarray([0, 1])
        self.n_features_in_ = values.shape[1]
        self.scaler_ = StandardScaler().fit(values[:, : self.continuous_count])
        design = np.column_stack(
            (self.scaler_.transform(values[:, : self.continuous_count]), values[:, -2:])
        )
        self.training_prevalence_ = float(targets.mean())
        self.pair_indices_ = training_pairs(targets, self.pair_cap)
        self.diagnostics_ = {
            "objective": "MEAN_BCE_PLUS_ALPHA_MEAN_PAIRWISE_LOGISTIC_PLUS_LAMBDA_HALF_L2_NONINTERCEPT",
            "training_rows": len(targets),
            "alpha": self.alpha,
            "penalty": self.penalty,
            "pair_cap": self.pair_cap,
            "pair_count": len(self.pair_indices_),
            "pair_seed": SEED,
            "pair_sampling": "ALL_PAIRS_IF_WITHIN_CAP_ELSE_SORTED_UNIFORM_WITHOUT_REPLACEMENT",
            "pair_indices": self.pair_indices_.tolist(),
            "intercept_penalised": False,
            "heldout_rows_used": False,
            "calibration": "JOINT_TRAINING_BERNOULLI_TERM_NO_HELDOUT_CALIBRATOR",
            "exact_top25_objective": False,
        }
        self.constant_probability_ = (
            self.training_prevalence_ if (len(np.unique(targets)) == 1) else None
        )
        if (self.constant_probability_ is not None):
            self.coef_ = np.zeros((1, values.shape[1]))
            self.intercept_ = np.asarray([0.0])
            self.diagnostics_.update(
                status = "CONSTANT_TRAINING_CLASS",
                optimiser = None,
                converged = True,
                iterations = 0,
            )
            return self
        differences = (
            design[self.pair_indices_[:, 0]] - design[self.pair_indices_[:, 1]]
        )
        initial = np.zeros(values.shape[1] + 1)
        initial[-1] = math.log(
            self.training_prevalence_ / (1 - self.training_prevalence_)
        )
        objective = lambda theta: composite_loss_gradient(
            theta, design, targets, differences, self.alpha, self.penalty
        )
        initial_loss = objective(initial)[0]
        result = minimize(
            objective,
            initial,
            method = "L-BFGS-B",
            jac = True,
            options = {
                "maxiter": self.max_iter,
                "maxls": 50,
                "ftol": 1e-12,
                "gtol": 1e-07,
            },
        )
        final_loss, final_gradient = objective(result.x)
        gradient_max = float(np.max(np.abs(final_gradient)))
        converged = bool(
            result.success
            and np.isfinite(result.x).all()
            and math.isfinite(final_loss)
            and (gradient_max <= 1e-05)
            and (final_loss <= initial_loss + 1e-10)
        )
        self.diagnostics_.update(
            status = "COMPLETE" if (converged) else "FAILED_CONVERGENCE",
            optimiser = "SCIPY_L_BFGS_B",
            solver_status = int(result.status),
            solver_message = str(result.message),
            converged = converged,
            iterations = int(result.nit),
            function_evaluations = int(result.nfev),
            initial_objective = initial_loss,
            final_objective = final_loss,
            gradient_max_absolute = gradient_max,
            max_iterations = self.max_iter,
            ftol = 1e-12,
            gtol = 1e-07,
            acceptance_gradient = 1e-05,
        )
        if (not converged):
            raise RuntimeError(
                "Composite optimiser failed declared convergence: "
                + json.dumps(self.diagnostics_, allow_nan = False)
            )
        self.coef_ = result.x[:-1].reshape(1, -1)
        self.intercept_ = np.asarray([result.x[-1]])
        return self

    def predict_proba(self, x):
        import numpy as np
        from scipy.special import expit

        require(
            hasattr(self, "diagnostics_") and self.diagnostics_["converged"],
            "Complete fitted estimator required",
        )
        values = arrays(x, self.continuous_count)
        if (self.constant_probability_ is not None):
            positive = np.full(len(values), self.constant_probability_)
        elif (not len(values)):
            positive = np.empty(0)
        else:
            design = np.column_stack(
                (
                    self.scaler_.transform(values[:, : self.continuous_count]),
                    values[:, -2:],
                )
            )
            positive = expit(design @ self.coef_[0] + self.intercept_[0])
        require(np.isfinite(positive).all(), "Nonfinite fitted probability")
        return np.column_stack((1 - positive, positive))

class MethodExperts:
    def __init__(self, family, parameters, continuous_count, training_rows):
        self.family = family
        self.parameters = dict(parameters)
        self.continuous_count = continuous_count
        self.training_rows = training_rows

    def fit(self, x, y):
        import numpy as np
        from sklearn.ensemble import HistGradientBoostingClassifier
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler

        values, targets = arrays(x, self.continuous_count, y)
        require(
            len(values) == self.training_rows,
            "Declared usable training population differs",
        )
        require(
            self.family in ("method_expert_logistic", "method_expert_hist"),
            "Unfrozen expert family",
        )
        self.classes_ = np.asarray([0, 1])
        self.n_features_in_ = values.shape[1]
        self.experts_ = {}
        self.diagnostics_ = {
            "family": self.family,
            "training_rows": len(values),
            "method_routing": "ONLY_DECLARED_FINAL_TWO_METHOD_INDICATORS",
            "pooled_model_fitted": False,
            "heldout_rows_used": False,
            "experts": {},
        }
        methods = method_indices(values)
        for (index, name) in (enumerate(METHODS)):
            selected = np.flatnonzero(methods == index)
            if (not len(selected)):
                self.diagnostics_["experts"][name] = {
                    "training_rows": 0,
                    "status": "NO_TRAINING_ROWS_PREDICTION_FORBIDDEN",
                }
                continue
            local_x, local_y = (
                values[selected, : self.continuous_count],
                targets[selected],
            )
            scaler = StandardScaler().fit(local_x)
            transformed = scaler.transform(local_x)
            constant = float(local_y.mean()) if (len(np.unique(local_y)) == 1) else None
            model, effective_c = (None, None)
            if (constant is None):
                if (self.family == "method_expert_logistic"):
                    effective_c = 1 / (len(selected) * self.parameters["penalty"])
                    model = LogisticRegression(
                        C = effective_c,
                        solver = "lbfgs",
                        max_iter = self.parameters["max_iter"],
                        tol = 1e-08,
                        fit_intercept = True,
                        random_state = SEED,
                    )
                else:
                    model = HistGradientBoostingClassifier(
                        max_leaf_nodes = self.parameters["leaves"],
                        max_iter = self.parameters["iterations"],
                        min_samples_leaf = self.parameters["min_leaf"],
                        learning_rate = 0.06,
                        l2_regularization = 1.0,
                        early_stopping = False,
                        random_state = SEED,
                    )
                model.fit(transformed, local_y.astype(int))
            self.experts_[index] = {
                "scaler": scaler,
                "model": model,
                "constant_probability": constant,
            }
            self.diagnostics_["experts"][name] = {
                "training_rows": len(selected),
                "training_prevalence": float(local_y.mean()),
                "training_positions": selected.tolist(),
                "status": "CONSTANT_TRAINING_CLASS"
                if (constant is not None)
                else "FITTED",
                "effective_C": effective_c,
                "early_stopping": False,
                "method_indicators_used_for_routing_only": True,
            }
        return self

    def predict_proba(self, x):
        import numpy as np

        require(hasattr(self, "experts_"), "Fitted method experts required")
        values = arrays(x, self.continuous_count)
        methods = method_indices(values)
        positive = np.empty(len(values))
        for (index) in (np.unique(methods)):
            require(
                int(index) in self.experts_,
                "Prediction method absent from training; no pooled fallback allowed",
            )
            selected = np.flatnonzero(methods == index)
            expert = self.experts_[int(index)]
            if (expert["constant_probability"] is not None):
                positive[selected] = expert["constant_probability"]
            else:
                transformed = expert["scaler"].transform(
                    values[selected, : self.continuous_count]
                )
                model = expert["model"]
                positive[selected] = model.predict_proba(transformed)[
                    :, list(model.classes_).index(1)
                ]
        require(
            np.isfinite(positive).all() and ((positive >= 0) & (positive <= 1)).all(),
            "Invalid expert probability",
        )
        return np.column_stack((1 - positive, positive))

def build_estimator(spec, training_rows, continuous_count):
    require(
        spec in specifications(),
        "Only the exact frozen eight-specification grid is supported",
    )
    require(
        type(training_rows) is int and training_rows > 0,
        "Positive complete training count required",
    )
    if (spec["family"] == "composite_pairwise_logistic"):
        return CompositePairwiseLogistic(
            continuous_count, training_rows, **spec["params"]
        )
    return MethodExperts(
        spec["family"], spec["params"], continuous_count, training_rows
    )
