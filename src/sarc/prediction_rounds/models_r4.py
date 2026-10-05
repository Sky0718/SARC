import warnings

SEED = 20260930
METHODS = ("DawnRank", "PersonaDrive", "PRODIGY")
CLASS_MAPPING = (0, 1, 2)

def require(value, message):
    if (not value):
        raise ValueError(message)

def specifications():
    values = []
    for (penalty) in ((0.003, 0.03)):
        values.append(
            (
                "signed_multinomial_l" + str(penalty),
                "signed_multinomial_logistic",
                {"penalty": penalty, "max_iter": 1000, "tol": 1e-08},
            )
        )
    for (family, prefix) in ((
        ("signed_hist", "signed_hist"),
        ("method_expert_signed_hist", "expert_signed_hist"),
    )):
        for (leaves, minimum) in (((7, 64), (15, 128))):
            values.append(
                (
                    prefix + "_" + str(leaves),
                    family,
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

def class_from_signed_fraction(numerator, denominator):
    require(
        type(numerator) is int and type(denominator) is int and (denominator > 0),
        "Exact signed rational integers and positive denominator required",
    )
    return 0 if (numerator <= -denominator) else 2 if (numerator >= denominator) else 1

def arrays(x, continuous_count, auxiliary = None):
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
    if (auxiliary is None):
        return values
    labels = np.asarray(auxiliary)
    require(
        labels.ndim == 1 and len(labels) == len(values) and (len(labels) > 0),
        "Complete nonempty training auxiliary labels required",
    )
    require(
        labels.dtype.kind in ("i", "u") and np.isin(labels, CLASS_MAPPING).all(),
        "Training labels must be integer signed classes zero, one or two",
    )
    return (values, labels.astype(int))

def method_indices(x):
    import numpy as np

    return np.where(x[:, -2] == 1, 0, np.where(x[:, -1] == 1, 1, 2))

def class_counts(labels):
    import numpy as np

    return {
        str(category): int(np.count_nonzero(labels == category))
        for (category) in (CLASS_MAPPING)
    }

def validate_three_class(probabilities, rows):
    import numpy as np

    values = np.asarray(probabilities, dtype = float)
    require(values.shape == (rows, 3), "Exact three-class probability axis required")
    require(
        np.isfinite(values).all() and ((values >= 0) & (values <= 1)).all(),
        "Invalid signed-class probabilities",
    )
    require(
        np.allclose(values.sum(axis = 1), 1.0, rtol = 1e-12, atol = 1e-12),
        "Signed-class probabilities do not sum to one; repair forbidden",
    )
    return values

class SignedResponseEstimator:
    def __init__(self, spec, training_rows, continuous_count):
        self.spec = {**spec, "params": dict(spec["params"])}
        self.training_rows = training_rows
        self.continuous_count = continuous_count
        self._fit_complete = False

    def _fit_component(self, values, labels, positions, expert):
        import numpy as np
        from sklearn.ensemble import HistGradientBoostingClassifier
        from sklearn.exceptions import ConvergenceWarning
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler

        family, parameters = (self.spec["family"], self.spec["params"])
        observed = np.unique(labels).astype(int).tolist()
        diagnostic = {
            "training_rows": len(labels),
            "observed_classes": observed,
            "training_class_counts": class_counts(labels),
            "training_class_mapping": list(CLASS_MAPPING),
            "heldout_rows_used": False,
            "early_stopping": False,
            "effective_C": None,
            "continuous_scaling": "NONE",
            "library_class_mapping": observed,
            "warnings": [],
            "status": "FITTING",
        }
        if (expert):
            diagnostic.update(
                training_positions = positions.tolist(), no_pooled_fallback = True
            )
        scaler = None
        transformed = values[:, : self.continuous_count] if (expert) else values
        if (family == "signed_multinomial_logistic"):
            scaler = StandardScaler().fit(values[:, : self.continuous_count])
            transformed = np.column_stack(
                (scaler.transform(values[:, : self.continuous_count]), values[:, -2:])
            )
            diagnostic.update(
                continuous_scaling = "TRAINING_ONLY_STANDARD_SCALER_METHOD_BITS_PASSTHROUGH",
                scaler_training_rows = int(scaler.n_samples_seen_),
                effective_C = 1 / (len(labels) * parameters["penalty"]),
            )
        component = {
            "scaler": scaler,
            "model": None,
            "observed_classes": observed,
            "constant_class": observed[0] if (len(observed) == 1) else None,
            "expert": expert,
            "diagnostics": diagnostic,
        }
        if (component["constant_class"] is not None):
            diagnostic.update(
                status = "CONSTANT_TRAINING_CLASS", converged = True, iterations = 0
            )
            return component
        if (family == "signed_multinomial_logistic"):
            learner = LogisticRegression(
                C = diagnostic["effective_C"],
                solver = "lbfgs",
                max_iter = parameters["max_iter"],
                tol = parameters["tol"],
                fit_intercept = True,
                class_weight = None,
                random_state = SEED,
            )
        else:
            learner = HistGradientBoostingClassifier(
                max_leaf_nodes = parameters["leaves"],
                max_iter = parameters["iterations"],
                min_samples_leaf = parameters["min_leaf"],
                learning_rate = 0.06,
                l2_regularization = 1.0,
                early_stopping = False,
                random_state = SEED,
            )
        with warnings.catch_warnings(record = True) as caught:
            warnings.simplefilter("always")
            learner.fit(transformed, labels)
        diagnostic["warnings"] = [
            {"category": item.category.__name__, "message": str(item.message)}
            for (item) in (caught)
        ]
        if (any((issubclass(item.category, ConvergenceWarning) for (item) in (caught)))):
            diagnostic.update(status = "FAILED_CONVERGENCE_WARNING", converged = False)
            raise RuntimeError(
                "Preserved signed-response convergence warning: "
                + str(diagnostic["warnings"])
            )
        actual_classes = np.asarray(learner.classes_)
        require(
            actual_classes.dtype.kind in ("i", "u")
            and actual_classes.tolist() == observed,
            "Learner class axis differs from observed training classes",
        )
        if (family == "signed_multinomial_logistic"):
            require(
                np.isfinite(learner.coef_).all()
                and np.isfinite(learner.intercept_).all(),
                "Nonfinite fitted multinomial coefficients",
            )
            diagnostic.update(
                iterations = np.asarray(learner.n_iter_).astype(int).tolist(),
                optimisation = "LIBRARY_LBFGS_NO_CONVERGENCE_WARNING",
                converged = True,
            )
        else:
            require(
                int(learner.n_iter_) == parameters["iterations"],
                "Fixed histogram iteration count changed",
            )
            diagnostic.update(
                iterations = int(learner.n_iter_),
                optimisation = "FIXED_ITERATIONS_NO_EARLY_STOPPING",
                converged = None,
            )
        diagnostic.update(
            status = "FITTED", library_class_mapping = actual_classes.tolist()
        )
        component["model"] = learner
        return component

    def fit(self, x, auxiliary):
        import numpy as np

        self._fit_complete = False
        require(
            self.spec in specifications(),
            "Only exact frozen six-specification grid is supported",
        )
        values, labels = arrays(x, self.continuous_count, auxiliary)
        require(
            type(self.training_rows) is int and len(labels) == self.training_rows,
            "Declared usable training population differs",
        )
        expert = self.spec["family"] == "method_expert_signed_hist"
        self.classes_ = np.asarray([0, 1])
        self.three_classes_ = np.asarray(CLASS_MAPPING)
        self.n_features_in_ = values.shape[1]
        self.diagnostics_ = {
            "family": self.spec["family"],
            "training_rows": len(labels),
            "training_class_mapping": list(CLASS_MAPPING),
            "observed_classes": np.unique(labels).astype(int).tolist(),
            "training_class_counts": class_counts(labels),
            "heldout_rows_used": False,
            "early_stopping": False,
            "pooled_model_fitted": False,
            "probability_collapse": "[P_CLASS_1,P_CLASS_0_PLUS_P_CLASS_2]",
            "auxiliary_labels_are_predictors": False,
            "missing_training_classes_probability": 0,
            "status": "FITTING",
        }
        self.components_ = {}
        if (expert):
            self.diagnostics_.update(
                method_routing = "ONLY_DECLARED_FINAL_TWO_METHOD_INDICATORS",
                experts = {},
                no_pooled_fallback = True,
            )
            routes = method_indices(values)
            for (index, name) in (enumerate(METHODS)):
                positions = np.flatnonzero(routes == index)
                if (not len(positions)):
                    self.diagnostics_["experts"][name] = {
                        "training_rows": 0,
                        "training_positions": [],
                        "observed_classes": [],
                        "training_class_counts": {
                            str(category): 0 for (category) in (CLASS_MAPPING)
                        },
                        "training_class_mapping": list(CLASS_MAPPING),
                        "heldout_rows_used": False,
                        "early_stopping": False,
                        "no_pooled_fallback": True,
                        "status": "NO_TRAINING_ROWS_PREDICTION_FORBIDDEN",
                    }
                    continue
                component = self._fit_component(
                    values[positions], labels[positions], positions, True
                )
                self.components_[index] = component
                self.diagnostics_["experts"][name] = component["diagnostics"]
        else:
            component = self._fit_component(
                values, labels, np.arange(len(labels)), False
            )
            self.components_["pooled"] = component
            self.diagnostics_["component"] = component["diagnostics"]
            self.diagnostics_["pooled_model_fitted"] = component["model"] is not None
        self.diagnostics_["status"] = "COMPLETE"
        self._fit_complete = True
        return self

    def _component_probabilities(self, component, values):
        import numpy as np

        result = np.zeros((len(values), 3), dtype = float)
        if (component["constant_class"] is not None):
            result[:, component["constant_class"]] = 1.0
        elif (len(values)):
            transformed = (
                values[:, : self.continuous_count] if (component["expert"]) else values
            )
            if (component["scaler"] is not None):
                transformed = np.column_stack(
                    (
                        component["scaler"].transform(
                            values[:, : self.continuous_count]
                        ),
                        values[:, -2:],
                    )
                )
            model = component["model"]
            actual_classes = np.asarray(model.classes_)
            require(
                actual_classes.dtype.kind in ("i", "u")
                and actual_classes.tolist() == component["observed_classes"],
                "Prediction class mapping changed from training",
            )
            probabilities = np.asarray(model.predict_proba(transformed), dtype = float)
            require(
                probabilities.shape
                == (len(values), len(component["observed_classes"])),
                "Learner probability shape differs from class mapping",
            )
            result[:, actual_classes] = probabilities
        validate_three_class(result, len(values))
        absent = [
            category
            for (category) in (CLASS_MAPPING)
            if (category not in component["observed_classes"])
        ]
        require(
            not absent or np.all(result[:, absent] == 0),
            "Missing training class acquired nonzero probability",
        )
        return result

    def predict_three_class(self, x):
        import numpy as np

        require(
            self._fit_complete, "Complete fitted signed-response estimator required"
        )
        values = arrays(x, self.continuous_count)
        if ("pooled" in self.components_):
            return self._component_probabilities(self.components_["pooled"], values)
        result = np.zeros((len(values), 3), dtype = float)
        routes = method_indices(values)
        for (index) in (np.unique(routes)):
            require(
                int(index) in self.components_,
                "Prediction method absent from training; no pooled fallback allowed",
            )
            positions = np.flatnonzero(routes == index)
            result[positions] = self._component_probabilities(
                self.components_[int(index)], values[positions]
            )
        return validate_three_class(result, len(values))

    def predict_proba(self, x):
        import numpy as np

        three = self.predict_three_class(x)
        return np.column_stack((three[:, 1], three[:, 0] + three[:, 2]))

def build_estimator(spec, training_rows, continuous_count):
    require(
        spec in specifications(),
        "Only exact frozen six-specification grid is supported",
    )
    require(
        type(training_rows) is int and training_rows > 0,
        "Positive complete usable training count required",
    )
    require(
        type(continuous_count) is int and continuous_count > 0,
        "Positive declared continuous feature count required",
    )
    return SignedResponseEstimator(spec, training_rows, continuous_count)
