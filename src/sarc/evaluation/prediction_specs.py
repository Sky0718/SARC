import math
from fractions import Fraction

SEED = 20260930
SELECTION_RULE = "MAX_EQUAL_INNER_CONTEXT_CAPTURE25_THEN_MIN_BRIER_THEN_SPEC_ID"
FEATURES = (
    "log_edit_count",
    "mutated_edit_fraction",
    "deg_edit_fraction",
    "between_module_edit_fraction",
    "mean_endpoint_log_degree",
    "mean_consumed_cost_change",
)
TRACKS = {
    "strict_input_only": list(FEATURES),
    "edit_volume_only": [FEATURES[0]],
    "degree_incidence_only": [FEATURES[1], FEATURES[4]],
    "additive_volume_degree": [FEATURES[0], FEATURES[1], FEATURES[4]],
    "old_output_margin": ["old_native_rank_margin"],
    "input_plus_margin": [*FEATURES, "old_native_rank_margin"],
}

def specifications():
    values = [("prevalence", "prevalence", {})]
    values += [
        ("ridge_" + str(value), "ridge_logistic", {"penalty": value})
        for (value) in ((0.003, 0.01, 0.03, 0.3, 3.0, 30.0))
    ]
    values += [
        ("quadratic_" + str(value), "quadratic_logistic", {"penalty": value})
        for (value) in ((0.03, 0.3, 3.0))
    ]
    values += [
        ("spline_3", "spline_logistic", {"knots": 3, "penalty": 0.03}),
        ("spline_5", "spline_logistic", {"knots": 5, "penalty": 0.3}),
    ]
    values += [
        (
            "hist_7",
            "hist_gradient_boosting",
            {"leaves": 7, "iterations": 150, "min_leaf": 64},
        ),
        (
            "hist_15",
            "hist_gradient_boosting",
            {"leaves": 15, "iterations": 150, "min_leaf": 128},
        ),
    ]
    values += [
        (
            family + "_" + str(depth),
            family,
            {"depth": depth, "min_leaf": 32 if (depth == 6) else 64, "trees": 96},
        )
        for (family) in (("extra_trees", "random_forest"))
        for (depth) in ((6, 10))
    ]
    values += [
        (
            "rbf_" + str(gamma),
            "approximate_rbf_logistic",
            {"gamma": gamma, "components": 128, "penalty": 0.03},
        )
        for (gamma) in ((0.1, 1.0))
    ]
    values += [
        ("xgb_" + str(depth), "xgboost", {"depth": depth, "iterations": 150})
        for (depth) in ((2, 4))
    ]
    values += [
        ("cat_" + str(depth), "catboost", {"depth": depth, "iterations": 150})
        for (depth) in ((3, 5))
    ]
    return [
        {
            "spec_id": identity,
            "family": family,
            "params": params,
            "complexity_order": index,
        }
        for ((index, (identity, family, params))) in (enumerate(values))
    ]

def number(value):
    if (isinstance(value, dict) and {"numerator", "denominator"} <= set(value)):
        value = Fraction(value["numerator"], value["denominator"])
    value = float(value)
    if (not math.isfinite(value)):
        raise ValueError("Nonfinite feature")
    return value

def vector(row, track):
    values = []
    for (feature) in (TRACKS[track]):
        value = row["features"].get(feature)
        if (feature == "mean_consumed_cost_change" and row["method"] != "PRODIGY"):
            if (
                value is not None
                or row["features"].get("cost_structurally_not_applicable") is not True
            ):
                raise ValueError("Binary method structural cost contract differs")
            value = 0.0
        if (value is None):
            return None
        values.append(number(value))
    return values + [
        float(row["method"] == "DawnRank"),
        float(row["method"] == "PersonaDrive"),
    ]

def build_estimator(spec, training_rows, continuous_count):
    from sklearn.compose import ColumnTransformer
    from sklearn.ensemble import (
        ExtraTreesClassifier,
        HistGradientBoostingClassifier,
        RandomForestClassifier,
    )
    from sklearn.kernel_approximation import RBFSampler
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import (
        PolynomialFeatures,
        SplineTransformer,
        StandardScaler,
    )

    family, parameters = (spec["family"], spec["params"])
    if (family == "prevalence"):
        return None
    if (family.endswith("logistic")):
        continuous = list(range(continuous_count))
        method_columns = [continuous_count, continuous_count + 1]
        steps = [
            (
                "initial_scale",
                ColumnTransformer(
                    [
                        ("continuous", StandardScaler(), continuous),
                        ("method", "passthrough", method_columns),
                    ],
                    remainder = "drop",
                ),
            )
        ]
        if (family == "quadratic_logistic"):
            steps.extend(
                [
                    ("interactions", PolynomialFeatures(degree = 2, include_bias = False)),
                    ("interaction_scale", StandardScaler()),
                ]
            )
        elif (family == "spline_logistic"):
            steps.append(
                (
                    "splines",
                    ColumnTransformer(
                        [
                            (
                                "continuous",
                                SplineTransformer(
                                    n_knots = parameters["knots"],
                                    degree = 3,
                                    include_bias = False,
                                    extrapolation = "linear",
                                ),
                                continuous,
                            ),
                            ("method", "passthrough", method_columns),
                        ],
                        remainder = "drop",
                    ),
                )
            )
        elif (family == "approximate_rbf_logistic"):
            steps.append(
                (
                    "rbf",
                    RBFSampler(
                        gamma = parameters["gamma"],
                        n_components = parameters["components"],
                        random_state = SEED,
                    ),
                )
            )
        steps.append(
            (
                "model",
                LogisticRegression(
                    C = 1.0 / (training_rows * parameters["penalty"]),
                    solver = "lbfgs",
                    max_iter = 1000,
                    tol = 1e-08,
                    fit_intercept = True,
                    random_state = SEED,
                ),
            )
        )
        return Pipeline(steps)
    if (family == "hist_gradient_boosting"):
        return HistGradientBoostingClassifier(
            max_leaf_nodes = parameters["leaves"],
            max_iter = parameters["iterations"],
            min_samples_leaf = parameters["min_leaf"],
            learning_rate = 0.06,
            l2_regularization = 1.0,
            early_stopping = False,
            random_state = SEED,
        )
    if (family in ("extra_trees", "random_forest")):
        cls = (
            ExtraTreesClassifier
            if (family == "extra_trees")
            else RandomForestClassifier
        )
        return cls(
            n_estimators = parameters["trees"],
            max_depth = parameters["depth"],
            min_samples_leaf = parameters["min_leaf"],
            max_features = 1.0,
            bootstrap = family == "random_forest",
            class_weight = None,
            n_jobs = 1,
            random_state = SEED,
        )
    if (family == "xgboost"):
        from xgboost import XGBClassifier

        return XGBClassifier(
            n_estimators = parameters["iterations"],
            max_depth = parameters["depth"],
            learning_rate = 0.05,
            min_child_weight = 20,
            reg_lambda = 10,
            reg_alpha = 0,
            subsample = 1,
            colsample_bytree = 1,
            tree_method = "hist",
            objective = "binary:logistic",
            eval_metric = "logloss",
            n_jobs = 1,
            random_state = SEED,
            verbosity = 0,
        )
    if (family == "catboost"):
        from catboost import CatBoostClassifier

        return CatBoostClassifier(
            iterations = parameters["iterations"],
            depth = parameters["depth"],
            learning_rate = 0.05,
            l2_leaf_reg = 10,
            loss_function = "Logloss",
            random_seed = SEED,
            thread_count = 1,
            verbose = False,
            allow_writing_files = False,
            bootstrap_type = "No",
            random_strength = 0,
            use_best_model = False,
        )
    raise ValueError("Unfrozen family")
