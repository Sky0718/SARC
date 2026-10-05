from . import specs_r1 as original
from .mechanism import FEATURE_NAMES

SEED = original.SEED
FEATURES = original.FEATURES
METHODS = original.METHODS
MASTERS = original.MASTERS
COHORTS = original.COHORTS
TRANSITIONS = original.TRANSITIONS
EXPECTED_RUNTIME = original.EXPECTED_RUNTIME
EXPECTED_PYTHON = original.EXPECTED_PYTHON
SELECTION_RULE = original.SELECTION_RULE
number = original.number
LOCALISATION = (
    "added_mutant_contact_edges",
    "removed_mutant_contact_edges",
    "mutant_candidates_touched_fraction",
    "mean_mutant_relative_edit_burden",
    "max_mutant_relative_edit_burden",
    "mean_mutant_deg_onehop_edit_burden",
)
STRICT = [*FEATURES, *LOCALISATION]
TRACKS = {
    "mechanism_only": list(FEATURE_NAMES),
    "strict_localisation_mechanism": [*STRICT, *FEATURE_NAMES],
    "strict_localisation_direct_support": [
        *STRICT,
        *FEATURE_NAMES[:2],
        *FEATURE_NAMES[10:],
    ],
    "strict_localisation_reweight_interaction_support": [
        *STRICT,
        *FEATURE_NAMES[2:6],
        *FEATURE_NAMES[8:],
    ],
}

def vector(row, track):
    if (row["method"] not in METHODS or track not in TRACKS):
        raise ValueError("Unknown method or undeclared feature track")
    values = []
    for (name) in (TRACKS[track]):
        value = row["features"].get(name)
        if (name == "mean_consumed_cost_change" and row["method"] != "PRODIGY"):
            if (
                value is not None
                or row["features"].get("cost_structurally_not_applicable") is not True
            ):
                raise ValueError("Unchanged structural cost contract required")
            value = 0.0
        if (value is None):
            return None
        values.append(number(value))
    return values + [
        float(row["method"] == "DawnRank"),
        float(row["method"] == "PersonaDrive"),
    ]

def specifications():
    values = [
        ("r5_ridge_" + str(value), "ridge_logistic", {"penalty": value})
        for (value) in ((0.03, 3.0))
    ]
    values += [
        (
            "r5_extra_trees_" + str(value),
            "extra_trees",
            {"trees": 500, "min_leaf": value},
        )
        for (value) in ((16, 64))
    ]
    values += [
        (
            "r5_hist_" + str(leaves),
            "hist_gradient_boosting",
            {"leaves": leaves, "iterations": 150, "min_leaf": minimum},
        )
        for (leaves, minimum) in (((7, 64), (15, 128)))
    ]
    return [
        {
            "spec_id": identity,
            "family": family,
            "params": parameters,
            "complexity_order": index,
        }
        for (index, (identity, family, parameters)) in (enumerate(values))
    ]

def build_estimator(spec, training_rows, continuous_count):
    from sklearn.compose import ColumnTransformer
    from sklearn.ensemble import ExtraTreesClassifier, HistGradientBoostingClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    if (
        spec not in specifications()
        or training_rows < 1
        or continuous_count not in {len(values) for (values) in (TRACKS.values())}
    ):
        raise ValueError(
            "Exact new finite specification and valid training dimensions required"
        )
    parameters = spec["params"]
    if (spec["family"] == "ridge_logistic"):
        transform = ColumnTransformer(
            [
                ("continuous", StandardScaler(), list(range(continuous_count))),
                ("method", "passthrough", [continuous_count, continuous_count + 1]),
            ],
            remainder = "drop",
        )
        model = LogisticRegression(
            C = 1.0 / (training_rows * parameters["penalty"]),
            solver = "lbfgs",
            max_iter = 3000,
            tol = 1e-08,
            fit_intercept = True,
            random_state = SEED,
        )
        return Pipeline([("initial_scale", transform), ("model", model)])
    if (spec["family"] == "extra_trees"):
        return ExtraTreesClassifier(
            n_estimators = 500,
            max_depth = None,
            min_samples_leaf = parameters["min_leaf"],
            max_features = 1.0,
            bootstrap = False,
            class_weight = None,
            n_jobs = 1,
            random_state = SEED,
        )
    return HistGradientBoostingClassifier(
        max_leaf_nodes = parameters["leaves"],
        max_iter = 150,
        min_samples_leaf = parameters["min_leaf"],
        learning_rate = 0.06,
        l2_regularization = 1.0,
        early_stopping = False,
        random_state = SEED,
    )
