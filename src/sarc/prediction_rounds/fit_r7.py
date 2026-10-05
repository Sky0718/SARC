import copy
import json
import time
import warnings
from . import specs_r1 as original
from . import models_r7 as residual
from ..evaluation import prediction as base

TRACKS = residual.TRACKS

def baseline_spec(depth):
    return next(
        (
            spec
            for (spec) in (original.specifications())
            if (spec["spec_id"] == "extra_trees_" + str(depth))
        )
    )

def specifications():
    result = []
    for (depth) in ((6, 10)):
        for (ridge) in ((0.03, 0.3)):
            for (budget) in ((0.0, 0.5)):
                result.append(
                    {
                        "spec_id": "r7_depth_"
                        + str(depth)
                        + "_ridge_"
                        + str(ridge)
                        + "_budget_"
                        + str(budget),
                        "family": "residual_budget" if (budget) else "residual_bce",
                        "params": {
                            "depth": depth,
                            "ridge": ridge,
                            "budget_weight": budget,
                        },
                    }
                )
        result.append(
            {
                "spec_id": "r7_anchor_depth_" + str(depth),
                "family": "anchor",
                "params": {"depth": depth},
            }
        )
        result.append(
            {
                "spec_id": "r7_same_feature_extra_trees_" + str(depth),
                "family": "same_feature_extra_trees",
                "params": {"depth": depth},
            }
        )
    return [{**spec, "complexity_order": index} for (index, spec) in (enumerate(result))]

def feature_vector(row, names):
    values = []
    for (name) in (names):
        value = row["features"].get(name)
        if (name == "mean_consumed_cost_change" and row["method"] != "PRODIGY"):
            base.require(
                value is None
                and row["features"].get("cost_structurally_not_applicable") is True,
                "Original binary-method cost semantics changed",
            )
            value = 0.0
        if (value is None):
            return None
        values.append(original.number(value))
    return values + [
        float(row["method"] == "DawnRank"),
        float(row["method"] == "PersonaDrive"),
    ]

def vector(row, track):
    return feature_vector(
        row, TRACKS["old_inputs"] if (track == "input_plus_margin") else TRACKS[track]
    )

def predict_tree(estimator, rows, names):
    import numpy as np
    from threadpoolctl import threadpool_limits

    vectors = [feature_vector(row, names) for (row) in (rows)]
    base.require(
        all((value is not None for (value) in (vectors))),
        "Predictable full feature support required",
    )
    with threadpool_limits(limits = 1):
        probabilities = (
            np.full(len(rows), estimator["constant_probability"])
            if (isinstance(estimator, dict))
            else estimator.predict_proba(np.asarray(vectors))[
                :, list(estimator.classes_).index(1)
            ]
        )
    base.require(
        np.isfinite(probabilities).all()
        and ((probabilities >= 0) & (probabilities <= 1)).all(),
        "Invalid tree component probabilities",
    )
    return probabilities.tolist()

def predict_saved_residual(diagnostics, testing, matrix, test_anchor):
    import numpy as np
    from scipy.special import expit, logit

    count = len(diagnostics["specification"]["feature_names"])
    array = np.asarray(matrix, dtype = float)
    scale = diagnostics["scaler"]
    continuous = (array[:, :count] - np.asarray(scale["mean"])) / np.asarray(
        scale["scale"]
    )
    bits = array[:, count:]
    design = np.column_stack(
        (
            np.ones(len(testing)),
            continuous,
            bits,
            continuous * bits[:, [0]],
            continuous * bits[:, [1]],
        )
    )
    probabilities = expit(
        logit(
            np.clip(
                test_anchor["scores"],
                residual.PROBABILITY_CLIP,
                1.0 - residual.PROBABILITY_CLIP,
            )
        )
        + design @ np.asarray(diagnostics["coefficients"])
    )
    base.require(
        np.isfinite(probabilities).all(), "Invalid cached residual probabilities"
    )
    return probabilities.tolist()

def cache_identity(kind, training, names, spec, data_id, partition_order = None):
    residual.row_ids(training)
    base.require(
        isinstance(data_id, str) and data_id,
        "Explicit training-carrier identity required",
    )
    vectors = [feature_vector(row, names) for (row) in (training)]
    base.require(
        all(value is not None for (value) in (vectors)),
        "Complete current training features required",
    )
    base.require(
        all(
            type(row["target"]) in (int, float) and row["target"] in (0, 1)
            for (row) in (training)
        ),
        "Complete binary training labels required",
    )
    value = {
        "kind": kind,
        "data_id": data_id,
        "specification": spec,
        "feature_names": list(names),
        "training": [
            {
                "row_id": row["row_id"],
                "cohort": row["cohort"],
                "sample_id": row["sample_id"],
                "method": row["method"],
                "target": row["target"],
                "vector": vector,
            }
            for (row, vector) in (zip(training, vectors))
        ],
    }
    if (partition_order is not None):
        value["partition_order"] = partition_order
    return (
        "R7",
        json.dumps(value, sort_keys = True, separators = (",", ":"), allow_nan = False),
    )

def cached_record(cache, key):
    if (key not in cache):
        return None
    record = cache[key]
    base.require(
        isinstance(record, dict) and record.get("status") in ("COMPLETE", "FAILED"),
        "Invalid R7 cache record",
    )
    if (record["status"] == "FAILED"):
        raise RuntimeError(
            "Preserved R7 fit failure: " + json.dumps(record["failure"], sort_keys = True)
        )
    return record

def failure_record(error, started, observed_warnings):
    return {
        "status": "FAILED",
        "failure": {"type": type(error).__name__, "message": str(error)},
        "elapsed_seconds": time.monotonic() - started,
        "warnings": [
            {"category": warning.category.__name__, "message": str(warning.message)}
            for (warning) in (observed_warnings)
        ],
    }

def cached_estimator(training, names, spec, cache, data_id):
    import numpy as np
    from threadpoolctl import threadpool_limits

    base.require(
        spec == baseline_spec(spec["params"]["depth"]),
        "Exact 96-tree R2 specification required",
    )
    key = cache_identity("TREE_COMPONENT", training, names, spec, data_id)
    stored = cached_record(cache, key)
    if (stored is not None):
        return stored["estimator"], {
            **copy.deepcopy(stored["diagnostics"]),
            "fit_cache_hit": True,
            "actual_new_fit": False,
        }
    started = time.monotonic()
    observed = []
    try:
        vectors = [feature_vector(row, names) for (row) in (training)]
        labels = np.asarray([row["target"] for (row) in (training)], dtype = int)
        with (
            warnings.catch_warnings(record = True) as observed,
            threadpool_limits(limits = 1),
        ):
            warnings.simplefilter("always")
            estimator = (
                {"constant_probability": float(labels.mean())}
                if (len(set(labels)) == 1)
                else original.build_estimator(spec, len(training), len(names))
            )
            if (not isinstance(estimator, dict)):
                estimator.fit(np.asarray(vectors), labels)
        diagnostics = {
            "status": "COMPLETE",
            "specification": spec,
            "feature_names": list(names),
            "training_row_ids": residual.row_ids(training),
            "training_rows": len(training),
            "data_id": data_id,
            "random_seed": original.SEED,
            "threads": 1,
            "elapsed_seconds": time.monotonic() - started,
            "warnings": [
                {"category": warning.category.__name__, "message": str(warning.message)}
                for (warning) in (observed)
            ],
            "fit_cache_hit": False,
            "actual_new_fit": True,
        }
        cache[key] = {
            "status": "COMPLETE",
            "estimator": estimator,
            "diagnostics": copy.deepcopy(diagnostics),
        }
        return estimator, diagnostics
    except Exception as error:
        cache[key] = failure_record(error, started, observed)
        raise

def anchor_identity(spec, data_id):
    return {
        "baseline_spec": spec,
        "feature_names": [
            *residual.LOCALISATION,
            "method_DawnRank",
            "method_PersonaDrive",
        ],
        "source_id": "R2_LOCALISATION_EXTRA_TREES_96_SEED_20260930",
        "data_id": data_id,
    }

def anchor_predict(training, testing, identity, cache, partition_order):
    residual.validate_identity(identity)
    estimator, receipt = cached_estimator(
        training,
        residual.LOCALISATION,
        identity["baseline_spec"],
        cache,
        identity["data_id"],
    )
    package = {
        "scores": predict_tree(estimator, testing, residual.LOCALISATION),
        "row_ids": residual.row_ids(testing),
        "training_row_ids": residual.row_ids(training),
        "identity": copy.deepcopy(identity),
        "mode": "outer_training_fit",
        "fit_receipt": receipt,
    }
    residual.validate_anchor(package, testing, training, False, partition_order)
    return package

def cached_residual(
    training, testing, track, spec, test_anchor, partition_order, cache, data_id
):
    from threadpoolctl import threadpool_limits

    base.require(track in TRACKS, "Frozen R7 feature track required")
    key = cache_identity(
        "RESIDUAL_COMPONENT", training, TRACKS[track], spec, data_id, partition_order
    )
    matrix = [vector(row, track) for (row) in (testing)]
    base.require(
        all(value is not None for (value) in (matrix)),
        "Complete prediction feature support required",
    )
    stored = cached_record(cache, key)
    if (stored is not None):
        diagnostics = copy.deepcopy(stored["diagnostics"])
        base.require(
            diagnostics["anchor_identity"] == test_anchor["identity"],
            "Cached residual anchor identity differs",
        )
        predictions = predict_saved_residual(diagnostics, testing, matrix, test_anchor)
        diagnostics.update(test_rows = len(testing), fit_cache_hit = True)
        return predictions, diagnostics
    started = time.monotonic()
    observed = []
    try:
        callback = lambda fitrows, predictrows, identity: anchor_predict(
            fitrows, predictrows, identity, cache, partition_order
        )
        training_anchor = residual.build_crossfit_anchors(
            training, test_anchor["identity"], callback, partition_order
        )
        parameters = spec["params"]
        residual_spec = next(
            item
            for (item) in (residual.specifications(TRACKS[track]))
            if (item["alpha"] == 1
            and item["ridge"] == parameters["ridge"]
            and item["budget_weight"] == parameters["budget_weight"])
        )
        with (
            warnings.catch_warnings(record = True) as observed,
            threadpool_limits(limits = 1),
        ):
            warnings.simplefilter("always")
            result = residual.fit_predict(
                training,
                testing,
                [vector(row, track) for (row) in (training)],
                matrix,
                [row["target"] for (row) in (training)],
                training_anchor,
                test_anchor,
                residual_spec,
                partition_order,
            )
        diagnostics = {
            **result["diagnostics"],
            "training_anchor_fit_receipts": training_anchor["fit_receipts"],
            "optimiser_configuration": {
                "method": "L-BFGS-B",
                "maxiter": 300,
                "ftol": 1e-12,
                "gtol": 1e-7,
                "maxls": 50,
                "initial": "ZERO",
                "bounds": [-10, 10],
                "ridge_includes_intercept": True,
            },
            "elapsed_seconds": time.monotonic() - started,
            "warnings": [
                {"category": warning.category.__name__, "message": str(warning.message)}
                for (warning) in (observed)
            ],
            "fit_cache_hit": False,
        }
        cache[key] = {"status": "COMPLETE", "diagnostics": copy.deepcopy(diagnostics)}
        return result["score"], diagnostics
    except Exception as error:
        cache[key] = failure_record(error, started, observed)
        raise

def fit(training, testing, track, spec, partition_order, cache, data_id):
    base.require(isinstance(cache, dict), "Shared in-memory training cache required")
    base.require(
        spec in specifications() and track in TRACKS,
        "Exact finite R7 specification and feature track required",
    )
    residual.row_ids(training)
    residual.row_ids(testing)
    base.require(
        not set(residual.row_ids(training)).intersection(residual.row_ids(testing)),
        "Training and testing row identities overlap",
    )
    base.require(
        all(vector(row, track) is not None for (row) in (training + testing)),
        "Use complete declared feature support",
    )
    tree_spec = baseline_spec(spec["params"]["depth"])
    diagnostics = {
        "specification": copy.deepcopy(spec),
        "feature_names": [*TRACKS[track], "method_DawnRank", "method_PersonaDrive"],
        "preprocessing_fit_row_ids": residual.row_ids(training),
        "training_only_preprocessing": True,
        "threads": 1,
        "random_seed": original.SEED,
    }
    if (spec["family"] == "same_feature_extra_trees"):
        estimator, component = cached_estimator(
            training, TRACKS[track], tree_spec, cache, data_id
        )
        probabilities = predict_tree(estimator, testing, TRACKS[track])
        diagnostics.update(
            component = component, fit_cache_hit = component["fit_cache_hit"]
        )
    else:
        test_anchor = anchor_predict(
            training,
            testing,
            anchor_identity(tree_spec, data_id),
            cache,
            partition_order,
        )
        diagnostics["evaluation_anchor"] = {
            key: value for (key, value) in (test_anchor.items()) if (key != "scores")
        }
        if (spec["family"] == "anchor"):
            probabilities = test_anchor["scores"]
            diagnostics.update(
                fit_cache_hit = test_anchor["fit_receipt"]["fit_cache_hit"],
                zero_residual_exact_identity = True,
            )
        else:
            probabilities, details = cached_residual(
                training,
                testing,
                track,
                spec,
                test_anchor,
                partition_order,
                cache,
                data_id,
            )
            diagnostics.update(
                fit_cache_hit = details["fit_cache_hit"], residual_diagnostics = details
            )
    return probabilities, diagnostics
