import json
import math
import time
import warnings
from fractions import Fraction
from ..evaluation.prediction import (
    require,
    canonical_identity,
    task_metrics,
    unavailable_predictions,
)
from . import models_r4, models_r6, fit_r7

def safe_diagnostics(value):
    try:
        json.dumps(value, allow_nan = False)
        return value
    except (TypeError, ValueError, OverflowError):
        if (isinstance(value, dict)):
            return {str(key): safe_diagnostics(item) for (key, item) in (value.items())}
        if (isinstance(value, (list, tuple))):
            return [safe_diagnostics(item) for (item) in (value)]
        if (isinstance(value, float) and not math.isfinite(value)):
            return {
                "nonfinite_diagnostic_repr": repr(value),
                "type": type(value).__name__,
            }
        if (hasattr(value, "tolist")):
            return safe_diagnostics(value.tolist())
        if (hasattr(value, "item")):
            return safe_diagnostics(value.item())
        return {"non_json_diagnostic_repr": repr(value), "type": type(value).__name__}

def auxiliary_class(row):
    require(
        row["target"] in (0, 1) and row["availability"] == "AVAILABLE",
        "Only available training targets can supply auxiliary classes",
    )
    raw = row["event_details"]["signed_mean_delta"]
    require(
        isinstance(raw, dict)
        and type(raw["numerator"]) is int
        and type(raw["denominator"]) is int
        and raw["denominator"] > 0,
        "Exact signed rational numerator and positive denominator required",
    )
    delta = Fraction(raw["numerator"], raw["denominator"])
    label = 0 if (delta <= -1) else 2 if (delta >= 1) else 1
    require(
        int(label != 1) == row["target"],
        "Auxiliary class must collapse to immutable binary target",
    )
    return label

class Fitter:
    def __init__(self, specs, partition_order = None, data_id = None):
        self.specs = specs
        self.partition_order = partition_order
        self.data_id = data_id
        self.last_failure = None

    def fit_model(self, spec, track, usable, vectors, cache):
        import numpy as np
        from sklearn.exceptions import ConvergenceWarning
        from threadpoolctl import threadpool_limits

        specs = self.specs
        auxiliary = (
            [auxiliary_class(row) for (row) in (usable)] if (specs.round_id == "R4") else None
        )
        context = {
            "round": specs.round_id,
            "specification": spec,
            "track": track,
            "training": [
                [row["row_id"], row["target"], vectors[row["row_id"]]] for (row) in (usable)
            ],
            "auxiliary": auxiliary,
            "partition_order": self.partition_order,
        }
        key = canonical_identity(context)
        if (key in cache):
            saved = cache[key]
            self.last_failure = saved if (saved["status"] != "COMPLETE") else None
            require(
                saved["status"] == "COMPLETE",
                "Preserved failed fit: " + str(saved.get("failure")),
            )
            return saved["model"], saved["diagnostics"], True
        model, diagnostics = None, {}
        try:
            x = np.asarray([vectors[row["row_id"]] for (row) in (usable)])
            y = np.asarray([row["target"] for (row) in (usable)], dtype = int)
            with (
                threadpool_limits(limits = 1),
                warnings.catch_warnings(record = True) as caught,
            ):
                warnings.simplefilter("always")
                if (specs.round_id == "R6"):
                    model, diagnostics = models_r6.fit_calibrated(
                        spec, track, usable, vectors, self.partition_order
                    )
                elif (specs.round_id == "R4"):
                    model = specs.build_estimator(
                        spec, len(usable), len(specs.TRACKS[track])
                    )
                    model.fit(x, np.asarray(auxiliary, dtype = int))
                    diagnostics = {
                        "auxiliary_classes": auxiliary,
                        "auxiliary_training_row_ids": [row["row_id"] for (row) in (usable)],
                        "signed_response": model.diagnostics_,
                    }
                else:
                    model = specs.build_estimator(
                        spec, len(usable), len(specs.TRACKS[track])
                    )
                    if (model is None or len(set(y)) == 1):
                        model = {"constant_probability": float(y.mean())}
                    else:
                        model.fit(x, y)
                        diagnostics = getattr(model, "diagnostics_", {})
                messages = [
                    {"category": item.category.__name__, "message": str(item.message)}
                    for (item) in (caught)
                ]
                if (any(
                    issubclass(item.category, ConvergenceWarning) for (item) in (caught)
                )):
                    raise RuntimeError(
                        "Preserved convergence warning: " + str(messages)
                    )
            diagnostics = safe_diagnostics({"model": diagnostics, "warnings": messages})
            cache[key] = {
                "status": "COMPLETE",
                "model": model,
                "diagnostics": diagnostics,
            }
            return model, diagnostics, False
        except Exception as error:
            cache[key] = {
                "status": "FAILED",
                "failure": {"type": type(error).__name__, "message": str(error)},
                "diagnostics": safe_diagnostics(
                    getattr(model, "diagnostics_", diagnostics)
                ),
            }
            self.last_failure = cache[key]
            raise

    def __call__(self, index, task, spec, train_ids, test_ids, cache):
        import numpy as np
        from threadpoolctl import threadpool_limits

        started, cpu = time.monotonic(), time.process_time()
        self.last_failure = None
        specs = self.specs
        train, test = (
            [index[key] for (key) in (train_ids)],
            [index[key] for (key) in (test_ids)],
        )
        vectors = {
            row["row_id"]: specs.vector(row, task["track"]) for (row) in (train + test)
        }
        usable = [
            row
            for (row) in (train)
            if (row["target"] is not None and vectors[row["row_id"]] is not None)
        ]
        predictable = [row for (row) in (test) if (vectors[row["row_id"]] is not None)]
        metadata = {
            "specification": spec,
            "feature_names": specs.TRACKS[task["track"]]
            + ["method_DawnRank", "method_PersonaDrive"],
            "preprocessing_fit_row_ids": [row["row_id"] for (row) in (usable)],
            "training_only_preprocessing": True,
            "threads": 1,
            "random_seed": specs.SEED,
            "tuning": "INNER_FIXED_GRID"
            if (task["phase"] == "INNER")
            else "INNER_ONLY_NESTED_SELECTION_NO_OUTER_LABEL_TUNING",
        }
        result = {
            **task,
            "spec_id": spec["spec_id"],
            "train_row_ids": train_ids,
            "test_row_ids": test_ids,
            "train_usable_row_ids": [row["row_id"] for (row) in (usable)],
            "fit_metadata": metadata,
            "warnings": [],
            "failure": None,
        }
        try:
            require(
                spec in specs.specifications() and usable,
                "Exact finite candidate and nonempty usable training population required",
            )
            if (specs.round_id == "R7"):
                require(predictable, "Nonempty predictable population required")
                probabilities, diagnostics = fit_r7.fit(
                    usable,
                    predictable,
                    task["track"],
                    spec,
                    self.partition_order,
                    cache,
                    self.data_id,
                )
                metadata["diagnostics"] = safe_diagnostics(diagnostics)
            else:
                model, diagnostics, reused = self.fit_model(
                    spec, task["track"], usable, vectors, cache
                )
                metadata.update(diagnostics = diagnostics, fit_cache_hit = reused)
                matrix = np.asarray(
                    [vectors[row["row_id"]] for (row) in (predictable)]
                ).reshape(len(predictable), len(specs.TRACKS[task["track"]]) + 2)
                with threadpool_limits(limits = 1):
                    probabilities = (
                        np.full(len(predictable), model["constant_probability"])
                        if (isinstance(model, dict))
                        else model.predict_proba(matrix)[
                            :, list(model.classes_).index(1)
                        ]
                        if (predictable)
                        else np.asarray([])
                    )
                    if (specs.round_id == "R4"):
                        three = model.predict_three_class(matrix)
                        metadata["three_class_probabilities"] = [
                            {"row_id": row["row_id"], "probabilities": values.tolist()}
                            for (row, values) in (zip(predictable, three))
                        ]
                    if (specs.round_id == "R6"):
                        metadata["prediction_components"] = model.predict_components(
                            matrix
                        )
            values = np.asarray(probabilities)
            require(
                values.shape == (len(predictable),)
                and np.isfinite(values).all()
                and ((values >= 0) & (values <= 1)).all(),
                "Invalid fitted probability vector",
            )
            pmap = {
                row["row_id"]: float(probability)
                for (row, probability) in (zip(predictable, values))
            }
            result.update(
                status = "COMPLETE",
                predictions = [
                    {
                        "row_id": row["row_id"],
                        "status": "AVAILABLE"
                        if (row["row_id"] in pmap)
                        else "UNAVAILABLE",
                        "probability": pmap.get(row["row_id"]),
                        "reason": None
                        if (row["row_id"] in pmap)
                        else "MODEL_INPUT_UNAVAILABLE",
                    }
                    for (row) in (test)
                ],
            )
        except Exception as error:
            result.update(
                status = "FAILED",
                predictions = unavailable_predictions(test, "FAILED_TRIAL"),
                failure = {"type": type(error).__name__, "message": str(error)},
            )
            if (self.last_failure is not None):
                metadata["preserved_fit_failure"] = safe_diagnostics(self.last_failure)
        result["metrics"] = task_metrics(test, result["predictions"])
        result["timing"] = {
            "elapsed_seconds": time.monotonic() - started,
            "process_cpu_seconds": time.process_time() - cpu,
        }
        return result
