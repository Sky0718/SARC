import json
import math
from collections import Counter
from . import localisation as specs

def require(value, message):
    if (not value):
        raise ValueError(message)

def canonical_identity(value):
    return json.dumps(value, sort_keys = True, separators = (",", ":"), allow_nan = False)

def validate_rows(rows, production = True):
    require(
        len({row["row_id"] for (row) in (rows)}) == len(rows), "Duplicate logical rows"
    )
    if (production):
        require(len(rows) == 13824, "Full development carrier required")
    for (row) in (rows):
        require(
            row["cohort"] in specs.COHORTS
            and row["transition"] in specs.TRANSITIONS
            and (row["method"] in specs.METHODS)
            and (row["control_seed"] in (None, 271828, 314159)),
            "Protected or out-of-scope development row",
        )
        require(
            row["protected_values_opened"] is False
            and row["master_seeds"]
            == (list(specs.MASTERS) if (row["method"] == "PRODIGY") else [None]),
            "Protected access or altered algorithm seed axis",
        )
        require(
            row["target"] in (None, 0, 1)
            and (row["availability"] == "AVAILABLE") == (row["target"] is not None),
            "Target/availability differs",
        )
        require(
            row["row_id"]
            == "::".join(
                (row["cohort"], row["network_id"], row["method"], row["sample_id"])
            ),
            "Logical row identity differs",
        )
    if (production):
        for (cohort) in (specs.COHORTS):
            selected = [row for (row) in (rows) if (row["cohort"] == cohort)]
            require(
                len(selected) == 6912
                and len({row["sample_id"] for (row) in (selected)}) == 36
                and (len({row["network_id"] for (row) in (selected)}) == 64),
                "Complete biological or graph population differs",
            )
            require(
                set(Counter((row["sample_id"] for (row) in (selected))).values())
                == {192},
                "Unequal complete sample row coverage",
            )

def make_manifest(rows, production = True, families = None, tracks = None, schemes = None):
    validate_rows(rows, production)
    folds, inner_tasks, outer_tasks = ([], [], [])
    candidates = [
        item
        for (item) in (specs.specifications())
        if (families is None or item["family"] in families)
    ]
    families = sorted({item["family"] for (item) in (candidates)})
    for (axis, groups, inner_axis, inner_groups) in ((
        ("cohort", specs.COHORTS, "transition", specs.TRANSITIONS),
        ("transition", specs.TRANSITIONS, "cohort", specs.COHORTS),
    )):
        scheme = "LOCO" if (axis == "cohort") else "TRANSITION_EXPLORATORY"
        if (schemes is not None and scheme not in schemes):
            continue
        for (heldout) in (groups):
            fold_id = scheme + "::" + heldout
            train = [row for (row) in (rows) if (row[axis] != heldout)]
            test = [row for (row) in (rows) if (row[axis] == heldout)]
            inner = [
                {
                    "inner_fold_id": inner_axis + "::" + group,
                    "heldout_group": group,
                    "train_row_ids": [
                        row["row_id"] for (row) in (train) if (row[inner_axis] != group)
                    ],
                    "test_row_ids": [
                        row["row_id"] for (row) in (train) if (row[inner_axis] == group)
                    ],
                }
                for (group) in (inner_groups)
            ]
            require(
                train
                and test
                and all(
                    (
                        item["train_row_ids"] and item["test_row_ids"]
                        for (item) in (inner)
                    )
                ),
                "Incomplete grouped nested split",
            )
            folds.append(
                {
                    "fold_id": fold_id,
                    "scheme": scheme,
                    "axis": axis,
                    "heldout_group": heldout,
                    "inner_axis": inner_axis,
                    "train_row_ids": [row["row_id"] for (row) in (train)],
                    "test_row_ids": [row["row_id"] for (row) in (test)],
                    "inner_folds": inner,
                }
            )
            for (track) in (specs.TRACKS):
                if (tracks is not None and track not in tracks):
                    continue
                for (spec) in (candidates):
                    for (item) in (inner):
                        fields = {
                            "phase": "INNER",
                            "spec_id": spec["spec_id"],
                            "family": spec["family"],
                            "track": track,
                            "fold_id": fold_id,
                            "inner_fold_id": item["inner_fold_id"],
                        }
                        inner_tasks.append(
                            {"task_id": canonical_identity(fields), **fields}
                        )
                for (family) in (families):
                    fields = {
                        "phase": "OUTER",
                        "family": family,
                        "track": track,
                        "fold_id": fold_id,
                    }
                    outer_tasks.append(
                        {"task_id": canonical_identity(fields), **fields}
                    )
    return {
        "schema": "B_UPGRADE_NESTED_DEV_01",
        "expected_rows": len(rows),
        "row_ids": [row["row_id"] for (row) in (rows)],
        "tracks": {
            name: value
            for ((name, value)) in (specs.TRACKS.items())
            if (tracks is None or name in tracks)
        },
        "specs": candidates,
        "folds": folds,
        "inner_tasks": inner_tasks,
        "outer_tasks": outer_tasks,
        "selection_rule": specs.SELECTION_RULE,
        "all_candidates_and_failures_retained": True,
        "protected_access_allowed": False,
    }

def metric(rows, predictions, ids = None):
    import numpy as np
    from sklearn.metrics import average_precision_score

    pmap = {item["row_id"]: item for (item) in (predictions)}
    requested = set(ids) if (ids is not None) else {row["row_id"] for (row) in (rows)}
    population = [row for (row) in (rows) if (row["row_id"] in requested)]
    retained = [
        row
        for (row) in (population)
        if (
            row["target"] is not None and pmap[row["row_id"]]["probability"] is not None
        )
    ]
    retained_ids = {row["row_id"] for (row) in (retained)}
    omitted = [
        row["row_id"] for (row) in (population) if (row["row_id"] not in retained_ids)
    ]
    if (not retained):
        return {
            "total_rows": len(population),
            "available_rows": 0,
            "unavailable_row_ids": omitted,
            "selected_row_ids": [],
            "selected_count": 0,
            "total_events": 0,
            "captured_events": 0,
            "capture": None,
            "brier": None,
            "average_precision": None,
            "calibration": [],
        }
    y = np.asarray([row["target"] for (row) in (retained)], dtype = float)
    p = np.asarray(
        [pmap[row["row_id"]]["probability"] for (row) in (retained)], dtype = float
    )
    require(
        np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all(),
        "Invalid finite probability",
    )
    order = sorted(
        range(len(retained)),
        key = lambda index: (-p[index], retained[index]["row_id"].encode("utf-8")),
    )
    count = math.ceil(0.25 * len(retained))
    event_count, captured = (int(y.sum()), int(y[order[:count]].sum()))
    bins = []
    for (index) in (range(10)):
        mask = (p >= index / 10) & (p < (index + 1) / 10 if (index < 9) else p <= 1)
        bins.append(
            {
                "lower": index / 10,
                "upper": (index + 1) / 10,
                "n": int(mask.sum()),
                "mean_probability": float(p[mask].mean()) if (mask.any()) else None,
                "event_fraction": float(y[mask].mean()) if (mask.any()) else None,
            }
        )
    return {
        "total_rows": len(population),
        "available_rows": len(retained),
        "unavailable_row_ids": omitted,
        "selected_row_ids": [retained[index]["row_id"] for (index) in (order[:count])],
        "selected_count": count,
        "total_events": event_count,
        "captured_events": captured,
        "capture": captured / event_count if (event_count) else None,
        "brier": float(((p - y) ** 2).mean()),
        "average_precision": float(average_precision_score(y, p))
        if (event_count)
        else None,
        "calibration": bins,
    }

def task_metrics(rows, predictions):
    by_cohort = {
        cohort: metric(
            [row for (row) in (rows) if (row["cohort"] == cohort)], predictions
        )
        for (cohort) in (specs.COHORTS)
        if (any((row["cohort"] == cohort for (row) in (rows))))
    }
    captures = [item["capture"] for (item) in (by_cohort.values())]
    briers = [item["brier"] for (item) in (by_cohort.values())]
    return {
        "by_cohort": by_cohort,
        "mean_capture": sum(captures) / len(captures)
        if (captures and all((value is not None for (value) in (captures))))
        else None,
        "mean_brier": sum(briers) / len(briers)
        if (briers and all((value is not None for (value) in (briers))))
        else None,
    }

def unavailable_predictions(rows, reason):
    return [
        {
            "row_id": row["row_id"],
            "status": "UNAVAILABLE",
            "probability": None,
            "reason": reason,
        }
        for (row) in (rows)
    ]

def choose_spec(trials):
    eligible = [
        item
        for (item) in (trials)
        if (
            item["status"] == "COMPLETE"
            and item["mean_capture"] is not None
            and (item["mean_brier"] is not None)
        )
    ]
    return (
        min(
            eligible,
            key = lambda item: (
                -item["mean_capture"],
                item["mean_brier"],
                item["spec_id"],
            ),
        )["spec_id"]
        if (eligible)
        else None
    )

def fit_task(index, task, spec, train_ids, test_ids, fit_cache = None):
    import warnings
    import numpy as np
    from sklearn.exceptions import ConvergenceWarning
    from threadpoolctl import threadpool_limits

    train = [index[identity] for (identity) in (train_ids)]
    test = [index[identity] for (identity) in (test_ids)]
    vectors = {
        row["row_id"]: specs.vector(row, task["track"]) for (row) in (train + test)
    }
    usable = [
        row
        for (row) in (train)
        if (row["target"] is not None and vectors[row["row_id"]] is not None)
    ]
    result = {
        **task,
        "spec_id": spec["spec_id"],
        "train_row_ids": train_ids,
        "test_row_ids": test_ids,
        "train_usable_row_ids": [row["row_id"] for (row) in (usable)],
        "predictions": [],
        "warnings": [],
        "failure": None,
    }
    result["fit_metadata"] = {
        "preprocessing_fit_row_ids": [],
        "tuning": "INNER_FIXED_GRID"
        if (task["phase"] == "INNER")
        else "INNER_ONLY_NESTED_SELECTION_NO_OUTER_LABEL_TUNING",
        "feature_names": specs.TRACKS[task["track"]]
        + ["method_DawnRank", "method_PersonaDrive"],
        "specification": spec,
        "effective_C": 1 / (len(usable) * spec["params"]["penalty"])
        if (usable and "penalty" in spec["params"])
        else None,
        "training_only_preprocessing": True,
        "threads": 1,
        "random_seed": specs.SEED,
        "random_internal_validation": False,
    }
    cache = {} if (fit_cache is None) else fit_cache
    key = canonical_identity(
        {
            "specification": spec,
            "track": task["track"],
            "training": [
                {
                    "row_id": row["row_id"],
                    "target": row["target"],
                    "vector": vectors[row["row_id"]],
                }
                for (row) in (usable)
            ],
        }
    )
    try:
        if (key not in cache):
            warnings_seen = []
            try:
                require(usable, "No complete training rows")
                x = np.asarray(
                    [vectors[row["row_id"]] for (row) in (usable)], dtype = float
                )
                y = np.asarray([row["target"] for (row) in (usable)], dtype = int)
                with (
                    warnings.catch_warnings(record = True) as caught,
                    threadpool_limits(limits = 1),
                ):
                    warnings.simplefilter("always")
                    model = specs.build_estimator(
                        spec, len(usable), len(specs.TRACKS[task["track"]])
                    )
                    constant = model is None or len(set(y)) == 1
                    if (constant):
                        model = {"constant_probability": float(y.mean())}
                    else:
                        model.fit(x, y)
                    warnings_seen = [
                        {
                            "category": item.category.__name__,
                            "message": str(item.message)[:500],
                        }
                        for (item) in (caught)
                    ]
                    if (any(
                        (
                            issubclass(item.category, ConvergenceWarning)
                            for (item) in (caught)
                        )
                    )):
                        raise RuntimeError(
                            "Numerical convergence warning is a preserved failed trial"
                        )
                cache[key] = {
                    "status": "COMPLETE",
                    "model": model,
                    "warnings": warnings_seen,
                    "constant": constant,
                    "prevalence": float(y.mean()),
                    "failure": None,
                }
            except Exception as error:
                cache[key] = {
                    "status": "FAILED",
                    "model": None,
                    "warnings": warnings_seen,
                    "failure": {
                        "type": type(error).__name__,
                        "message": str(error)[:700],
                    },
                }
        fit = cache[key]
        result["warnings"] = fit["warnings"]
        require(
            fit["status"] == "COMPLETE", "Preserved fit failure: " + str(fit["failure"])
        )
        model = fit["model"]
        result["fit_metadata"].update(
            preprocessing_fit_row_ids = result["train_usable_row_ids"],
            constant_training_class_or_prevalence = fit["constant"],
            training_prevalence = fit["prevalence"],
        )
        predictable = [row for (row) in (test) if (vectors[row["row_id"]] is not None)]
        with threadpool_limits(limits = 1):
            if (isinstance(model, dict)):
                probabilities = np.full(len(predictable), model["constant_probability"])
            else:
                probabilities = (
                    model.predict_proba(
                        np.asarray(
                            [vectors[row["row_id"]] for (row) in (predictable)],
                            dtype = float,
                        )
                    )[:, list(model.classes_).index(1)]
                    if (predictable)
                    else np.asarray([])
                )
        require(
            np.isfinite(probabilities).all()
            and ((probabilities >= 0) & (probabilities <= 1)).all(),
            "Invalid fitted probability",
        )
        pmap = {
            row["row_id"]: float(value)
            for ((row, value)) in (zip(predictable, probabilities))
        }
        result["predictions"] = [
            {
                "row_id": row["row_id"],
                "status": "AVAILABLE" if (row["row_id"] in pmap) else "UNAVAILABLE",
                "probability": pmap.get(row["row_id"]),
                "reason": None
                if (row["row_id"] in pmap)
                else "MODEL_INPUT_UNAVAILABLE",
            }
            for (row) in (test)
        ]
        result["status"] = "COMPLETE"
    except Exception as error:
        result.update(
            status = "FAILED",
            predictions = unavailable_predictions(test, "FAILED_TRIAL"),
            failure = {"type": type(error).__name__, "message": str(error)[:700]},
        )
    result["metrics"] = task_metrics(test, result["predictions"])
    return result

def fit_selected(
    rows,
    production = True,
    families = ("extra_trees",),
    tracks = ("localisation_only",),
    schemes = ("LOCO",),
):
    manifest = make_manifest(rows, production, families, tracks, schemes)
    require(
        manifest["inner_tasks"] and manifest["outer_tasks"],
        "Empty requested model or validation scope",
    )
    index = {row["row_id"]: row for (row) in (rows)}
    by_spec = {item["spec_id"]: item for (item) in (manifest["specs"])}
    by_fold = {item["fold_id"]: item for (item) in (manifest["folds"])}
    cache, inner = ({}, {})
    for (task) in (manifest["inner_tasks"]):
        fold = next(
            (
                item
                for (item) in (by_fold[task["fold_id"]]["inner_folds"])
                if (item["inner_fold_id"] == task["inner_fold_id"])
            )
        )
        inner[task["task_id"]] = fit_task(
            index,
            task,
            by_spec[task["spec_id"]],
            fold["train_row_ids"],
            fold["test_row_ids"],
            cache,
        )
    outer, selections = ([], {})
    for (task) in (manifest["outer_tasks"]):
        trials = []
        for (spec) in (manifest["specs"]):
            if (spec["family"] != task["family"]):
                continue
            actual = [
                inner[item["task_id"]]
                for (item) in (manifest["inner_tasks"])
                if (
                    item["fold_id"] == task["fold_id"]
                    and item["track"] == task["track"]
                    and (item["spec_id"] == spec["spec_id"])
                )
            ]
            complete = len(actual) == 2 and all(
                (item["status"] == "COMPLETE" for (item) in (actual))
            )
            capture = (
                [item["metrics"]["mean_capture"] for (item) in (actual)]
                if (complete)
                else []
            )
            brier = (
                [item["metrics"]["mean_brier"] for (item) in (actual)]
                if (complete)
                else []
            )
            trials.append(
                {
                    "spec_id": spec["spec_id"],
                    "inner_task_ids": [item["task_id"] for (item) in (actual)],
                    "status": "COMPLETE" if (complete) else "FAILED_OR_UNAVAILABLE",
                    "mean_capture": sum(capture) / 2
                    if (complete and all((value is not None for (value) in (capture))))
                    else None,
                    "mean_brier": sum(brier) / 2
                    if (complete and all((value is not None for (value) in (brier))))
                    else None,
                }
            )
        selection = {
            "selected_spec_id": choose_spec(trials),
            "trials": trials,
            "rule": specs.SELECTION_RULE,
        }
        selections[task["task_id"]] = selection
        fold = by_fold[task["fold_id"]]
        if (selection["selected_spec_id"] is None):
            test = [index[identity] for (identity) in (fold["test_row_ids"])]
            predictions = unavailable_predictions(test, "NO_ELIGIBLE_INNER_CANDIDATE")
            result = {
                **task,
                "status": "UNAVAILABLE",
                "selection": selection,
                "predictions": predictions,
                "metrics": task_metrics(test, predictions),
                "failure": {"type": "NO_ELIGIBLE_INNER_CANDIDATE"},
            }
        else:
            result = fit_task(
                index,
                {**task, "selection": selection},
                by_spec[selection["selected_spec_id"]],
                fold["train_row_ids"],
                fold["test_row_ids"],
                cache,
            )
        outer.append(result)
    return {
        "manifest": manifest,
        "inner": list(inner.values()),
        "selections": selections,
        "outer": outer,
        "exposure": "Exploratory model development; group-held-out predictions do not remove exposure from adaptive method selection",
    }
