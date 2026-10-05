from ..evaluation.prediction import (
    require,
    validate_rows,
    canonical_identity,
    choose_spec,
    unavailable_predictions,
    task_metrics,
)
from ..core.io import save_json

def make_manifest(rows, specs):
    production = True
    families = tracks = schemes = None
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
                    (item["train_row_ids"] and item["test_row_ids"] for (item) in (inner))
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
            for (name, value) in (specs.TRACKS.items())
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

def fit_round(rows, specs, fit, output):
    manifest = make_manifest(rows, specs)
    require(
        manifest["inner_tasks"] and manifest["outer_tasks"],
        "Empty requested model or validation scope",
    )
    index = {row["row_id"]: row for (row) in (rows)}
    by_spec = {item["spec_id"]: item for (item) in (manifest["specs"])}
    by_fold = {item["fold_id"]: item for (item) in (manifest["folds"])}
    cache, inner = ({}, {})
    for (task_number, task) in (enumerate(manifest["inner_tasks"])):
        fold = next(
            (
                item
                for (item) in (by_fold[task["fold_id"]]["inner_folds"])
                if (item["inner_fold_id"] == task["inner_fold_id"])
            )
        )
        inner[task["task_id"]] = fit(
            index,
            task,
            by_spec[task["spec_id"]],
            fold["train_row_ids"],
            fold["test_row_ids"],
            cache,
        )
        save_json(
            output / ("inner_" + str(task_number).zfill(4) + ".json"),
            inner[task["task_id"]],
        )
    outer, selections = ([], {})
    for (task_number, task) in (enumerate(manifest["outer_tasks"])):
        trials = []
        for (spec) in (manifest["specs"]):
            if (spec["family"] != task["family"]):
                continue
            actual = [
                inner[item["task_id"]]
                for (item) in (manifest["inner_tasks"])
                if (item["fold_id"] == task["fold_id"]
                and item["track"] == task["track"]
                and (item["spec_id"] == spec["spec_id"]))
            ]
            complete = len(actual) == 2 and all(
                (item["status"] == "COMPLETE" for (item) in (actual))
            )
            capture = (
                [item["metrics"]["mean_capture"] for (item) in (actual)] if (complete) else []
            )
            brier = (
                [item["metrics"]["mean_brier"] for (item) in (actual)] if (complete) else []
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
            result = fit(
                index,
                {**task, "selection": selection},
                by_spec[selection["selected_spec_id"]],
                fold["train_row_ids"],
                fold["test_row_ids"],
                cache,
            )
        outer.append(result)
        save_json(output / ("outer_" + str(task_number).zfill(4) + ".json"), result)
    return {
        "manifest": manifest,
        "inner": list(inner.values()),
        "selections": selections,
        "outer": outer,
        "exposure": "Exploratory model development; group-held-out predictions do not remove exposure from adaptive method selection",
    }
