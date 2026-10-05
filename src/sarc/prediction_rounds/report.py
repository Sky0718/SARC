import math
from ..evaluation.prediction import require, metric, choose_spec
from .adapters import adapter, ASSISTED_TRACKS

def merge_predictions(models, states, identity, result):
    target = models.setdefault(identity, {})
    for (prediction) in (result["predictions"]):
        require(prediction["row_id"] not in target, "Repeated out-of-fold row")
        target[prediction["row_id"]] = prediction
    states[identity] = (
        "COMPLETE"
        if (result["status"] == "COMPLETE"
        and states.get(identity, "COMPLETE") == "COMPLETE")
        else "FAILED_OR_UNAVAILABLE"
    )

def models_from_result(rows, result, round_id):
    specs = adapter(round_id)
    manifest = result["manifest"]
    require(
        manifest["specs"] == specs.specifications()
        and manifest["tracks"] == specs.TRACKS,
        "Complete frozen round specifications and tracks required",
    )
    folds = {fold["fold_id"]: fold for (fold) in (manifest["folds"])}
    require(len(folds) == 4, "Both complete grouped validation schemes required")
    require(
        len(result["inner"]) == len(manifest["inner_tasks"])
        and len(result["outer"]) == len(manifest["outer_tasks"]),
        "All finite task records including failures required",
    )
    inner = {record["task_id"]: record for (record) in (result["inner"])}
    outer = {record["task_id"]: record for (record) in (result["outer"])}
    require(
        set(inner) == {task["task_id"] for (task) in (manifest["inner_tasks"])}
        and set(outer) == {task["task_id"] for (task) in (manifest["outer_tasks"])},
        "Missing or duplicate task identity",
    )
    models, states, ledger = {}, {}, []
    for (task) in (manifest["outer_tasks"]):
        record = outer[task["task_id"]]
        require(
            all(record[key] == task[key] for (key) in (("track", "family", "fold_id"))),
            "Outer task identity differs",
        )
        require(
            [item["row_id"] for (item) in (record["predictions"])]
            == folds[task["fold_id"]]["test_row_ids"],
            "Full ordered outer population required",
        )
        identity = (folds[task["fold_id"]]["scheme"], task["track"], task["family"])
        merge_predictions(models, states, identity, record)
    if (round_id not in ("R1", "R2")):
        for (track) in (manifest["tracks"]):
            for (fold_id, fold) in (folds.items()):
                trials = []
                for (spec) in (manifest["specs"]):
                    tasks = [
                        task
                        for (task) in (manifest["inner_tasks"])
                        if (task["track"] == track
                        and task["fold_id"] == fold_id
                        and task["spec_id"] == spec["spec_id"])
                    ]
                    require(
                        len(tasks) == 2, "Each candidate requires both inner contexts"
                    )
                    values = [inner[task["task_id"]] for (task) in (tasks)]
                    complete = all(
                        value["status"] == "COMPLETE"
                        and value["metrics"]["mean_capture"] is not None
                        and value["metrics"]["mean_brier"] is not None
                        for (value) in (values)
                    )
                    trials.append(
                        {
                            "spec_id": spec["spec_id"],
                            "family": spec["family"],
                            "inner_task_ids": [task["task_id"] for (task) in (tasks)],
                            "status": "COMPLETE"
                            if (complete)
                            else "FAILED_OR_UNAVAILABLE",
                            "mean_capture": math.fsum(
                                value["metrics"]["mean_capture"] for (value) in (values)
                            )
                            / 2
                            if (complete)
                            else None,
                            "mean_brier": math.fsum(
                                value["metrics"]["mean_brier"] for (value) in (values)
                            )
                            / 2
                            if (complete)
                            else None,
                        }
                    )
                winner_id = choose_spec(trials)
                winner = next(
                    (trial for (trial) in (trials) if (trial["spec_id"] == winner_id)), None
                )
                task = (
                    next(
                        (
                            task
                            for (task) in (manifest["outer_tasks"])
                            if (task["track"] == track
                            and task["fold_id"] == fold_id
                            and task["family"] == winner["family"])
                        ),
                        None,
                    )
                    if (winner)
                    else None
                )
                record = outer[task["task_id"]] if (task) else None
                require(
                    record is None or record["spec_id"] == winner_id,
                    "Across-family inner winner differs from selected-family outer fit",
                )
                if (record is None):
                    record = {
                        "status": "UNAVAILABLE",
                        "predictions": [
                            {
                                "row_id": key,
                                "status": "UNAVAILABLE",
                                "probability": None,
                                "reason": "NO_ELIGIBLE_META_INNER_CANDIDATE",
                            }
                            for (key) in (fold["test_row_ids"])
                        ],
                    }
                identity = (fold["scheme"], track, "inner_selected_any")
                merge_predictions(models, states, identity, record)
                ledger.append(
                    {
                        "scheme": fold["scheme"],
                        "track": track,
                        "fold_id": fold_id,
                        "trials": trials,
                        "selected_spec_id": winner_id,
                        "selected_outer_task_id": task["task_id"] if (task) else None,
                        "status": record["status"],
                        "rule": specs.SELECTION_RULE,
                        "outer_labels_used_to_select_family": False,
                    }
                )
    expected = {row["row_id"] for (row) in (rows)}
    require(
        all(set(predictions) == expected for (predictions) in (models.values())),
        "Incomplete full out-of-fold denominator",
    )
    return models, states, ledger

def common_input_ids(rows, specs):
    tracks = (
        ["old_inputs"]
        if (specs.round_id == "R7")
        else ["strict_localisation_mechanism"]
        if (specs.round_id in ("R5", "R6"))
        else list(specs.TRACKS)
    )
    return {
        row["row_id"]
        for (row) in (rows)
        if (all(specs.vector(row, track) is not None for (track) in (tracks)))
    }

def subgroup_metric(rows, predictions, support):
    full = metric(rows, list(predictions.values()), support)
    persona = [row for (row) in (rows) if (row["method"] == "PersonaDrive")]
    selected = set(full["selected_row_ids"])
    applicable = [row for (row) in (persona) if (row["row_id"] in support)]
    retained = [row for (row) in (applicable) if (row["row_id"] in selected)]
    events = sum(row["target"] for (row) in (applicable))
    captured = sum(row["target"] for (row) in (retained))
    conditional = metric(persona, list(predictions.values()), support)
    return {
        "method": "PersonaDrive",
        "full_method_rows": len(persona),
        "available_rows": len(applicable),
        "global_all_method_selected_count": full["selected_count"],
        "selected_persona_count_within_global_budget": len(retained),
        "total_events": events,
        "captured_events_within_global_budget": captured,
        "capture": captured / events if (events) else None,
        "brier": conditional["brier"],
        "average_precision": conditional["average_precision"],
        "calibration": conditional["calibration"],
        "global_selected_row_ids": full["selected_row_ids"],
        "subgroup_budget_reallocated": False,
    }

def summaries(rows, round_id, models, states):
    specs = adapter(round_id)
    common = common_input_ids(rows, specs)
    result = []
    for (identity, predictions) in (models.items()):
        groups = {}
        for (cohort) in (specs.COHORTS):
            unit = [row for (row) in (rows) if (row["cohort"] == cohort)]
            native = [
                row["row_id"]
                for (row) in (unit)
                if (row["target"] is not None
                and predictions[row["row_id"]]["probability"] is not None)
            ]
            paired = [key for (key) in (native) if (key in common)]
            groups[cohort] = {
                "full_cohort_rows": len(unit),
                "native_support_row_ids": native,
                "common_input_row_ids": paired,
                "native_metrics": metric(unit, list(predictions.values()), native)
                if (states[identity] == "COMPLETE")
                else None,
                "common_input_metrics": metric(unit, list(predictions.values()), paired)
                if (states[identity] == "COMPLETE")
                else None,
            }
            if (round_id in ("R5", "R6")):
                groups[cohort]["persona_drive_subgroup"] = (
                    subgroup_metric(unit, predictions, set(native))
                    if (states[identity] == "COMPLETE")
                    else None
                )
        result.append(
            {
                "round": round_id,
                "scheme": identity[0],
                "track": identity[1],
                "family": identity[2],
                "status": states[identity],
                "by_cohort": groups,
            }
        )
    return result

def applicable_baseline(current, candidate, baseline_round, baseline):
    same_information = (
        candidate[1] in ASSISTED_TRACKS or baseline[1] not in ASSISTED_TRACKS
    )
    if (baseline_round != current):
        return same_information
    tracks = adapter(current).TRACKS
    if (current == "R7"):
        return baseline[2] in ("anchor", "same_feature_extra_trees") and set(
            tracks[baseline[1]]
        ) <= set(tracks[candidate[1]])
    if (current in ("R3", "R4", "R5") and baseline[2] != "inner_selected_any"):
        return False
    if (current == "R6" and baseline[2] != "inner_selected_any"):
        return (
            set(tracks[baseline[1]]) < set(tracks[candidate[1]])
            or baseline[1] == candidate[1]
            and baseline[2] == "xgb_binary_logistic"
        )
    return same_information and set(tracks[baseline[1]]) < set(tracks[candidate[1]])

def comparisons(rows, round_id, current, pools):
    specs = adapter(round_id)
    candidates, states = current
    panels, gates, cached, supports = [], [], {}, {}

    def measure(label, identity, predictions, cohort, support):
        key = (label, identity, cohort, tuple(support))
        if (key not in cached):
            cached[key] = metric(
                [row for (row) in (rows) if (row["cohort"] == cohort)],
                list(predictions.values()),
                support,
            )
        return cached[key]

    for (identity, predictions) in (candidates.items()):
        relevant = []
        for (label, models, baseline_states) in (pools):
            for (baseline_id, baseline) in (models.items()):
                if (
                    baseline_id[0] != identity[0]
                    or label == round_id
                    and baseline_id == identity
                ):
                    continue
                applicable = applicable_baseline(round_id, identity, label, baseline_id)
                panel = {
                    "candidate": list(identity),
                    "baseline": [label, *baseline_id],
                    "applicable_gate_baseline": applicable,
                    "by_cohort": {},
                    "status": "COMPLETE"
                    if (states[identity] == baseline_states[baseline_id] == "COMPLETE")
                    else "FAILED_OR_UNAVAILABLE",
                }
                full_support = True
                for (cohort) in (specs.COHORTS):
                    unit = [row for (row) in (rows) if (row["cohort"] == cohort)]
                    support = [
                        row["row_id"]
                        for (row) in (unit)
                        if (row["target"] is not None
                        and predictions[row["row_id"]]["probability"] is not None
                        and baseline[row["row_id"]]["probability"] is not None)
                    ]
                    support_key = (cohort, tuple(support))
                    if (support_key not in supports):
                        retained = set(support)
                        supports[support_key] = {
                            "id": len(supports),
                            "cohort": cohort,
                            "row_ids": support,
                            "omitted_row_ids": [
                                row["row_id"]
                                for (row) in (unit)
                                if (row["row_id"] not in retained)
                            ],
                        }
                    full_support &= len(support) == sum(
                        row["target"] is not None for (row) in (unit)
                    )
                    left = (
                        measure(round_id, identity, predictions, cohort, support)
                        if (panel["status"] == "COMPLETE")
                        else None
                    )
                    right = (
                        measure(label, baseline_id, baseline, cohort, support)
                        if (panel["status"] == "COMPLETE")
                        else None
                    )
                    valid = (
                        left is not None
                        and left["capture"] is not None
                        and right["capture"] is not None
                    )
                    panel["by_cohort"][cohort] = {
                        "full_rows": len(unit),
                        "support_count": len(support),
                        "support_id": supports[support_key]["id"],
                        "candidate": left,
                        "baseline": right,
                        "capture_gain": left["capture"] - right["capture"]
                        if (valid)
                        else None,
                        "brier_change": left["brier"] - right["brier"]
                        if (valid)
                        else None,
                        "average_precision_change": left["average_precision"]
                        - right["average_precision"]
                        if (valid)
                        else None,
                    }
                    if (round_id in ("R5", "R6")):
                        panel["by_cohort"][cohort]["persona_drive_subgroup"] = (
                            {
                                "candidate": subgroup_metric(
                                    unit, predictions, set(support)
                                ),
                                "baseline": subgroup_metric(
                                    unit, baseline, set(support)
                                ),
                                "used_for_primary_gate": False,
                            }
                            if (valid)
                            else None
                        )
                valid = all(
                    value["capture_gain"] is not None
                    for (value) in (panel["by_cohort"].values())
                )
                for (name) in (("capture_gain", "brier_change")):
                    panel["mean_" + name] = (
                        math.fsum(value[name] for (value) in (panel["by_cohort"].values()))
                        / 2
                        if (valid)
                        else None
                    )
                panel["mean_baseline_capture"] = (
                    math.fsum(
                        value["baseline"]["capture"]
                        for (value) in (panel["by_cohort"].values())
                    )
                    / 2
                    if (valid)
                    else None
                )
                panel["mean_baseline_brier"] = (
                    math.fsum(
                        value["baseline"]["brier"]
                        for (value) in (panel["by_cohort"].values())
                    )
                    / 2
                    if (valid)
                    else None
                )
                panel["full_paired_available_support"] = bool(full_support)
                panel["panel_index"] = len(panels)
                panels.append(panel)
                if (applicable):
                    relevant.append(panel)
        if (
            identity[0] == "LOCO"
            and round_id not in ("R1", "R2")
            and identity[2] == "inner_selected_any"
        ):
            eligible = [
                panel
                for (panel) in (relevant)
                if (panel["mean_baseline_capture"] is not None)
            ]
            strongest = (
                min(
                    eligible,
                    key = lambda panel: (
                        -panel["mean_baseline_capture"],
                        panel["mean_baseline_brier"],
                        tuple(panel["baseline"]),
                    ),
                )
                if (eligible)
                else None
            )
            conditions = {
                "candidate_complete": states[identity] == "COMPLETE",
                "all_applicable_baselines_complete": bool(relevant)
                and all(panel["status"] == "COMPLETE" for (panel) in (relevant)),
                "mean_capture_gain_at_least_0_10": strongest is not None
                and strongest["mean_capture_gain"] >= 0.10,
                "positive_capture_gain_each_cohort": strongest is not None
                and all(
                    value["capture_gain"] > 0
                    for (value) in (strongest["by_cohort"].values())
                ),
                "mean_brier_worsening_at_most_0_005": strongest is not None
                and strongest["mean_brier_change"] <= 0.005,
                "brier_worsening_at_most_0_005_each_cohort": strongest is not None
                and all(
                    value["brier_change"] <= 0.005
                    for (value) in (strongest["by_cohort"].values())
                ),
            }
            if (round_id == "R7"):
                conditions["full_paired_support"] = bool(relevant) and all(
                    panel["full_paired_available_support"] for (panel) in (relevant)
                )
            gates.append(
                {
                    "candidate": list(identity),
                    "strongest_applicable_baseline_panel_index": strongest[
                        "panel_index"
                    ]
                    if (strongest)
                    else None,
                    "conditions": conditions,
                    "passes_development_research_target": all(conditions.values()),
                    "independent_confirmation": False,
                    "adaptive_exploratory_only": True,
                }
            )
    return {
        "panels": panels,
        "gates": gates,
        "support_ledgers": list(supports.values()),
        "adaptive_development_not_confirmation": True,
    }

def summarise(rows, round_id, result, prior):
    expected_prior = ["R" + str(index) for (index) in (range(1, int(round_id[1:])))]
    require(
        list(prior) == expected_prior,
        "All preceding finite rounds, in order, are required for the complete strong-baseline comparison",
    )
    pools, ledgers = [], {}
    for (label, value) in ([*prior.items(), (round_id, result)]):
        models, states, ledger = models_from_result(rows, value, label)
        pools.append((label, models, states))
        ledgers[label] = ledger
    current = pools[-1]
    return {
        "round": round_id,
        "summaries": summaries(rows, round_id, current[1], current[2]),
        "meta_selection_ledger": ledgers[round_id],
        "comparison": comparisons(rows, round_id, current[1:], pools),
        "failed_tasks": [
            record
            for (record) in (result["inner"] + result["outer"])
            if (record["status"] != "COMPLETE")
        ],
        "exposure": "All adaptive round results remain exploratory; no independent confirmation is claimed",
    }
