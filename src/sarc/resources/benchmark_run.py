import json

from .benchmark_oracle import expected_answers
from .comparators import STRATEGIES
from .contracts import MappingProposal, PanelSnapshot

def scenario_inputs(scenario):
    left = PanelSnapshot("OpenTargets", "26.06", scenario.panel_id, scenario.baseline)
    if (scenario.identity_decidable):
        right_members = tuple(
            ((member, value) for ((_, member, value)) in (scenario.followup))
        )
        mappings = tuple(
            (
                MappingProposal(source, target, "EXPLICIT_ONE_TO_ONE_REPLACEMENT")
                for ((source, target)) in (scenario.mappings)
            )
        )
    else:
        right_members = tuple(
            (
                (f"unresolved_raw_{row_id}", value)
                for ((row_id, _, value)) in (scenario.followup)
            )
        )
        mappings = ()
    right = PanelSnapshot(
        "OpenTargets",
        "controlled",
        scenario.panel_id,
        right_members,
        scenario.identity_decidable,
    )
    return (left, right, mappings)

def assert_endpoint(actual, expected, context):
    for (key) in (("estimable", "value", "support_denominator")):
        if (actual[key] != expected[key]):
            raise AssertionError(f"{context}: {key} differs")
    if (expected["estimable"]):
        for (key) in (("baseline_set", "followup_set")):
            if (key in expected and actual[key] != expected[key]):
                raise AssertionError(f"{context}: {key} differs")

def evaluate_truth(report, truth):
    assert_endpoint(
        report["fixed"]["leader"],
        truth["tasks"]["fixed_leading_set_change"],
        "fixed leader",
    )
    assert_endpoint(
        report["native"]["leader"],
        truth["tasks"]["native_leading_set_change"],
        "native leader",
    )
    assert_endpoint(
        report["membership"],
        truth["tasks"]["persistent_versus_roster_membership_description"],
        "membership",
    )
    checks = 3
    if (truth["identity_decidable"]):
        if (report["support"]["fixed"] != truth["fixed_count"]):
            raise AssertionError("fixed support differs")
        for (top_k, expected) in (truth["top_k"].items()):
            for (view) in (("native", "fixed")):
                assert_endpoint(
                    report[view]["top_k"][top_k], expected[view], f"{view} top {top_k}"
                )
                if (expected[view]["estimable"]):
                    jaccard = report[view]["top_k"][top_k]["jaccard"]
                    if ((jaccard["numerator"], jaccard["denominator"]) != (
                        expected[view]["jaccard_numerator"],
                        expected[view]["jaccard_denominator"],
                    )):
                        raise AssertionError("Jaccard support differs")
                checks += 1
        if (truth["rank"] is not None):
            for (metric, expected) in (truth["rank"].items()):
                actual = report["fixed"]["rank"].get(metric)
                if (expected is None and actual is not None):
                    raise AssertionError("unavailable correlation was filled")
                if (expected is not None and (
                    actual is None or abs(actual - expected) > 1e-12
                )):
                    raise AssertionError(f"independent {metric} differs")
                checks += 1
    return checks

def task_records(instance_id, panel_id, scenario_name, phase, report, truth):
    rows = []
    for (strategy) in (STRATEGIES):
        for (task, expected) in (truth["tasks"].items()):
            actual = report["strategies"][strategy]["answers"][task]
            answered = actual["state"] == "ANSWERED"
            correct = (
                answered
                and expected["estimable"]
                and (actual["answer"] == expected["value"])
            )
            mismatch = (
                actual["state"] != "NOT_REPORTED"
                and actual["support_denominator"] != expected["support_denominator"]
            )
            rows.append(
                {
                    "instance_id": instance_id,
                    "source_panel_id": panel_id,
                    "phase": phase,
                    "scenario": scenario_name,
                    "strategy": strategy,
                    "task": task,
                    "oracle_eligible": expected["estimable"],
                    "answer_state": actual["state"],
                    "correct_answer": correct,
                    "denominator_mismatch": mismatch,
                    "expected_support_denominator": expected["support_denominator"],
                    "reported_support_denominator": actual["support_denominator"],
                    "truth_boolean": expected["value"]
                    if (isinstance(expected["value"], bool))
                    else None,
                    "answer_boolean": actual["answer"]
                    if (isinstance(actual["answer"], bool))
                    else None,
                    "answer_json": json.dumps(
                        actual["answer"],
                        sort_keys = True,
                        separators = (",", ":"),
                        ensure_ascii = False,
                        allow_nan = False,
                    ),
                    "reason": actual["reason"],
                }
            )
    return rows

def compact_report(report):
    result = {
        key: report[key]
        for (key) in (
            (
                "resource",
                "releases",
                "panel_id",
                "support",
                "native",
                "fixed",
                "wall_seconds",
            )
        )
    }
    return result
