import json
from collections import Counter

COUNTS = {
    "B": 180,
    "A": 2880,
    "E": 1872,
    "FIXED_POLICY_SELECTION": 3216,
    "C_CANDIDATE": 1296,
    "C_SELECTION": 1368,
    "C_REVIEW": 21456,
}
B_VIEWS = (
    "selected_candidate_recomputation",
    "declared_required_method_work",
    "baseline_results_cached_incremental_work",
)
N_VIEWS = (
    "required_method_or_evaluation_information_work",
    "selected_choice_method_work",
)

def require(value, message):
    if (not value):
        raise ValueError(message)

def same(left, right):
    return json.dumps(left, sort_keys = True, allow_nan = False) == json.dumps(
        right, sort_keys = True, allow_nan = False
    )

def pointer(binding, parts):
    return {"binding": binding, "pointer": parts}

def dims(scenario, fields):
    return {
        name: None if (scenario[name] is None) else scenario[name]["dimensions"]
        for (name) in (fields)
    }

def build(descriptive, nonb, b_packets, refs, checkpoint = lambda: None):
    nr = {row["catalog_id"]: (i, row) for ((i, row)) in (enumerate(nonb["rows"]))}
    ns = {row["scenario_id"]: (i, row) for ((i, row)) in (enumerate(nonb["scenarios"]))}
    require(
        len(nr) == len(nonb["rows"]) and len(ns) == len(nonb["scenarios"]),
        "Duplicate non-B identity",
    )
    br = {}
    for (packet) in (b_packets):
        require(
            packet["document"]["status"]
            == "FIRST_POLICY_SPECIFIC_B_METHOD_WORK_UNION_NOT_END_TO_END_FRONTIER",
            "Accepted B policy schema required",
        )
        for (i, row) in (enumerate(packet["document"]["rows"])):
            identity = (key(packet["accepted_conversion"]), row["row_id"])
            require(identity not in br, "Duplicate B carrier/row identity")
            br[identity] = (packet, i, row)
    (rows, used_n, used_s, used_b) = ([], set(), set(), set())
    for (i, old) in (enumerate(descriptive["rows"])):
        checkpoint()
        row = {
            name: old[name]
            for (name) in (
                (
                    "catalog_id",
                    "family",
                    "original_row_id",
                    "panel_id",
                    "axes",
                    "applicable",
                    "descriptive_point",
                )
            )
        }
        row["accepted_descriptive_row"] = pointer(refs["readout"], ["rows", i])
        row["original_summary_retained_by_reference"] = True
        row["contract_policy_support_budget_and_exposure_retained_by_reference"] = True
        row["end_to_end_cost"] = None
        row["end_to_end_cost_complete"] = False
        row["utility_cost_classification"] = "UNCLASSIFIED"
        row["deployable_saving_claimed"] = False
        row["partial_method_work_display_is_not_makespan_or_end_to_end"] = True
        if (old["family"] == "B"):
            identity = (key(old["source_binding"]), old["original_row_id"])
            require(
                identity in br and identity not in used_b, "Missing or duplicate B join"
            )
            used_b.add(identity)
            (packet, j, cost) = br[identity]
            require(
                same(cost["comparison_key"], old["contract"]["comparison"])
                and cost["information_regime"] == old["contract"]["information_regime"]
                and (cost["family"] == old["policy"]["policy_id"]),
                "B context/support/budget/information policy mismatch",
            )
            effects = cost["unchanged_effects_and_denominators"]
            require(
                all(
                    (
                        effects[name] == value
                        for ((name, value)) in (old["contract"]["counts"].items())
                    )
                ),
                "B denominator changed",
            )
            require(
                old["axes"]
                == [
                    {
                        "name": "event_capture",
                        "direction": "max",
                        "value": effects["event_capture"],
                    },
                    {"name": "brier", "direction": "min", "value": effects["brier"]},
                ],
                "B accepted effect values changed",
            )
            require(
                cost["policy_end_to_end_cost_complete"] is False
                and cost["deployable_saving_claimed"] is False,
                "B partial cost relabelled complete",
            )
            row["accepted_cost_row"] = pointer(packet["binding"], ["rows", j])
            row["accepted_cost_scenario"] = row["accepted_cost_row"]
            row["partial_method_work_dimensions"] = dims(cost, B_VIEWS)
            row["unknown_cost_evidence"] = {
                "ancestor_ledger": pointer(
                    packet["binding"], ["rows", j, "ancestor_category_ledger"]
                ),
                "information_boundary_resolved": cost["information_boundary_resolved"],
                "unresolved_required_source_ids": cost[
                    "unresolved_required_source_ids"
                ],
                "selected_actions_without_cost_ids": cost[
                    "selected_actions_without_cost_ids"
                ],
            }
            row["observed_study_not_attributed_to_policy"] = packet["document"][
                "observed_included_study_work"
            ]
        else:
            require(
                old["catalog_id"] in nr and old["catalog_id"] not in used_n,
                "Missing or duplicate non-B join",
            )
            used_n.add(old["catalog_id"])
            (j, cost) = nr[old["catalog_id"]]
            require(
                cost["family"] == old["family"]
                and cost["original_row_id"] == old["original_row_id"]
                and (key(cost["source_binding"]) == key(old["source_binding"]))
                and (cost["source_pointer"] == old["source_pointer"]),
                "Non-B inherited row identity mismatch",
            )
            require(
                same(cost["comparison_contract"], old["contract"]["comparison"])
                and same(cost["policy_contract"], old["policy"]),
                "Non-B context/support/budget/policy mismatch",
            )
            require(
                cost["scenario_id"] in ns and cost["end_to_end_cost_complete"] is False,
                "Exact non-B scenario required",
            )
            used_s.add(cost["scenario_id"])
            (s, scenario) = ns[cost["scenario_id"]]
            require(
                scenario["end_to_end_cost"] is None
                and scenario["whole_cohort_execution_not_divided_by_stratum_or_budget"]
                is True,
                "Non-B complete or fractional cost invented",
            )
            row["accepted_cost_row"] = pointer(refs["non_b_lineage"], ["rows", j])
            row["accepted_cost_scenario"] = pointer(
                refs["non_b_lineage"], ["scenarios", s]
            )
            row["partial_method_work_dimensions"] = dims(scenario, N_VIEWS)
            row["unknown_cost_evidence"] = {
                "ancestor_artifacts": pointer(
                    refs["non_b_lineage"], ["rows", j, "required_ancestor_artifacts"]
                ),
                "unknown_ancestor_measurements": cost["unknown_ancestor_measurements"],
                "mapping_gaps": scenario["mapping_gaps"],
                "physical_mapping_complete": scenario["physical_mapping_complete"],
            }
        rows.append(row)
    require(
        used_n == set(nr) and used_b == set(br) and (used_s == set(ns)),
        "Orphaned accepted cost row/scenario",
    )
    require(
        len({row["catalog_id"] for (row) in (rows)}) == len(rows),
        "Duplicate descriptive row",
    )
    return {
        "status": "COMPLETE_ROW_COST_REFERENCE_JOIN",
        "rows": rows,
        "counts": dict(Counter((row["family"] for (row) in (rows)))),
        "panels": descriptive["panels"],
        "inherited_descriptive_classification_counts": descriptive[
            "classification_counts"
        ],
        "non_b_scenarios": len(ns),
        "B_cost_rows": len(br),
        "input_bindings": refs,
        "accepted_ancestor_and_process_ledgers_retained_by_reference": True,
        "end_to_end_cost_complete": False,
        "cost_dominance_counts": {"UNCLASSIFIED": len(rows)},
        "limitations": [
            "Required unmeasured categories keep all end-to-end totals unavailable, not zero.",
            "Displayed method-process work is partial, not runtime makespan or deployable saving.",
            "Accepted within-contract descriptive point classifications remain unchanged, not statistical superiority.",
            "Scientific applicability and missing cost categories remain explicit.",
        ],
    }

def key(reference):
    value = reference["source_id"]
    require(
        isinstance(value, str) and bool(value) and value.strip() == value,
        "An explicit scientific source identity is required",
    )
    return value
