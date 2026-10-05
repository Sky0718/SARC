from __future__ import annotations
from collections import Counter
from typing import Any
from typing import Mapping
from typing import Sequence

def build_turnover_decomposition(
    source_rows: Sequence[Mapping[str, Any]], configuration: Mapping[str, Any]
) -> list[dict[str, Any]]:
    from .ranking_1 import parse_boolean

    reason_classes = list(configuration["top_sets"]["turnover_reason_classes"])
    output: list[dict[str, Any]] = []
    for (pair) in ([str(item["label"]) for (item) in (configuration["release_pairs"])]):
        for (definition) in (["TOP_5", "TOP_10", "TOP_20", "TOP_RANKED"]):
            for (population, eligibility_field) in ([
                ("PRIMARY_FIXED_ANCHORED", "primary_fixed_anchored_eligible"),
                ("ALL_NATIVE_SUPPLEMENTARY", "all_native_supplementary_eligible"),
            ]):
                eligible = [
                    row
                    for (row) in (source_rows)
                    if (
                        row["release_pair"] == pair
                        and row["set_definition"] == definition
                        and parse_boolean(row[eligibility_field])
                        and parse_boolean(row["estimable"])
                    )
                ]
                changed = [
                    row for (row) in (eligible) if (parse_boolean(row["changed"]))
                ]
                counts = Counter((str(row["change_reason"]) for (row) in (eligible)))
                for (reason) in (reason_classes):
                    count = counts.get(reason, 0)
                    changed_count = sum(
                        (
                            str(row["change_reason"]) == reason
                            and parse_boolean(row["changed"])
                            for (row) in (eligible)
                        )
                    )
                    output.append(
                        {
                            "release_pair": pair,
                            "set_definition": definition,
                            "population": population,
                            "reason_class": reason,
                            "estimable_panel_count": len(eligible),
                            "changed_panel_count": len(changed),
                            "reason_panel_count": count,
                            "reason_panel_fraction": None
                            if (not eligible)
                            else count / len(eligible),
                            "reason_changed_panel_fraction": None
                            if (not changed)
                            else changed_count / len(changed),
                            "descriptive_not_causal": True,
                        }
                    )
    return output

def assessment_row(
    scope: str,
    criterion: str,
    numerator: int | None,
    denominator: int | None,
    value: float | bool | None,
    comparison: str,
    threshold: float | bool | None,
    passed: bool | None,
    validation_role: str,
    notes: str,
) -> dict[str, Any]:
    return {
        "release_pair": scope,
        "population": "PRIMARY_FIXED_DISEASE_PANELS",
        "criterion": criterion,
        "numerator": numerator,
        "denominator": denominator,
        "value": value,
        "comparison": comparison,
        "threshold": threshold,
        "passed": passed,
        "validation_role": validation_role,
        "notes": notes,
    }

def proportion_count(values: Sequence[bool]) -> tuple[int, int, float | None]:
    numerator = sum(values)
    denominator = len(values)
    return (
        numerator,
        denominator,
        None if (denominator == 0) else numerator / denominator,
    )

def build_sensitivity_signal_assessment(
    disease_rows: Sequence[Mapping[str, Any]],
    therapeutic_rows: Sequence[Mapping[str, Any]],
    configuration: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, bool]]:
    from .ranking_1 import (
        evaluate_signal_assessment,
        evaluate_therapeutic_area_heterogeneity,
        finite_float,
        parse_boolean,
    )

    output: list[dict[str, Any]] = []
    validation_by_pair: dict[str, bool] = {}
    signal_by_pair: dict[str, bool] = {}
    threshold_changed = finite_float(
        configuration["signal_assessment"]["changed_tie_aware_top_ranked_set_fraction"]
    )
    threshold_overlap = finite_float(
        configuration["signal_assessment"]["native_top_10_jaccard_below"]
    )
    threshold_overlap_fraction = finite_float(
        configuration["signal_assessment"]["native_top_10_low_overlap_fraction"]
    )
    threshold_displacement = finite_float(
        configuration["signal_assessment"][
            "fixed_panel_p90_normalised_displacement_at_least"
        ]
    )
    threshold_displacement_fraction = finite_float(
        configuration["signal_assessment"]["fixed_panel_high_displacement_fraction"]
    )
    pair_labels = [str(pair["label"]) for (pair) in (configuration["release_pairs"])]
    scopes = [*pair_labels, "ALL_ADJACENT_PAIRS"]
    for (scope) in (scopes):
        rows = [
            row
            for (row) in (disease_rows)
            if (scope == "ALL_ADJACENT_PAIRS" or row["release_pair"] == scope)
        ]
        fixed_changed_raw = [
            parse_boolean(row["fixed_top_ranked_set_changed"]) for (row) in (rows)
        ]
        fixed_changed_mapping_excluded = [
            parse_boolean(row["fixed_top_ranked_set_changed"])
            and str(row["fixed_top_ranked_change_reason"]) != "MAPPING_AFFECTED"
            for (row) in (rows)
        ]
        native_rows = [
            row
            for (row) in (rows)
            if (row.get("native_top_10_jaccard") not in {None, ""})
        ]
        native_low = [
            finite_float(row["native_top_10_jaccard"]) < threshold_overlap
            for (row) in (native_rows)
        ]
        fixed_low = [
            finite_float(row["fixed_top_10_jaccard"]) < threshold_overlap
            for (row) in (rows)
        ]
        high_displacement = [
            finite_float(row["p90_normalised_rank_displacement"])
            >= threshold_displacement
            for (row) in (rows)
        ]
        native_changed = [
            row
            for (row) in (native_rows)
            if (parse_boolean(row["native_top_ranked_set_changed"]))
        ]
        mapping_context = [
            str(row["native_top_ranked_change_reason"]) in {"MAPPING_AFFECTED", "MIXED"}
            and int(row["native_top_ranked_mapping_affected_count"]) > 0
            for (row) in (native_changed)
        ]
        entrant_exit_context = [
            int(row["native_top_ranked_entrant_count"]) > 0
            or int(row["native_top_ranked_exit_count"]) > 0
            for (row) in (native_changed)
        ]
        (changed_n, changed_d, changed_value) = proportion_count(
            fixed_changed_mapping_excluded
        )
        (raw_changed_n, raw_changed_d, raw_changed_value) = proportion_count(
            fixed_changed_raw
        )
        (native_n, native_d, native_value) = proportion_count(native_low)
        (fixed_low_n, fixed_low_d, fixed_low_value) = proportion_count(fixed_low)
        (displacement_n, displacement_d, displacement_value) = proportion_count(
            high_displacement
        )
        (mapping_n, mapping_d, mapping_value) = proportion_count(mapping_context)
        (turnover_n, turnover_d, turnover_value) = proportion_count(
            entrant_exit_context
        )
        if (scope == "ALL_ADJACENT_PAIRS"):
            pair_heterogeneity = [
                validation_by_pair.get(pair, False) for (pair) in (pair_labels)
            ]
            heterogeneity = {
                "well_populated_area_count": sum(pair_heterogeneity),
                "fixed_changed_fraction_range": None,
                "native_low_overlap_fraction_range": None,
                "fixed_qualifying_changed_area_count": None,
                "native_qualifying_changed_area_count": None,
                "qualifying_changed_area_count": None,
                "fixed_metric_passed": None,
                "native_metric_passed": None,
                "passing_metric": "ANY_ADJACENT_PAIR"
                if (any(pair_heterogeneity))
                else "NONE",
                "passed": any(pair_heterogeneity),
            }
        else:
            heterogeneity = evaluate_therapeutic_area_heterogeneity(
                [row for (row) in (therapeutic_rows) if (row["release_pair"] == scope)],
                configuration,
            )
            validation_by_pair[scope] = bool(heterogeneity["passed"])
        if (scope == "ALL_ADJACENT_PAIRS"):
            heterogeneity_numerator = sum(
                (validation_by_pair.get(pair, False) for (pair) in (pair_labels))
            )
            heterogeneity_denominator = len(pair_labels)
        elif (heterogeneity["passing_metric"] == "FIXED_TOP_RANKED_SET_CHANGED_FRACTION"):
            heterogeneity_numerator = heterogeneity[
                "fixed_qualifying_changed_area_count"
            ]
            heterogeneity_denominator = heterogeneity["fixed_well_populated_area_count"]
        elif (heterogeneity["passing_metric"] == "NATIVE_TOP_10_LOW_OVERLAP_FRACTION"):
            heterogeneity_numerator = heterogeneity[
                "native_qualifying_changed_area_count"
            ]
            heterogeneity_denominator = heterogeneity[
                "native_well_populated_area_count"
            ]
        else:
            heterogeneity_numerator = None
            heterogeneity_denominator = None
        prerequisite = evaluate_signal_assessment(
            changed_value,
            native_value,
            displacement_value,
            bool(heterogeneity["passed"]),
            configuration,
        )
        if (scope != "ALL_ADJACENT_PAIRS"):
            signal_by_pair[scope] = bool(
                prerequisite["at_least_one_prespecified_signal"]
            )
            combined_signal_numerator = sum(
                (
                    prerequisite[name]
                    for (name) in (
                        [
                            "changed_tie_aware_top_ranked_set_fraction",
                            "native_top_10_low_overlap_fraction",
                            "fixed_panel_high_displacement_fraction",
                            "therapeutic_area_heterogeneity",
                        ]
                    )
                )
            )
            combined_signal_denominator = 4
            combined_signal_value = bool(
                prerequisite["at_least_one_prespecified_signal"]
            )
        else:
            combined_signal_numerator = sum(
                (signal_by_pair.get(pair, False) for (pair) in (pair_labels))
            )
            combined_signal_denominator = len(pair_labels)
            combined_signal_value = any(
                (signal_by_pair.get(pair, False) for (pair) in (pair_labels))
            )
        output.extend(
            [
                assessment_row(
                    scope,
                    "changed_tie_aware_fixed_top_ranked_set_fraction_mapping_only_excluded",
                    changed_n,
                    changed_d,
                    changed_value,
                    ">=",
                    threshold_changed,
                    prerequisite["changed_tie_aware_top_ranked_set_fraction"],
                    "PRESPECIFIED_SIGNAL",
                    "Mapping-only changes are excluded from the numerator and retained in the denominator",
                ),
                assessment_row(
                    scope,
                    "changed_tie_aware_fixed_top_ranked_set_fraction_raw",
                    raw_changed_n,
                    raw_changed_d,
                    raw_changed_value,
                    "DESCRIPTIVE",
                    None,
                    None,
                    "REQUIRED_REPORT",
                    "Tie-aware fixed-support top-ranked target sets",
                ),
                assessment_row(
                    scope,
                    "native_top_10_jaccard_below_0_90_fraction",
                    native_n,
                    native_d,
                    native_value,
                    ">=",
                    threshold_overlap_fraction,
                    prerequisite["native_top_10_low_overlap_fraction"],
                    "PRESPECIFIED_SIGNAL",
                    "Jaccard cutoff is strictly below 0.90",
                ),
                assessment_row(
                    scope,
                    "fixed_top_10_jaccard_below_0_90_fraction",
                    fixed_low_n,
                    fixed_low_d,
                    fixed_low_value,
                    "DESCRIPTIVE",
                    threshold_overlap,
                    None,
                    "REQUIRED_REPORT",
                    "Fixed-support counterpart reported outside the signal decision",
                ),
                assessment_row(
                    scope,
                    "fixed_panel_p90_displacement_at_least_0_10_fraction",
                    displacement_n,
                    displacement_d,
                    displacement_value,
                    ">=",
                    threshold_displacement_fraction,
                    prerequisite["fixed_panel_high_displacement_fraction"],
                    "PRESPECIFIED_SIGNAL",
                    "Panel P90 displacement cutoff is at least 0.10",
                ),
                assessment_row(
                    scope,
                    "therapeutic_area_heterogeneity",
                    heterogeneity_numerator,
                    heterogeneity_denominator,
                    bool(heterogeneity["passed"]),
                    "FROZEN_RULE",
                    True,
                    prerequisite["therapeutic_area_heterogeneity"],
                    "PRESPECIFIED_SIGNAL",
                    f"Uses one metric at a time for both the range and changed-panel-count rule; passing_metric={heterogeneity['passing_metric']};fixed_well_populated_areas={heterogeneity.get('fixed_well_populated_area_count')};native_well_populated_areas={heterogeneity.get('native_well_populated_area_count')};fixed_qualifying_areas={heterogeneity['fixed_qualifying_changed_area_count']};native_qualifying_areas={heterogeneity['native_qualifying_changed_area_count']};fixed_fraction_range={heterogeneity['fixed_changed_fraction_range']};native_fraction_range={heterogeneity['native_low_overlap_fraction_range']}",
                ),
                assessment_row(
                    scope,
                    "changed_native_top_ranked_panels_with_mapping_context_fraction",
                    mapping_n,
                    mapping_d,
                    mapping_value,
                    "DESCRIPTIVE",
                    None,
                    None,
                    "REQUIRED_REPORT",
                    "Contextual association only; not causal attribution",
                ),
                assessment_row(
                    scope,
                    "changed_native_top_ranked_panels_with_entrant_or_exit_fraction",
                    turnover_n,
                    turnover_d,
                    turnover_value,
                    "DESCRIPTIVE",
                    None,
                    None,
                    "REQUIRED_REPORT",
                    "Release-native membership context",
                ),
                assessment_row(
                    scope,
                    "at_least_one_prespecified_signal",
                    combined_signal_numerator,
                    combined_signal_denominator,
                    combined_signal_value,
                    "AT_LEAST_ONE",
                    True,
                    combined_signal_value,
                    "SIGNAL_DECISION",
                    "Data-release validation condition, not a biological threshold",
                ),
            ]
        )
    decision = {
        "at_least_one_prespecified_signal": any(
            (
                parse_boolean(row["passed"])
                for (row) in (output)
                if (
                    row["release_pair"] in pair_labels
                    and row["criterion"] == "at_least_one_prespecified_signal"
                )
            )
        )
    }
    return (output, decision)
