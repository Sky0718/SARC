from __future__ import annotations
import math
from typing import Any
from typing import Mapping
from typing import Sequence
from .ranking_types import DuckDbDataset, RankingSensitivityBuildError

def context_class_for_target(context: Mapping[str, Any]) -> str:
    from .ranking_1 import parse_boolean

    mapping = parse_boolean(context.get("mapping_affected", False))
    record_category = str(context.get("record_category", "")).upper()
    transition_class = str(context.get("transition_class", "")).upper()
    old_present = parse_boolean(context.get("old_present", False))
    new_present = parse_boolean(context.get("new_present", False))
    if (
        mapping
        or record_category == "MAPPING_AFFECTED"
        or "MAPPING" in transition_class
    ):
        return "MAPPING_AFFECTED"
    if (not old_present and new_present):
        return "ENTRANT"
    if (old_present and (not new_present)):
        return "EXIT"
    if (old_present and new_present):
        return "PERSISTENT_PAIR_RERANKING"
    return "UNRESOLVED"

def make_top_set_member_rows(
    analysis_type: str,
    release_pair: str,
    disease: str,
    definition: str,
    cutoff: int | None,
    old_set: set[str] | None,
    new_set: set[str] | None,
    context: Mapping[str, Mapping[str, Any]],
) -> list[dict[str, Any]]:
    from .ranking_1 import membership_state, parse_boolean

    if (old_set is None or new_set is None):
        return []
    rows: list[dict[str, Any]] = []
    for (target) in (sorted(old_set | new_set)):
        target_context = context.get(target, {})
        rows.append(
            {
                "analysis_type": analysis_type,
                "release_pair": release_pair,
                "canonical_disease_id": disease,
                "set_definition": definition,
                "cutoff": cutoff,
                "canonical_target_id": target,
                "old_member": target in old_set,
                "new_member": target in new_set,
                "membership_state": membership_state(target, old_set, new_set),
                "context_class": context_class_for_target(target_context),
                "mapping_affected": parse_boolean(
                    target_context.get("mapping_affected", False)
                ),
            }
        )
    return rows

def calculate_fixed_panel(
    rows: Sequence[Mapping[str, Any]], cutoffs: Sequence[int], minimum_score: float
) -> tuple[
    dict[str, Any],
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[str, Mapping[str, Any]],
]:
    from .ranking_1 import (
        classify_turnover_reason,
        finite_float,
        jaccard_overlap,
        joined_identifiers,
        kendall_tau_b,
        linear_quantile,
        normalised_rank_displacements,
        parse_boolean,
        spearman_correlation,
        tie_aware_top_ranked_set,
        tie_aware_top_set,
    )

    targets = [str(row["canonical_target_id"]) for (row) in (rows)]
    if (len(targets) != len(set(targets))):
        raise RankingSensitivityBuildError(
            "Fixed-support panel contains duplicate canonical targets"
        )
    old_scores = [finite_float(row["old_score"]) for (row) in (rows)]
    new_scores = [finite_float(row["new_score"]) for (row) in (rows)]
    if (any((value <= minimum_score for (value) in (old_scores + new_scores)))):
        raise RankingSensitivityBuildError(
            "Fixed-support panel contains a non-positive score"
        )
    (old_ranks, new_ranks, displacement) = normalised_rank_displacements(
        old_scores, new_scores
    )
    old_map = dict(zip(targets, old_scores))
    new_map = dict(zip(targets, new_scores))
    context = {
        target: {
            "mapping_affected": parse_boolean(row.get("mapping_affected", False)),
            "record_category": row.get("record_category", ""),
            "transition_class": row.get("transition_class", ""),
            "old_present": True,
            "new_present": True,
        }
        for ((target, row)) in (zip(targets, rows))
    }
    metrics: dict[str, Any] = {
        "fixed_target_count": len(targets),
        "old_distinct_score_count": len(set(old_scores)),
        "new_distinct_score_count": len(set(new_scores)),
        "kendall_tau_b": kendall_tau_b(old_scores, new_scores),
        "spearman_correlation": spearman_correlation(old_scores, new_scores),
        "median_normalised_rank_displacement": linear_quantile(displacement, 0.5),
        "p90_normalised_rank_displacement": linear_quantile(displacement, 0.9),
        "maximum_normalised_rank_displacement": max(displacement),
    }
    member_rows: list[dict[str, Any]] = []
    for (cutoff) in (cutoffs):
        old_set = tie_aware_top_set(old_map, cutoff, minimum_score)
        new_set = tie_aware_top_set(new_map, cutoff, minimum_score)
        prefix = f"fixed_top_{cutoff}"
        metrics[f"{prefix}_jaccard"] = jaccard_overlap(old_set, new_set)
        metrics[f"old_{prefix}_effective_size"] = (
            None if (old_set is None) else len(old_set)
        )
        metrics[f"new_{prefix}_effective_size"] = (
            None if (new_set is None) else len(new_set)
        )
        metrics[f"{prefix}_changed"] = (
            None if (old_set is None or new_set is None) else old_set != new_set
        )
        member_rows.extend(
            make_top_set_member_rows(
                "FIXED_SUPPORT",
                str(rows[0]["release_pair"]),
                str(rows[0]["canonical_disease_id"]),
                f"TOP_{cutoff}",
                cutoff,
                old_set,
                new_set,
                context,
            )
        )
    old_top = tie_aware_top_ranked_set(old_map, minimum_score)
    new_top = tie_aware_top_ranked_set(new_map, minimum_score)
    top_reason = classify_turnover_reason(old_top, new_top, context)
    metrics.update(
        {
            "old_fixed_top_ranked_target_ids": joined_identifiers(old_top),
            "new_fixed_top_ranked_target_ids": joined_identifiers(new_top),
            "old_fixed_top_ranked_set_size": None
            if (old_top is None)
            else len(old_top),
            "new_fixed_top_ranked_set_size": None
            if (new_top is None)
            else len(new_top),
            "fixed_top_ranked_jaccard": jaccard_overlap(old_top, new_top),
            "fixed_top_ranked_set_changed": None
            if (old_top is None or new_top is None)
            else old_top != new_top,
            "fixed_top_ranked_sets_disjoint": None
            if (old_top is None or new_top is None)
            else not bool(old_top & new_top),
            "fixed_top_ranked_change_reason": top_reason,
            "fixed_top_ranked_mapping_only": top_reason == "MAPPING_AFFECTED",
        }
    )
    member_rows.extend(
        make_top_set_member_rows(
            "FIXED_SUPPORT",
            str(rows[0]["release_pair"]),
            str(rows[0]["canonical_disease_id"]),
            "TOP_RANKED",
            None,
            old_top,
            new_top,
            context,
        )
    )
    displacement_rows = [
        {
            "release_pair": row["release_pair"],
            "old_release": row["old_release"],
            "new_release": row["new_release"],
            "canonical_disease_id": row["canonical_disease_id"],
            "canonical_target_id": target,
            "old_score": old_score,
            "new_score": new_score,
            "old_rank": old_rank,
            "new_rank": new_rank,
            "normalised_absolute_rank_displacement": target_displacement,
            "stable_entity_member": parse_boolean(row["stable_entity_member"]),
            "mapping_affected": parse_boolean(row.get("mapping_affected", False)),
        }
        for ((
            row,
            target,
            old_score,
            new_score,
            old_rank,
            new_rank,
            target_displacement,
        )) in (
            zip(
                rows,
                targets,
                old_scores,
                new_scores,
                old_ranks,
                new_ranks,
                displacement,
            )
        )
    ]
    top_ranked_row = {
        "analysis_type": "FIXED_SUPPORT",
        "release_pair": rows[0]["release_pair"],
        "canonical_disease_id": rows[0]["canonical_disease_id"],
        "old_top_ranked_target_ids": joined_identifiers(old_top),
        "new_top_ranked_target_ids": joined_identifiers(new_top),
        "old_set_size": None if (old_top is None) else len(old_top),
        "new_set_size": None if (new_top is None) else len(new_top),
        "jaccard": jaccard_overlap(old_top, new_top),
        "changed": None if (old_top is None or new_top is None) else old_top != new_top,
        "disjoint": None
        if (old_top is None or new_top is None)
        else not bool(old_top & new_top),
        "change_reason": top_reason,
        "descriptive_not_causal": True,
    }
    return (
        metrics,
        displacement_rows,
        member_rows,
        context | {"__top_ranked__": top_ranked_row},
    )

def fixed_registry_rows(
    connection: Any, configuration: Mapping[str, Any]
) -> dict[tuple[str, str], dict[str, Any]]:
    from .ranking_1 import query_dicts, sql_literal

    support = sql_literal(configuration["fixed_support"]["support_name"])
    rows = query_dicts(
        connection,
        f"SELECT * FROM fixed_ranking_panel_registry WHERE support_name = {support} ORDER BY release_pair, canonical_disease_id",
    )
    result: dict[tuple[str, str], dict[str, Any]] = {}
    for (row) in (rows):
        key = (str(row["release_pair"]), str(row["canonical_disease_id"]))
        if (key in result):
            raise RankingSensitivityBuildError(
                f"Duplicate fixed-panel registry key: {key}"
            )
        result[key] = row
    return result

def build_fixed_outputs(
    connection: Any, configuration: Mapping[str, Any], top_member_dataset: DuckDbDataset
) -> tuple[
    list[dict[str, Any]],
    DuckDbDataset,
    list[dict[str, Any]],
    dict[tuple[str, str], dict[str, Any]],
    dict[str, int],
]:
    from .ranking_export import output_specifications
    from .ranking_1 import (
        append_duckdb_dataset_rows,
        create_duckdb_dataset,
        finite_float,
        parse_boolean,
        query_grouped_dicts,
        sql_literal,
    )

    registry = fixed_registry_rows(connection, configuration)
    support = sql_literal(configuration["fixed_support"]["support_name"])
    query = f"\n        WITH context AS (\n            SELECT\n                release_pair,\n                canonical_disease_id,\n                canonical_target_id,\n                bool_or(COALESCE(raw_identifier_mapping_affected, false)) AS mapping_affected,\n                max(COALESCE(record_category, '')) AS record_category,\n                max(COALESCE(transition_class, '')) AS transition_class\n            FROM entrant_exit_table\n            GROUP BY release_pair, canonical_disease_id, canonical_target_id\n        )\n        SELECT\n            m.*,\n            COALESCE(c.mapping_affected, false) AS mapping_affected,\n            COALESCE(c.record_category, '') AS record_category,\n            COALESCE(c.transition_class, '') AS transition_class\n        FROM fixed_ranking_panel_members m\n        LEFT JOIN context c USING (release_pair, canonical_disease_id, canonical_target_id)\n        WHERE m.support_name = {support}\n        ORDER BY m.release_pair, m.canonical_disease_id, m.canonical_target_id\n    "
    cutoffs = [int(value) for (value) in (configuration["top_sets"]["cutoffs"])]
    minimum_score = finite_float(configuration["top_sets"]["minimum_positive_score"])
    primary = int(configuration["fixed_support"]["primary_minimum_targets"])
    thresholds = sorted(
        {
            primary,
            *[
                int(value)
                for (value) in (
                    configuration["fixed_support"]["sensitivity_minimum_targets"]
                )
            ],
        }
    )
    metric_rows: list[dict[str, Any]] = []
    displacement_specification = output_specifications()[
        "fixed_target_rank_displacement"
    ]
    displacement_dataset = create_duckdb_dataset(
        connection,
        "ranking_sensitivity_fixed_target_rank_displacement_data",
        displacement_specification["fields"],
        displacement_specification["types"],
        int(configuration["resources"]["bulk_insert_batch_rows"]),
    )
    top_ranked_rows: list[dict[str, Any]] = []
    top_member_fields = output_specifications()["top_set_members"]["fields"]
    primary_metrics: dict[tuple[str, str], dict[str, Any]] = {}
    diagnostics = {
        "fixed_missing_registry_count": 0,
        "fixed_registry_without_members_count": 0,
        "fixed_roster_count_mismatch_count": 0,
        "fixed_registry_eligibility_mismatch_count": 0,
        "fixed_stable_registry_mismatch_count": 0,
        "fixed_invalid_evidence_count_count": 0,
        "fixed_release_pair_contract_mismatch_count": 0,
        "fixed_nonfinite_metric_count": 0,
        "fixed_displacement_out_of_bounds_count": 0,
    }
    seen_registry_keys: set[tuple[str, str]] = set()
    pair_contract = {
        str(pair["label"]): (str(pair["old_release"]), str(pair["new_release"]))
        for (pair) in (configuration["release_pairs"])
    }
    for (key, rows) in (query_grouped_dicts(
        connection, query, ["release_pair", "canonical_disease_id"]
    )):
        registry_key = (str(key[0]), str(key[1]))
        registry_row = registry.get(registry_key)
        if (registry_row is None):
            diagnostics["fixed_missing_registry_count"] += 1
            continue
        seen_registry_keys.add(registry_key)
        expected_releases = pair_contract.get(registry_key[0])
        if (expected_releases is None or any(
            (
                (str(row["old_release"]), str(row["new_release"])) != expected_releases
                for (row) in (rows)
            )
        )):
            diagnostics["fixed_release_pair_contract_mismatch_count"] += 1
            continue
        target_ids = [str(row["canonical_target_id"]) for (row) in (rows)]
        if (len(target_ids) != len(set(target_ids))):
            raise RankingSensitivityBuildError(
                "Fixed-support panel contains duplicate canonical targets"
            )
        old_scores = [finite_float(row["old_score"]) for (row) in (rows)]
        new_scores = [finite_float(row["new_score"]) for (row) in (rows)]
        if (any((value <= minimum_score for (value) in (old_scores + new_scores)))):
            raise RankingSensitivityBuildError(
                "Fixed-support panel contains a non-positive score"
            )
        target_count = len(rows)
        old_distinct = len(set(old_scores))
        new_distinct = len(set(new_scores))
        if (int(registry_row["fixed_target_count"]) != target_count):
            diagnostics["fixed_roster_count_mismatch_count"] += 1
        minimum_distinct = int(
            configuration["fixed_support"]["minimum_distinct_scores_per_release"]
        )
        eligible_by_threshold = {
            threshold: target_count >= threshold
            and old_distinct >= minimum_distinct
            and (new_distinct >= minimum_distinct)
            for (threshold) in (thresholds)
        }
        expected_exclusion = (
            "FEWER_THAN_MINIMUM_FIXED_TARGETS"
            if (target_count < primary)
            else "DEGENERATE_OLD_SCORE_DISTRIBUTION"
            if (old_distinct < minimum_distinct)
            else "DEGENERATE_NEW_SCORE_DISTRIBUTION"
            if (new_distinct < minimum_distinct)
            else ""
        )
        eligibility_mismatch = (
            parse_boolean(registry_row["eligible_min_20"]) != eligible_by_threshold[20]
            or parse_boolean(registry_row["eligible_primary_30"])
            != eligible_by_threshold[30]
            or parse_boolean(registry_row["eligible_min_50"])
            != eligible_by_threshold[50]
            or (
                parse_boolean(registry_row["primary_panel_eligible"])
                != eligible_by_threshold[primary]
            )
            or (str(registry_row["exclusion_reason"] or "") != expected_exclusion)
            or (int(registry_row["old_distinct_score_count"]) != old_distinct)
            or (int(registry_row["new_distinct_score_count"]) != new_distinct)
        )
        if (eligibility_mismatch):
            diagnostics["fixed_registry_eligibility_mismatch_count"] += 1
        stable_rows = [
            row for (row) in (rows) if (parse_boolean(row["stable_entity_member"]))
        ]
        stable_count = len(stable_rows)
        stable_old_distinct = len(
            {finite_float(row["old_score"]) for (row) in (stable_rows)}
        )
        stable_new_distinct = len(
            {finite_float(row["new_score"]) for (row) in (stable_rows)}
        )
        expected_stable = (
            stable_count >= primary
            and stable_old_distinct >= minimum_distinct
            and (stable_new_distinct >= minimum_distinct)
        )
        if (
            int(registry_row["stable_fixed_target_count"]) != stable_count
            or int(registry_row["stable_old_distinct_score_count"])
            != stable_old_distinct
            or int(registry_row["stable_new_distinct_score_count"])
            != stable_new_distinct
            or (parse_boolean(registry_row["stable_panel_eligible"]) != expected_stable)
        ):
            diagnostics["fixed_stable_registry_mismatch_count"] += 1
        diagnostics["fixed_invalid_evidence_count_count"] += sum(
            (
                row["old_evidence_count"] is None
                or row["new_evidence_count"] is None
                or int(row["old_evidence_count"]) < 1
                or (int(row["new_evidence_count"]) < 1)
                for (row) in (rows)
            )
        )
        if (not eligible_by_threshold[20]):
            continue
        (metrics, targets, members, context) = calculate_fixed_panel(
            rows, cutoffs, minimum_score
        )
        expected_primary = eligible_by_threshold[primary]
        if (any(
            (
                value is None or not math.isfinite(float(value))
                for (value) in (
                    [metrics["kendall_tau_b"], metrics["spearman_correlation"]]
                )
                if (expected_primary)
            )
        )):
            diagnostics["fixed_nonfinite_metric_count"] += 1
        diagnostics["fixed_displacement_out_of_bounds_count"] += sum(
            (
                not 0.0
                <= finite_float(row["normalised_absolute_rank_displacement"])
                <= 1.0
                for (row) in (targets)
            )
        )
        base = {
            "release_pair": key[0],
            "old_release": rows[0]["old_release"],
            "new_release": rows[0]["new_release"],
            "canonical_disease_id": key[1],
            "canonical_label": registry_row["canonical_label"],
            "canonical_therapeutic_area_ids": registry_row[
                "canonical_therapeutic_area_ids"
            ],
            "stable_panel_eligible": expected_stable,
        }
        for (threshold) in (thresholds):
            if (eligible_by_threshold[threshold]):
                metric_rows.append(
                    base
                    | metrics
                    | {
                        "panel_population": "PRIMARY_FIXED"
                        if (threshold == primary)
                        else "PANEL_SIZE_SENSITIVITY",
                        "minimum_target_threshold": threshold,
                    }
                )
        if (eligible_by_threshold[20]):
            enriched_targets = [
                row
                | {
                    "eligible_min_20": True,
                    "eligible_primary_30": eligible_by_threshold[30],
                    "eligible_min_50": eligible_by_threshold[50],
                }
                for (row) in (targets)
            ]
            append_duckdb_dataset_rows(
                connection,
                displacement_dataset,
                displacement_specification["fields"],
                enriched_targets,
            )
            append_duckdb_dataset_rows(
                connection, top_member_dataset, top_member_fields, members
            )
            top_ranked_rows.append(dict(context["__top_ranked__"]))
        if (eligible_by_threshold[primary]):
            primary_metrics[str(key[0]), str(key[1])] = base | metrics
    diagnostics["fixed_registry_without_members_count"] = len(
        set(registry) - seen_registry_keys
    )
    expected_primary_keys = {
        key
        for ((key, row)) in (registry.items())
        if (parse_boolean(row["primary_panel_eligible"]))
    }
    diagnostics["fixed_primary_panel_missing_metrics_count"] = len(
        expected_primary_keys - set(primary_metrics)
    )
    return (
        metric_rows,
        displacement_dataset,
        top_ranked_rows,
        primary_metrics,
        diagnostics,
    )
