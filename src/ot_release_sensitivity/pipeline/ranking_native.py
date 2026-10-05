from __future__ import annotations
from typing import Any
from typing import Mapping
from typing import Sequence
from .ranking_types import DuckDbDataset, RankingSensitivityBuildError

def native_registry_rows(connection: Any) -> dict[tuple[str, str], dict[str, Any]]:
    from .ranking_1 import query_dicts

    rows = query_dicts(
        connection,
        "SELECT * FROM release_native_panel_registry ORDER BY release, canonical_disease_id",
    )
    result: dict[tuple[str, str], dict[str, Any]] = {}
    for (row) in (rows):
        key = (str(row["release"]), str(row["canonical_disease_id"]))
        if (key in result):
            raise RankingSensitivityBuildError(
                f"Duplicate release-native registry key: {key}"
            )
        result[key] = row
    return result

def native_pair_query(pair: Mapping[str, Any]) -> str:
    from .ranking_1 import sql_literal

    old_release = sql_literal(pair["old_release"])
    new_release = sql_literal(pair["new_release"])
    release_pair = sql_literal(pair["label"])
    return f"\n        WITH paired_diseases AS (\n            SELECT o.canonical_disease_id\n            FROM release_native_panel_registry o\n            JOIN release_native_panel_registry n USING (canonical_disease_id)\n            WHERE o.release = {old_release} AND n.release = {new_release}\n        ),\n        old_members AS (\n            SELECT m.canonical_disease_id, m.canonical_target_id, m.score, m.evidence_count, m.stable_entity_member\n            FROM release_native_panel_members m\n            JOIN paired_diseases p USING (canonical_disease_id)\n            WHERE m.release = {old_release}\n        ),\n        new_members AS (\n            SELECT m.canonical_disease_id, m.canonical_target_id, m.score, m.evidence_count, m.stable_entity_member\n            FROM release_native_panel_members m\n            JOIN paired_diseases p USING (canonical_disease_id)\n            WHERE m.release = {new_release}\n        ),\n        joined AS (\n            SELECT\n                COALESCE(o.canonical_disease_id, n.canonical_disease_id) AS canonical_disease_id,\n                COALESCE(o.canonical_target_id, n.canonical_target_id) AS canonical_target_id,\n                o.score AS old_score,\n                n.score AS new_score,\n                o.evidence_count AS old_evidence_count,\n                n.evidence_count AS new_evidence_count,\n                o.score IS NOT NULL AS old_present,\n                n.score IS NOT NULL AS new_present,\n                COALESCE(o.stable_entity_member, false) AS old_stable_entity_member,\n                COALESCE(n.stable_entity_member, false) AS new_stable_entity_member,\n                COALESCE(o.stable_entity_member, false) AND COALESCE(n.stable_entity_member, false) AS stable_entity_member\n            FROM old_members o\n            FULL OUTER JOIN new_members n USING (canonical_disease_id, canonical_target_id)\n        ),\n        context AS (\n            SELECT\n                canonical_disease_id,\n                canonical_target_id,\n                bool_or(COALESCE(raw_identifier_mapping_affected, false)) AS mapping_affected,\n                max(COALESCE(record_category, '')) AS record_category,\n                max(COALESCE(transition_class, '')) AS transition_class\n            FROM entrant_exit_table\n            WHERE release_pair = {release_pair}\n            GROUP BY canonical_disease_id, canonical_target_id\n        )\n        SELECT\n            {release_pair} AS release_pair,\n            {old_release} AS old_release,\n            {new_release} AS new_release,\n            j.*,\n            COALESCE(c.mapping_affected, false) AS mapping_affected,\n            COALESCE(c.record_category, '') AS record_category,\n            COALESCE(c.transition_class, '') AS transition_class\n        FROM joined j\n        LEFT JOIN context c USING (canonical_disease_id, canonical_target_id)\n        ORDER BY j.canonical_disease_id, j.canonical_target_id\n    "

def calculate_native_panel(
    rows: Sequence[Mapping[str, Any]], cutoffs: Sequence[int], minimum_score: float
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, Any]]:
    from .ranking_fixed import context_class_for_target
    from .ranking_fixed import make_top_set_member_rows
    from .ranking_1 import (
        classify_turnover_reason,
        finite_float,
        jaccard_overlap,
        joined_identifiers,
        parse_boolean,
        tie_aware_top_ranked_set,
        tie_aware_top_set,
    )

    old_scores = {
        str(row["canonical_target_id"]): finite_float(row["old_score"])
        for (row) in (rows)
        if (parse_boolean(row["old_present"]))
    }
    new_scores = {
        str(row["canonical_target_id"]): finite_float(row["new_score"])
        for (row) in (rows)
        if (parse_boolean(row["new_present"]))
    }
    if (any(
        (
            value <= minimum_score
            for (value) in ([*old_scores.values(), *new_scores.values()])
        )
    )):
        raise RankingSensitivityBuildError(
            "Release-native panel contains a non-positive score"
        )
    context = {
        str(row["canonical_target_id"]): {
            "mapping_affected": parse_boolean(row["mapping_affected"]),
            "record_category": row["record_category"],
            "transition_class": row["transition_class"],
            "old_present": parse_boolean(row["old_present"]),
            "new_present": parse_boolean(row["new_present"]),
        }
        for (row) in (rows)
    }
    metrics: dict[str, Any] = {
        "old_native_target_count": len(old_scores),
        "new_native_target_count": len(new_scores),
        "old_native_distinct_score_count": len(set(old_scores.values())),
        "new_native_distinct_score_count": len(set(new_scores.values())),
        "persistent_member_count": len(set(old_scores) & set(new_scores)),
        "entrant_member_count": len(set(new_scores) - set(old_scores)),
        "exit_member_count": len(set(old_scores) - set(new_scores)),
        "mapping_affected_member_count": sum(
            (
                parse_boolean(value.get("mapping_affected", False))
                for (value) in (context.values())
            )
        ),
    }
    member_rows: list[dict[str, Any]] = []
    set_records: dict[str, dict[str, Any]] = {}
    for (cutoff) in (cutoffs):
        old_set = tie_aware_top_set(old_scores, cutoff, minimum_score)
        new_set = tie_aware_top_set(new_scores, cutoff, minimum_score)
        reason = classify_turnover_reason(old_set, new_set, context)
        prefix = f"native_top_{cutoff}"
        estimable = old_set is not None and new_set is not None
        changed_targets = None if (not estimable) else old_set ^ new_set
        metrics.update(
            {
                f"{prefix}_jaccard": jaccard_overlap(old_set, new_set),
                f"old_{prefix}_effective_size": None
                if (old_set is None)
                else len(old_set),
                f"new_{prefix}_effective_size": None
                if (new_set is None)
                else len(new_set),
                f"{prefix}_changed": None if (not estimable) else old_set != new_set,
                f"{prefix}_change_reason": reason,
                f"{prefix}_entrant_count": None
                if (changed_targets is None)
                else sum(
                    (
                        context_class_for_target(context[target]) == "ENTRANT"
                        for (target) in (changed_targets)
                    )
                ),
                f"{prefix}_exit_count": None
                if (changed_targets is None)
                else sum(
                    (
                        context_class_for_target(context[target]) == "EXIT"
                        for (target) in (changed_targets)
                    )
                ),
                f"{prefix}_mapping_affected_count": None
                if (changed_targets is None)
                else sum(
                    (
                        context_class_for_target(context[target]) == "MAPPING_AFFECTED"
                        for (target) in (changed_targets)
                    )
                ),
            }
        )
        member_rows.extend(
            make_top_set_member_rows(
                "RELEASE_NATIVE",
                str(rows[0]["release_pair"]),
                str(rows[0]["canonical_disease_id"]),
                f"TOP_{cutoff}",
                cutoff,
                old_set,
                new_set,
                context,
            )
        )
        set_records[f"TOP_{cutoff}"] = {
            "old_set": old_set,
            "new_set": new_set,
            "reason": reason,
        }
    old_top = tie_aware_top_ranked_set(old_scores, minimum_score)
    new_top = tie_aware_top_ranked_set(new_scores, minimum_score)
    top_reason = classify_turnover_reason(old_top, new_top, context)
    top_estimable = old_top is not None and new_top is not None
    changed_top_targets = None if (not top_estimable) else old_top ^ new_top
    metrics.update(
        {
            "old_native_top_ranked_target_ids": joined_identifiers(old_top),
            "new_native_top_ranked_target_ids": joined_identifiers(new_top),
            "old_native_top_ranked_set_size": None
            if (old_top is None)
            else len(old_top),
            "new_native_top_ranked_set_size": None
            if (new_top is None)
            else len(new_top),
            "native_top_ranked_jaccard": jaccard_overlap(old_top, new_top),
            "native_top_ranked_set_changed": None
            if (not top_estimable)
            else old_top != new_top,
            "native_top_ranked_sets_disjoint": None
            if (not top_estimable)
            else not bool(old_top & new_top),
            "native_top_ranked_change_reason": top_reason,
            "native_top_ranked_entrant_count": None
            if (changed_top_targets is None)
            else sum(
                (
                    context_class_for_target(context[target]) == "ENTRANT"
                    for (target) in (changed_top_targets)
                )
            ),
            "native_top_ranked_exit_count": None
            if (changed_top_targets is None)
            else sum(
                (
                    context_class_for_target(context[target]) == "EXIT"
                    for (target) in (changed_top_targets)
                )
            ),
            "native_top_ranked_mapping_affected_count": None
            if (changed_top_targets is None)
            else sum(
                (
                    context_class_for_target(context[target]) == "MAPPING_AFFECTED"
                    for (target) in (changed_top_targets)
                )
            ),
        }
    )
    member_rows.extend(
        make_top_set_member_rows(
            "RELEASE_NATIVE",
            str(rows[0]["release_pair"]),
            str(rows[0]["canonical_disease_id"]),
            "TOP_RANKED",
            None,
            old_top,
            new_top,
            context,
        )
    )
    set_records["TOP_RANKED"] = {
        "old_set": old_top,
        "new_set": new_top,
        "reason": top_reason,
    }
    top_ranked_row = {
        "analysis_type": "RELEASE_NATIVE",
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
    return (metrics, member_rows, top_ranked_row | {"_set_records": set_records})

def build_native_outputs(
    connection: Any,
    configuration: Mapping[str, Any],
    primary_fixed: Mapping[tuple[str, str], Mapping[str, Any]],
    top_member_dataset: DuckDbDataset,
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[tuple[str, str], dict[str, Any]],
    list[dict[str, Any]],
    dict[str, int],
]:
    from .ranking_export import output_specifications
    from .ranking_1 import (
        append_duckdb_dataset_rows,
        finite_float,
        parse_boolean,
        query_grouped_dicts,
    )

    registry = native_registry_rows(connection)
    cutoffs = [int(value) for (value) in (configuration["top_sets"]["cutoffs"])]
    minimum_score = finite_float(configuration["top_sets"]["minimum_positive_score"])
    metrics_rows: list[dict[str, Any]] = []
    top_ranked_rows: list[dict[str, Any]] = []
    top_member_fields = output_specifications()["top_set_members"]["fields"]
    primary_native: dict[tuple[str, str], dict[str, Any]] = {}
    turnover_source_rows: list[dict[str, Any]] = []
    diagnostics = {
        "native_missing_registry_count": 0,
        "native_member_without_registry_count": int(
            connection.execute(
                "SELECT COUNT(*) FROM (SELECT DISTINCT release, canonical_disease_id FROM release_native_panel_members EXCEPT SELECT release, canonical_disease_id FROM release_native_panel_registry)"
            ).fetchone()[0]
        ),
        "native_registry_without_members_count": int(
            connection.execute(
                "SELECT COUNT(*) FROM (SELECT release, canonical_disease_id FROM release_native_panel_registry EXCEPT SELECT DISTINCT release, canonical_disease_id FROM release_native_panel_members)"
            ).fetchone()[0]
        ),
        "native_duplicate_member_key_count": int(
            connection.execute(
                "SELECT COUNT(*) FROM (SELECT release, canonical_disease_id, canonical_target_id FROM release_native_panel_members GROUP BY ALL HAVING COUNT(*) <> 1)"
            ).fetchone()[0]
        ),
        "native_invalid_evidence_count_count": int(
            connection.execute(
                "SELECT COUNT(*) FROM release_native_panel_members WHERE evidence_count IS NULL OR TRY_CAST(evidence_count AS BIGINT) IS NULL OR TRY_CAST(evidence_count AS BIGINT) < 1"
            ).fetchone()[0]
        ),
        "native_release_contract_mismatch_count": int(
            connection.execute(
                "SELECT COUNT(*) FROM (SELECT DISTINCT release FROM release_native_panel_members UNION SELECT DISTINCT release FROM release_native_panel_registry) WHERE release NOT IN ('25.12', '26.03', '26.06')"
            ).fetchone()[0]
        ),
        "native_roster_count_mismatch_count": 0,
        "native_registry_eligibility_mismatch_count": 0,
        "native_stable_registry_mismatch_count": 0,
        "native_reason_unresolved_changed_count": 0,
        "native_reason_partition_mismatch_count": 0,
    }
    for (pair) in (configuration["release_pairs"]):
        query = native_pair_query(pair)
        for (key_tuple, rows) in (query_grouped_dicts(
            connection, query, ["canonical_disease_id"]
        )):
            disease = str(key_tuple[0])
            old_registry = registry.get((str(pair["old_release"]), disease))
            new_registry = registry.get((str(pair["new_release"]), disease))
            if (old_registry is None or new_registry is None):
                diagnostics["native_missing_registry_count"] += 1
                continue
            (metrics, members, top_row) = calculate_native_panel(
                rows, cutoffs, minimum_score
            )
            if (
                int(old_registry["native_target_count"])
                != metrics["old_native_target_count"]
            ):
                diagnostics["native_roster_count_mismatch_count"] += 1
            if (
                int(new_registry["native_target_count"])
                != metrics["new_native_target_count"]
            ):
                diagnostics["native_roster_count_mismatch_count"] += 1
            old_stable_count = sum(
                (parse_boolean(row["old_stable_entity_member"]) for (row) in (rows))
            )
            new_stable_count = sum(
                (parse_boolean(row["new_stable_entity_member"]) for (row) in (rows))
            )
            if (
                int(old_registry["stable_native_target_count"]) != old_stable_count
                or int(new_registry["stable_native_target_count"]) != new_stable_count
            ):
                diagnostics["native_stable_registry_mismatch_count"] += 1
            registry_checks = [
                int(old_registry["distinct_score_count"])
                == metrics["old_native_distinct_score_count"],
                int(new_registry["distinct_score_count"])
                == metrics["new_native_distinct_score_count"],
                parse_boolean(old_registry["release_native_panel_eligible"])
                == (metrics["old_native_target_count"] > 0),
                parse_boolean(new_registry["release_native_panel_eligible"])
                == (metrics["new_native_target_count"] > 0),
            ]
            for (cutoff) in (cutoffs):
                registry_checks.extend(
                    [
                        parse_boolean(old_registry[f"top_{cutoff}_estimable"])
                        == (metrics["old_native_target_count"] >= cutoff),
                        parse_boolean(new_registry[f"top_{cutoff}_estimable"])
                        == (metrics["new_native_target_count"] >= cutoff),
                    ]
                )
            if (not all(registry_checks)):
                diagnostics["native_registry_eligibility_mismatch_count"] += 1
            native_eligible = (
                metrics["old_native_target_count"] > 0
                and metrics["new_native_target_count"] > 0
            )
            primary_eligible = (
                str(pair["label"]),
                disease,
            ) in primary_fixed and native_eligible
            label = old_registry["canonical_label"] or new_registry["canonical_label"]
            base = {
                "release_pair": pair["label"],
                "old_release": pair["old_release"],
                "new_release": pair["new_release"],
                "canonical_disease_id": disease,
                "canonical_label": label,
                "canonical_therapeutic_area_ids": old_registry[
                    "canonical_therapeutic_area_ids"
                ],
                "primary_fixed_anchored_eligible": primary_eligible,
                "all_native_supplementary_eligible": native_eligible,
            }
            combined = base | metrics
            metrics_rows.append(combined)
            if (native_eligible):
                append_duckdb_dataset_rows(
                    connection, top_member_dataset, top_member_fields, members
                )
                top_ranked_rows.append(
                    {
                        key: value
                        for ((key, value)) in (top_row.items())
                        if (key != "_set_records")
                    }
                )
                for (definition, record) in (top_row["_set_records"].items()):
                    estimable = (
                        record["old_set"] is not None and record["new_set"] is not None
                    )
                    changed = estimable and record["old_set"] != record["new_set"]
                    turnover_source_rows.append(
                        {
                            "release_pair": pair["label"],
                            "canonical_disease_id": disease,
                            "set_definition": definition,
                            "change_reason": record["reason"],
                            "estimable": estimable,
                            "changed": changed,
                            "primary_fixed_anchored_eligible": primary_eligible,
                            "all_native_supplementary_eligible": native_eligible,
                        }
                    )
                    if (estimable and changed and (record["reason"] == "UNRESOLVED")):
                        diagnostics["native_reason_unresolved_changed_count"] += 1
                    if (estimable and changed != (record["reason"] != "UNCHANGED")):
                        diagnostics["native_reason_partition_mismatch_count"] += 1
            if (primary_eligible):
                primary_native[str(pair["label"]), disease] = combined
    primary_fixed_keys = set(primary_fixed)
    primary_native_keys = set(primary_native)
    reported_primary_keys = {
        (str(row["release_pair"]), str(row["canonical_disease_id"]))
        for (row) in (metrics_rows)
        if (parse_boolean(row["primary_fixed_anchored_eligible"]))
    }
    diagnostics["native_primary_fixed_key_set_mismatch_count"] = len(
        primary_fixed_keys ^ primary_native_keys
    )
    diagnostics["native_primary_fixed_count_mismatch_count"] = abs(
        len(primary_fixed_keys) - len(primary_native_keys)
    )
    diagnostics["native_primary_eligibility_key_set_mismatch_count"] = len(
        primary_fixed_keys ^ reported_primary_keys
    )
    diagnostics["native_primary_top_set_estimability_mismatch_count"] = sum(
        (
            metrics is None
            or any(
                (
                    int(metrics[f"old_native_target_count"]) < cutoff
                    or int(metrics[f"new_native_target_count"]) < cutoff
                    or metrics.get(f"native_top_{cutoff}_jaccard") is None
                    for (cutoff) in (cutoffs)
                )
            )
            or metrics.get("native_top_ranked_jaccard") is None
            or (metrics.get("native_top_ranked_set_changed") is None)
            for (metrics) in (
                (primary_native.get(key) for (key) in (sorted(primary_fixed_keys)))
            )
        )
    )
    return (
        metrics_rows,
        top_ranked_rows,
        primary_native,
        turnover_source_rows,
        diagnostics,
    )
