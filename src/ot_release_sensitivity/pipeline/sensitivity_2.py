from __future__ import annotations
import csv
import math
from pathlib import Path
from typing import Any
from typing import Mapping
from typing import Sequence
from ot_release_sensitivity.pipeline.ranking import parse_boolean
from .sensitivity_types import CrossLayerSensitivityBuildError

def output_specifications() -> dict[str, dict[str, Any]]:
    score_fields = [
        "release_pair",
        "analysis_layer",
        "support_denominator",
        "eligible_target_pair_count",
        "eligible_disease_panel_count",
        "distinct_target_count",
        "median_signed_change",
        "median_absolute_change",
        "q90_absolute_change",
        "q95_absolute_change",
        "maximum_absolute_change",
        "exact_unchanged_fraction",
        "absolute_change_at_least_0_01_fraction",
        "absolute_change_at_least_0_05_fraction",
        "absolute_change_at_least_0_10_fraction",
    ]
    ranking_fields = [
        "analysis_layer",
        "release_pair",
        "old_release",
        "new_release",
        "canonical_disease_id",
        "canonical_label",
        "canonical_therapeutic_area_ids",
        "fixed_target_count",
        "old_distinct_score_count",
        "new_distinct_score_count",
        "old_positive_score_count",
        "new_positive_score_count",
        "kendall_tau_b",
        "spearman_correlation",
        "median_normalised_rank_displacement",
        "p90_normalised_rank_displacement",
        "maximum_normalised_rank_displacement",
        "stable_entity_member_fraction",
    ]
    for (cutoff) in ([5, 10, 20]):
        ranking_fields.extend(
            [
                f"top_{cutoff}_jaccard",
                f"old_top_{cutoff}_effective_size",
                f"new_top_{cutoff}_effective_size",
                f"top_{cutoff}_changed",
            ]
        )
    ranking_fields.extend(
        [
            "old_top_ranked_target_ids",
            "new_top_ranked_target_ids",
            "old_top_ranked_set_size",
            "new_top_ranked_set_size",
            "top_ranked_jaccard",
            "top_ranked_set_changed",
            "top_ranked_sets_disjoint",
            "support_denominator",
            "overall_primary_panel_available",
            "overall_top_ranked_set_changed",
            "changed_top_target_overlap_jaccard_with_overall",
            "corroborates_overall_top_change",
            "descriptive_not_causal",
        ]
    )
    native_fields = [
        "analysis_layer",
        "release_pair",
        "old_release",
        "new_release",
        "canonical_disease_id",
        "canonical_label",
        "canonical_therapeutic_area_ids",
        "old_native_target_count",
        "new_native_target_count",
        "persistent_member_count",
        "entrant_member_count",
        "exit_member_count",
    ]
    for (cutoff) in ([5, 10, 20]):
        native_fields.extend(
            [
                f"native_top_{cutoff}_jaccard",
                f"old_native_top_{cutoff}_effective_size",
                f"new_native_top_{cutoff}_effective_size",
                f"native_top_{cutoff}_changed",
                f"native_top_{cutoff}_change_reason",
            ]
        )
    native_fields.extend(
        [
            "old_native_top_ranked_target_ids",
            "new_native_top_ranked_target_ids",
            "native_top_ranked_jaccard",
            "native_top_ranked_set_changed",
            "native_top_ranked_change_reason",
            "support_denominator",
            "descriptive_not_causal",
        ]
    )
    return {
        "stable_entity_pair_revisions": {
            "fields": [
                "release_pair",
                "old_release",
                "new_release",
                "canonical_disease_id",
                "canonical_label",
                "canonical_therapeutic_area_ids",
                "canonical_target_id",
                "old_score",
                "new_score",
                "signed_change",
                "absolute_change",
                "old_rank",
                "new_rank",
                "normalised_absolute_rank_displacement",
                "fixed_target_count",
                "stable_entity_eligible",
                "support_denominator",
            ],
            "order": ["release_pair", "canonical_disease_id", "canonical_target_id"],
            "parquet": True,
            "table": "stable_entity_pair_revisions",
        },
        "stable_entity_score_summary": {
            "fields": score_fields,
            "order": ["release_pair"],
            "parquet": False,
        },
        "stable_entity_ranking_metrics": {
            "fields": ranking_fields,
            "order": ["release_pair", "canonical_disease_id"],
            "parquet": False,
        },
        "stable_entity_native_metrics": {
            "fields": native_fields,
            "order": ["release_pair", "canonical_disease_id"],
            "parquet": False,
        },
        "evidence_domain_pair_revisions": {
            "fields": [
                "endpoint_level",
                "endpoint_id",
                "release_pair",
                "old_release",
                "new_release",
                "canonical_disease_id",
                "canonical_label",
                "canonical_therapeutic_area_ids",
                "canonical_target_id",
                "old_score",
                "new_score",
                "signed_change",
                "absolute_change",
                "old_rank",
                "new_rank",
                "normalised_absolute_rank_displacement",
                "fixed_target_count",
                "stable_entity_eligible",
                "old_source_object_relative_path",
                "new_source_object_relative_path",
                "support_denominator",
            ],
            "order": ["release_pair", "canonical_disease_id", "canonical_target_id"],
            "parquet": True,
            "table": "evidence_domain_pair_revisions",
        },
        "evidence_domain_score_summary": {
            "fields": score_fields,
            "order": ["release_pair"],
            "parquet": False,
        },
        "evidence_domain_ranking_metrics": {
            "fields": ranking_fields,
            "order": ["release_pair", "canonical_disease_id"],
            "parquet": False,
        },
        "diagnostic_target_context": {
            "fields": [
                "release_pair",
                "old_release",
                "new_release",
                "canonical_disease_id",
                "canonical_target_id",
                "old_top_ranked_member",
                "new_top_ranked_member",
                "primary_top_ranked_membership",
                "diagnostic_family",
                "old_endpoint_id",
                "new_endpoint_id",
                "comparison_status",
                "old_score",
                "new_score",
                "old_evidence_count",
                "new_evidence_count",
                "diagnostic_presence",
                "signed_change",
                "absolute_change",
                "large_diagnostic_change",
                "old_source_object_relative_path",
                "new_source_object_relative_path",
                "diagnostic_only_not_causal",
            ],
            "order": [
                "release_pair",
                "canonical_disease_id",
                "canonical_target_id",
                "diagnostic_family",
            ],
            "parquet": True,
            "table": "diagnostic_target_context",
        },
        "literature_clinical_diagnostics": {
            "fields": [
                "release_pair",
                "old_release",
                "new_release",
                "canonical_disease_id",
                "diagnostic_family",
                "comparison_status",
                "changed_top_target_count",
                "target_count_with_diagnostic_support",
                "target_count_with_both_release_support",
                "target_count_with_presence_change",
                "target_count_with_large_comparable_change",
                "large_comparable_change_fraction",
                "evidence_group",
                "any_large_comparable_change",
                "diagnostic_only_not_causal",
            ],
            "order": ["release_pair", "canonical_disease_id", "diagnostic_family"],
            "parquet": False,
        },
        "release_note_alignment": {
            "fields": [
                "release",
                "release_pair_context",
                "category",
                "documented_change",
                "analytical_relevance",
                "source_evidence_status",
                "alignment_class",
                "observed_evidence",
                "mechanism_attribution",
                "source_url",
                "source_snapshot",
                "causality_boundary",
                "context_only_not_causal",
            ],
            "order": ["release", "category"],
            "parquet": False,
        },
        "sensitivity_comparison": {
            "fields": [
                "release_pair",
                "analysis_layer",
                "support_denominator",
                "denominator_definition",
                "eligible_disease_panel_count",
                "eligible_target_pair_count",
                "median_absolute_score_change",
                "q90_absolute_score_change",
                "median_kendall_tau_b",
                "top_10_change_numerator",
                "top_10_change_denominator",
                "top_10_change_fraction",
                "top_ranked_change_numerator",
                "top_ranked_change_denominator",
                "top_ranked_change_fraction",
                "mapping_affected_numerator",
                "mapping_affected_denominator",
                "mapping_affected_fraction",
                "mapping_affected_definition",
                "entrant_driven_numerator",
                "entrant_driven_denominator",
                "entrant_driven_fraction",
                "entrant_driven_definition",
                "therapeutic_area_coverage",
                "distinct_support_denominator",
            ],
            "order": ["release_pair", "analysis_layer"],
            "parquet": False,
        },
        "persistence_registry": {
            "fields": [
                "release_pair",
                "canonical_disease_id",
                "canonical_label",
                "canonical_therapeutic_area_ids",
                "stable_therapeutic_area_ids",
                "overall_primary_eligible",
                "stable_entity_eligible",
                "evidence_domain_eligible",
                "overall_top_ranked_set_changed",
                "stable_entity_top_ranked_set_changed",
                "evidence_domain_top_ranked_set_changed",
                "overall_kendall_tau_b",
                "stable_entity_kendall_tau_b",
                "evidence_domain_kendall_tau_b",
                "evidence_domain_corroborates_overall",
                "mapping_context_present",
                "primary_fixed_target_count",
                "stable_fixed_target_count",
                "entity_stability_control_removed_fraction",
                "diagnostic_datasource_context_available",
                "persistence_class",
                "descriptive_not_causal",
            ],
            "order": ["release_pair", "canonical_disease_id"],
            "parquet": False,
        },
    }

def serialise_csv_cell(value: Any) -> Any:
    if (value is None):
        return ""
    if (isinstance(value, bool)):
        return "true" if (value) else "false"
    if (isinstance(value, float)):
        if (not math.isfinite(value)):
            raise CrossLayerSensitivityBuildError(
                "Cannot serialise a non-finite Cross-layer sensitivity value"
            )
        return repr(value)
    return value

def deterministic_sort_key(
    row: Mapping[str, Any], fields: Sequence[str]
) -> tuple[str, ...]:
    return tuple(
        (
            "" if (row.get(field) is None) else str(row.get(field))
            for (field) in (fields)
        )
    )

def write_csv_dataset(
    destination: Path,
    rows: Sequence[Mapping[str, Any]],
    fields: Sequence[str],
    order: Sequence[str],
) -> None:
    destination.parent.mkdir(parents = True, exist_ok = True)
    ordered = sorted(rows, key = lambda row: deterministic_sort_key(row, order))
    with destination.open("w", encoding = "utf-8", newline = "") as handle:
        writer = csv.DictWriter(handle, fieldnames = list(fields), lineterminator = "\n")
        writer.writeheader()
        for (row) in (ordered):
            writer.writerow(
                {field: serialise_csv_cell(row.get(field)) for (field) in (fields)}
            )

def set_export_threads(connection: Any, configuration: Mapping[str, Any]) -> None:
    export_threads = int(configuration["resources"]["export_threads"])
    if (export_threads != 1):
        raise CrossLayerSensitivityBuildError(
            "Ordered Cross-layer sensitivity public exports require one DuckDB thread"
        )
    connection.execute(f"SET threads = {export_threads}")
    observed = int(
        connection.execute("SELECT current_setting('threads')").fetchone()[0]
    )
    if (observed != export_threads):
        raise CrossLayerSensitivityBuildError(
            "DuckDB did not apply the frozen Cross-layer sensitivity export-thread setting"
        )

def export_parquet_table(
    connection: Any,
    table: str,
    destination: Path,
    fields: Sequence[str],
    order: Sequence[str],
    configuration: Mapping[str, Any],
) -> None:
    from .sensitivity_1 import sql_identifier, sql_literal

    destination.parent.mkdir(parents = True, exist_ok = True)
    selection = ", ".join((sql_identifier(field) for (field) in (fields)))
    ordering = ", ".join((sql_identifier(field) for (field) in (order)))
    compression = str(configuration["determinism"]["parquet_compression"])
    row_group_size = int(configuration["determinism"]["parquet_row_group_size"])
    set_export_threads(connection, configuration)
    connection.execute(
        f"COPY (SELECT {selection} FROM {sql_identifier(table)} ORDER BY {ordering}) TO {sql_literal(destination.as_posix())} (FORMAT PARQUET, COMPRESSION {compression}, ROW_GROUP_SIZE {row_group_size})"
    )

def export_stage(
    connection: Any,
    stage: Path,
    datasets: Mapping[str, Sequence[Mapping[str, Any]]],
    configuration: Mapping[str, Any],
) -> dict[str, Path]:
    stage.mkdir(parents = True, exist_ok = True)
    paths: dict[str, Path] = {}
    for (name, specification) in (output_specifications().items()):
        destination = stage / Path(configuration["outputs"][name]).name
        if (parse_boolean(specification["parquet"])):
            export_parquet_table(
                connection,
                str(specification["table"]),
                destination,
                specification["fields"],
                specification["order"],
                configuration,
            )
        else:
            write_csv_dataset(
                destination,
                datasets[name],
                specification["fields"],
                specification["order"],
            )
        paths[name] = destination
    return paths

def build_cross_layer_sensitivity_datasets(
    connection: Any,
    project_root: Path,
    input_paths: Mapping[str, Path],
    configuration: Mapping[str, Any],
    diagnostic_paths: Mapping[str, Sequence[tuple[str, Path]]],
) -> dict[str, list[dict[str, Any]]]:
    from .sensitivity_1 import (
        build_diagnostic_summary_rows,
        build_fixed_metric_rows,
        build_persistence_rows,
        build_release_note_rows,
        build_sensitivity_comparison_rows,
        build_stable_native_rows,
        create_changed_target_table,
        create_diagnostic_tables,
        create_fixed_revision_tables,
        create_stable_area_table,
        query_dicts,
        score_summary_rows,
        sql_literal,
        stable_area_memberships,
    )

    create_stable_area_table(connection, configuration)
    create_fixed_revision_tables(connection, configuration)
    primary_rows = query_dicts(
        connection,
        "SELECT * FROM disease_ranking_panels ORDER BY release_pair, canonical_disease_id",
    )
    primary_by_key = {
        (str(row["release_pair"]), str(row["canonical_disease_id"])): row
        for (row) in (primary_rows)
    }
    stable_registry_rows = query_dicts(
        connection,
        f"SELECT * FROM fixed_ranking_panel_registry WHERE support_name = {sql_literal(configuration['ranking']['support_name'])} AND TRY_CAST(primary_panel_eligible AS BOOLEAN) ORDER BY release_pair, canonical_disease_id",
    )
    primary_score_rows = score_summary_rows(
        connection, "primary_pair_revisions", "PRIMARY_OVERALL_FIXED"
    )
    stable_score_rows = score_summary_rows(
        connection, "stable_entity_pair_revisions", "STABLE_ENTITY_FIXED"
    )
    domain_score_rows = score_summary_rows(
        connection, "evidence_domain_pair_revisions", "GENETIC_ASSOCIATION_FIXED"
    )
    stable_metric_rows = build_fixed_metric_rows(
        connection,
        "stable_entity_pair_revisions",
        configuration,
        "STABLE_ENTITY_FIXED",
        primary_by_key,
    )
    domain_metric_rows = build_fixed_metric_rows(
        connection,
        "evidence_domain_pair_revisions",
        configuration,
        "GENETIC_ASSOCIATION_FIXED",
        primary_by_key,
    )
    stable_native_rows = build_stable_native_rows(connection, configuration)
    create_changed_target_table(connection, primary_rows)
    diagnostics = create_diagnostic_tables(
        connection, project_root, diagnostic_paths, configuration
    )
    diagnostic_rows = build_diagnostic_summary_rows(connection)
    release_note_rows = build_release_note_rows(connection, configuration)
    area_memberships = stable_area_memberships(connection)
    comparison_rows = build_sensitivity_comparison_rows(
        primary_rows,
        stable_metric_rows,
        domain_metric_rows,
        stable_native_rows,
        primary_score_rows,
        stable_score_rows,
        domain_score_rows,
        area_memberships,
        configuration,
    )
    persistence_rows = build_persistence_rows(
        primary_rows,
        stable_metric_rows,
        domain_metric_rows,
        diagnostic_rows,
        stable_registry_rows,
        area_memberships,
    )
    datasets: dict[str, list[dict[str, Any]]] = {
        "stable_entity_pair_revisions": [],
        "stable_entity_score_summary": stable_score_rows,
        "stable_entity_ranking_metrics": stable_metric_rows,
        "stable_entity_native_metrics": stable_native_rows,
        "evidence_domain_pair_revisions": [],
        "evidence_domain_score_summary": domain_score_rows,
        "evidence_domain_ranking_metrics": domain_metric_rows,
        "diagnostic_target_context": [],
        "literature_clinical_diagnostics": diagnostic_rows,
        "release_note_alignment": release_note_rows,
        "sensitivity_comparison": comparison_rows,
        "persistence_registry": persistence_rows,
    }
    return datasets

def configure_connection(
    connection: Any,
    project_root: Path,
    configuration: Mapping[str, Any],
    build_identifier: str,
) -> Path:
    from .sensitivity_1 import resolve_project_path, sql_literal

    temporary_parent = resolve_project_path(
        project_root, str(configuration["resources"]["duckdb_temporary_directory"])
    )
    interim = (project_root / "data" / "interim").resolve()
    if (not temporary_parent.is_relative_to(interim)):
        raise CrossLayerSensitivityBuildError(
            "DuckDB temporary storage must remain inside data/interim"
        )
    temporary = (temporary_parent / build_identifier).resolve()
    if (not temporary.is_relative_to(temporary_parent)):
        raise CrossLayerSensitivityBuildError(
            "DuckDB build storage escapes its configured parent"
        )
    temporary.mkdir(parents = True, exist_ok = True)
    connection.execute(f"SET threads = {int(configuration['resources']['threads'])}")
    connection.execute(
        f"SET memory_limit = {sql_literal(configuration['resources']['memory_limit'])}"
    )
    connection.execute(f"SET temp_directory = {sql_literal(temporary.as_posix())}")
    connection.execute("SET preserve_insertion_order = false")
    return temporary
