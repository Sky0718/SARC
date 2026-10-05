from __future__ import annotations
import csv
import math
from pathlib import Path
from typing import Any
from typing import Mapping
from typing import Sequence
from .ranking_types import RankingSensitivityBuildError

def output_specifications() -> dict[str, dict[str, Any]]:
    fixed_metric_fields = [
        "release_pair",
        "old_release",
        "new_release",
        "canonical_disease_id",
        "canonical_label",
        "canonical_therapeutic_area_ids",
        "stable_panel_eligible",
        "fixed_target_count",
        "old_distinct_score_count",
        "new_distinct_score_count",
        "kendall_tau_b",
        "spearman_correlation",
        "median_normalised_rank_displacement",
        "p90_normalised_rank_displacement",
        "maximum_normalised_rank_displacement",
    ]
    for (cutoff) in ([5, 10, 20]):
        fixed_metric_fields.extend(
            [
                f"fixed_top_{cutoff}_jaccard",
                f"old_fixed_top_{cutoff}_effective_size",
                f"new_fixed_top_{cutoff}_effective_size",
                f"fixed_top_{cutoff}_changed",
            ]
        )
    fixed_metric_fields.extend(
        [
            "old_fixed_top_ranked_target_ids",
            "new_fixed_top_ranked_target_ids",
            "old_fixed_top_ranked_set_size",
            "new_fixed_top_ranked_set_size",
            "fixed_top_ranked_jaccard",
            "fixed_top_ranked_set_changed",
            "fixed_top_ranked_sets_disjoint",
            "fixed_top_ranked_change_reason",
            "fixed_top_ranked_mapping_only",
            "panel_population",
            "minimum_target_threshold",
        ]
    )
    native_fields = [
        "release_pair",
        "old_release",
        "new_release",
        "canonical_disease_id",
        "canonical_label",
        "canonical_therapeutic_area_ids",
        "primary_fixed_anchored_eligible",
        "all_native_supplementary_eligible",
        "old_native_target_count",
        "new_native_target_count",
        "old_native_distinct_score_count",
        "new_native_distinct_score_count",
        "persistent_member_count",
        "entrant_member_count",
        "exit_member_count",
        "mapping_affected_member_count",
    ]
    for (cutoff) in ([5, 10, 20]):
        native_fields.extend(
            [
                f"native_top_{cutoff}_jaccard",
                f"old_native_top_{cutoff}_effective_size",
                f"new_native_top_{cutoff}_effective_size",
                f"native_top_{cutoff}_changed",
                f"native_top_{cutoff}_change_reason",
                f"native_top_{cutoff}_entrant_count",
                f"native_top_{cutoff}_exit_count",
                f"native_top_{cutoff}_mapping_affected_count",
            ]
        )
    native_fields.extend(
        [
            "old_native_top_ranked_target_ids",
            "new_native_top_ranked_target_ids",
            "old_native_top_ranked_set_size",
            "new_native_top_ranked_set_size",
            "native_top_ranked_jaccard",
            "native_top_ranked_set_changed",
            "native_top_ranked_sets_disjoint",
            "native_top_ranked_change_reason",
            "native_top_ranked_entrant_count",
            "native_top_ranked_exit_count",
            "native_top_ranked_mapping_affected_count",
        ]
    )
    disease_fields = [
        "release_pair",
        "old_release",
        "new_release",
        "canonical_disease_id",
        "canonical_label",
        "canonical_therapeutic_area_ids",
        "fixed_target_count",
        "old_native_target_count",
        "new_native_target_count",
        "kendall_tau_b",
        "spearman_correlation",
        "median_normalised_rank_displacement",
        "p90_normalised_rank_displacement",
        "maximum_normalised_rank_displacement",
        "fixed_top_5_jaccard",
        "fixed_top_10_jaccard",
        "fixed_top_20_jaccard",
        "native_top_5_jaccard",
        "native_top_10_jaccard",
        "native_top_20_jaccard",
        "old_fixed_top_ranked_target_ids",
        "new_fixed_top_ranked_target_ids",
        "fixed_top_ranked_set_changed",
        "fixed_top_ranked_sets_disjoint",
        "fixed_top_ranked_change_reason",
        "old_native_top_ranked_target_ids",
        "new_native_top_ranked_target_ids",
        "native_top_ranked_set_changed",
        "native_top_ranked_sets_disjoint",
        "native_top_ranked_change_reason",
        "entrant_count",
        "exit_count",
        "mapping_affected_member_count",
        "native_top_ranked_entrant_count",
        "native_top_ranked_exit_count",
        "native_top_ranked_mapping_affected_count",
        "mapping_stability",
        "therapeutic_area_membership_status",
        "stable_entity_fixed_panel_eligible",
        "release_native_primary_eligible",
        "turnover_reason_is_descriptive_not_causal",
    ]
    return {
        "disease_ranking_panels": {
            "fields": disease_fields,
            "order": ["release_pair", "canonical_disease_id"],
            "parquet": False,
        },
        "fixed_target_rank_displacement": {
            "fields": [
                "release_pair",
                "old_release",
                "new_release",
                "canonical_disease_id",
                "canonical_target_id",
                "old_score",
                "new_score",
                "old_rank",
                "new_rank",
                "normalised_absolute_rank_displacement",
                "stable_entity_member",
                "mapping_affected",
                "eligible_min_20",
                "eligible_primary_30",
                "eligible_min_50",
            ],
            "types": [
                "VARCHAR",
                "VARCHAR",
                "VARCHAR",
                "VARCHAR",
                "VARCHAR",
                "DOUBLE",
                "DOUBLE",
                "DOUBLE",
                "DOUBLE",
                "DOUBLE",
                "BOOLEAN",
                "BOOLEAN",
                "BOOLEAN",
                "BOOLEAN",
                "BOOLEAN",
            ],
            "order": ["release_pair", "canonical_disease_id", "canonical_target_id"],
            "parquet": True,
        },
        "fixed_ranking_metrics": {
            "fields": fixed_metric_fields,
            "order": [
                "release_pair",
                "minimum_target_threshold",
                "canonical_disease_id",
            ],
            "parquet": False,
        },
        "release_native_top_set_metrics": {
            "fields": native_fields,
            "order": ["release_pair", "canonical_disease_id"],
            "parquet": False,
        },
        "top_set_members": {
            "fields": [
                "analysis_type",
                "release_pair",
                "canonical_disease_id",
                "set_definition",
                "cutoff",
                "canonical_target_id",
                "old_member",
                "new_member",
                "membership_state",
                "context_class",
                "mapping_affected",
            ],
            "types": [
                "VARCHAR",
                "VARCHAR",
                "VARCHAR",
                "VARCHAR",
                "BIGINT",
                "VARCHAR",
                "BOOLEAN",
                "BOOLEAN",
                "VARCHAR",
                "VARCHAR",
                "BOOLEAN",
            ],
            "order": [
                "analysis_type",
                "release_pair",
                "canonical_disease_id",
                "set_definition",
                "canonical_target_id",
            ],
            "parquet": True,
        },
        "top_ranked_target_sets": {
            "fields": [
                "analysis_type",
                "release_pair",
                "canonical_disease_id",
                "old_top_ranked_target_ids",
                "new_top_ranked_target_ids",
                "old_set_size",
                "new_set_size",
                "jaccard",
                "changed",
                "disjoint",
                "change_reason",
                "descriptive_not_causal",
            ],
            "order": ["analysis_type", "release_pair", "canonical_disease_id"],
            "parquet": False,
        },
        "turnover_reason_decomposition": {
            "fields": [
                "release_pair",
                "set_definition",
                "population",
                "reason_class",
                "estimable_panel_count",
                "changed_panel_count",
                "reason_panel_count",
                "reason_panel_fraction",
                "reason_changed_panel_fraction",
                "descriptive_not_causal",
            ],
            "order": ["release_pair", "set_definition", "population", "reason_class"],
            "parquet": False,
        },
        "therapeutic_area_ranking_summary": {
            "fields": [
                "release_pair",
                "therapeutic_area_id",
                "therapeutic_area_label",
                "eligible_disease_panel_count",
                "fractional_panel_count",
                "native_top_10_estimable_panel_count",
                "kendall_tau_b_q25",
                "median_kendall_tau_b",
                "kendall_tau_b_q75",
                "median_panel_p90_displacement",
                "q90_panel_p90_displacement",
                "fixed_top_10_change_fraction",
                "native_top_10_change_fraction",
                "fixed_top_ranked_set_changed_count",
                "fixed_top_ranked_set_changed_fraction",
                "native_top_10_low_overlap_count",
                "native_top_10_low_overlap_fraction",
                "native_top_ranked_set_change_fraction",
                "release_native_entrant_driven_change_fraction",
                "mapping_affected_change_fraction",
                "area_comparison_status",
                "area_rank",
            ],
            "order": ["release_pair", "therapeutic_area_id"],
            "parquet": False,
        },
        "sensitivity_signal_assessment": {
            "fields": [
                "release_pair",
                "population",
                "criterion",
                "numerator",
                "denominator",
                "value",
                "comparison",
                "threshold",
                "passed",
                "validation_role",
                "notes",
            ],
            "order": ["release_pair", "criterion"],
            "parquet": False,
        },
    }

def serialise_csv_cell(value: Any) -> Any:
    if (value is None):
        return ""
    if (isinstance(value, bool)):
        return str(value).lower()
    if (isinstance(value, float)):
        if (not math.isfinite(value)):
            raise RankingSensitivityBuildError("Cannot export a non-finite value")
        return format(value, ".17g")
    return value

def sort_key(row: Mapping[str, Any], fields: Sequence[str]) -> tuple[str, ...]:
    return tuple(
        (
            "" if (row.get(field) is None) else str(row.get(field))
            for (field) in (fields)
        )
    )

def write_csv_dataset(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    fields: Sequence[str],
    order: Sequence[str],
) -> None:
    path.parent.mkdir(parents = True, exist_ok = True)
    with path.open("w", encoding = "utf-8", newline = "") as handle:
        writer = csv.DictWriter(
            handle, fieldnames = fields, extrasaction = "raise", lineterminator = "\n"
        )
        writer.writeheader()
        for (row) in (sorted(rows, key = lambda item: sort_key(item, order))):
            writer.writerow(
                {field: serialise_csv_cell(row.get(field)) for (field) in (fields)}
            )
