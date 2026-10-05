from __future__ import annotations
from collections import defaultdict
from decimal import Decimal
from typing import Any
from typing import Mapping
from typing import Sequence
from .ranking_types import RankingSensitivityBuildError

def build_disease_registry(
    primary_fixed: Mapping[tuple[str, str], Mapping[str, Any]],
    primary_native: Mapping[tuple[str, str], Mapping[str, Any]],
    configuration: Mapping[str, Any],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for (key) in (sorted(primary_fixed)):
        fixed = primary_fixed[key]
        native = primary_native.get(key, {})
        mapping_context = (
            bool(fixed.get("fixed_top_ranked_mapping_only"))
            or int(native.get("mapping_affected_member_count", 0) or 0) > 0
        )
        row = {
            "release_pair": fixed["release_pair"],
            "old_release": fixed["old_release"],
            "new_release": fixed["new_release"],
            "canonical_disease_id": fixed["canonical_disease_id"],
            "canonical_label": fixed["canonical_label"],
            "canonical_therapeutic_area_ids": fixed["canonical_therapeutic_area_ids"],
            "fixed_target_count": fixed["fixed_target_count"],
            "old_native_target_count": native.get("old_native_target_count"),
            "new_native_target_count": native.get("new_native_target_count"),
            "kendall_tau_b": fixed["kendall_tau_b"],
            "spearman_correlation": fixed["spearman_correlation"],
            "median_normalised_rank_displacement": fixed[
                "median_normalised_rank_displacement"
            ],
            "p90_normalised_rank_displacement": fixed[
                "p90_normalised_rank_displacement"
            ],
            "maximum_normalised_rank_displacement": fixed[
                "maximum_normalised_rank_displacement"
            ],
            "fixed_top_5_jaccard": fixed["fixed_top_5_jaccard"],
            "fixed_top_10_jaccard": fixed["fixed_top_10_jaccard"],
            "fixed_top_20_jaccard": fixed["fixed_top_20_jaccard"],
            "native_top_5_jaccard": native.get("native_top_5_jaccard"),
            "native_top_10_jaccard": native.get("native_top_10_jaccard"),
            "native_top_20_jaccard": native.get("native_top_20_jaccard"),
            "old_fixed_top_ranked_target_ids": fixed["old_fixed_top_ranked_target_ids"],
            "new_fixed_top_ranked_target_ids": fixed["new_fixed_top_ranked_target_ids"],
            "fixed_top_ranked_set_changed": fixed["fixed_top_ranked_set_changed"],
            "fixed_top_ranked_sets_disjoint": fixed["fixed_top_ranked_sets_disjoint"],
            "fixed_top_ranked_change_reason": fixed["fixed_top_ranked_change_reason"],
            "old_native_top_ranked_target_ids": native.get(
                "old_native_top_ranked_target_ids"
            ),
            "new_native_top_ranked_target_ids": native.get(
                "new_native_top_ranked_target_ids"
            ),
            "native_top_ranked_set_changed": native.get(
                "native_top_ranked_set_changed"
            ),
            "native_top_ranked_sets_disjoint": native.get(
                "native_top_ranked_sets_disjoint"
            ),
            "native_top_ranked_change_reason": native.get(
                "native_top_ranked_change_reason"
            ),
            "entrant_count": native.get("entrant_member_count"),
            "exit_count": native.get("exit_member_count"),
            "mapping_affected_member_count": native.get(
                "mapping_affected_member_count"
            ),
            "native_top_ranked_entrant_count": native.get(
                "native_top_ranked_entrant_count"
            ),
            "native_top_ranked_exit_count": native.get("native_top_ranked_exit_count"),
            "native_top_ranked_mapping_affected_count": native.get(
                "native_top_ranked_mapping_affected_count"
            ),
            "mapping_stability": "MAPPING_CONTEXT_PRESENT"
            if (mapping_context)
            else "NO_MAPPING_CONTEXT_FLAG",
            "stable_entity_fixed_panel_eligible": fixed["stable_panel_eligible"],
            "release_native_primary_eligible": bool(native),
            "turnover_reason_is_descriptive_not_causal": True,
        }
        rows.append(row)
    return rows

def stable_area_memberships(
    connection: Any, configuration: Mapping[str, Any]
) -> tuple[
    dict[tuple[str, str], list[dict[str, Any]]],
    dict[tuple[str, str], str],
    dict[str, int],
]:
    from .ranking_1 import finite_decimal, parse_boolean, query_dicts

    noncanonical_contract_mismatch = int(
        connection.execute(
            "\n                SELECT COUNT(*)\n                FROM therapeutic_area_mapping\n                WHERE NULLIF(TRIM(CAST(canonical_disease_id AS VARCHAR)), '') IS NULL\n                  AND (\n                      COALESCE(UPPER(CAST(therapeutic_area_stability AS VARCHAR)), '') <> 'NOT_COMPARABLE'\n                      OR TRY_CAST(primary_area_contrast_eligible AS BOOLEAN) IS DISTINCT FROM false\n                      OR fractional_weight IS NOT NULL\n                  )\n            "
        ).fetchone()[0]
    )
    rows = query_dicts(
        connection,
        "\n            SELECT\n                release,\n                canonical_disease_id,\n                canonical_therapeutic_area_id,\n                canonical_therapeutic_area_label,\n                upper(CAST(therapeutic_area_stability AS VARCHAR)) AS therapeutic_area_stability,\n                primary_area_contrast_eligible,\n                CAST(fractional_weight AS DECIMAL(38,18)) AS fractional_weight\n            FROM therapeutic_area_mapping\n            WHERE NULLIF(TRIM(CAST(canonical_disease_id AS VARCHAR)), '') IS NOT NULL\n            ORDER BY release, canonical_disease_id, canonical_therapeutic_area_id\n        ",
    )
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for (row) in (rows):
        grouped[str(row["release"]), str(row["canonical_disease_id"])].append(row)
    areas = configuration["therapeutic_areas"]
    stable_status = str(areas["stable_area_status"])
    unstable_status = str(areas["unstable_area_status"])
    not_comparable_status = str(areas["not_comparable_area_status"])
    memberships: dict[tuple[str, str], list[dict[str, Any]]] = {}
    classifications: dict[tuple[str, str], str] = {}
    diagnostics = {
        "therapeutic_area_source_duplicate_membership_count": 0,
        "therapeutic_area_source_status_domain_mismatch_count": 0,
        "therapeutic_area_source_status_consistency_mismatch_count": 0,
        "therapeutic_area_source_eligibility_mismatch_count": 0,
        "therapeutic_area_source_noncanonical_row_contract_mismatch_count": noncanonical_contract_mismatch,
    }
    for (key, values) in (sorted(grouped.items())):
        area_ids = [str(value["canonical_therapeutic_area_id"]) for (value) in (values)]
        diagnostics["therapeutic_area_source_duplicate_membership_count"] += len(
            area_ids
        ) - len(set(area_ids))
        statuses = {
            str(value["therapeutic_area_stability"]).upper() for (value) in (values)
        }
        try:
            eligibility = {
                parse_boolean(value["primary_area_contrast_eligible"])
                for (value) in (values)
            }
        except ValueError as exc:
            raise RankingSensitivityBuildError(
                f"Invalid therapeutic-area eligibility for {key}"
            ) from exc
        diagnostics["therapeutic_area_source_status_domain_mismatch_count"] += sum(
            (
                status not in {"STABLE", "UNSTABLE", "NOT_COMPARABLE"}
                for (status) in (statuses)
            )
        )
        diagnostics["therapeutic_area_source_status_consistency_mismatch_count"] += (
            len(statuses) != 1
        )
        if (statuses == {"STABLE"} and eligibility == {True}):
            converted: list[dict[str, Any]] = []
            for (value) in (values):
                weight = finite_decimal(value["fractional_weight"])
                if (weight <= Decimal(0)):
                    raise RankingSensitivityBuildError(
                        "Therapeutic-area fractional weights must be positive"
                    )
                converted.append(dict(value) | {"fractional_weight": weight})
            memberships[key] = converted
            classifications[key] = stable_status
        elif (statuses == {"UNSTABLE"} and eligibility == {False}):
            classifications[key] = unstable_status
        elif (statuses == {"NOT_COMPARABLE"} and eligibility == {False}):
            classifications[key] = not_comparable_status
        else:
            diagnostics["therapeutic_area_source_eligibility_mismatch_count"] += 1
            classifications[key] = not_comparable_status
    return (memberships, classifications, diagnostics)

def build_therapeutic_area_summary(
    disease_rows: Sequence[Mapping[str, Any]],
    memberships: Mapping[tuple[str, str], Sequence[Mapping[str, Any]]],
    classifications: Mapping[tuple[str, str], str],
    configuration: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, int], dict[str, dict[str, Any]]]:
    from .ranking_1 import (
        finite_decimal,
        finite_float,
        parse_boolean,
        weighted_fraction,
        weighted_quantile,
    )

    grouped: dict[
        tuple[str, str, str], list[tuple[Mapping[str, Any], Decimal, str]]
    ] = defaultdict(list)
    diagnostics = {
        "therapeutic_area_panel_partition_mismatch_count": 0,
        "therapeutic_area_status_key_mismatch_count": 0,
        "therapeutic_area_status_count_mismatch_count": 0,
        "therapeutic_area_contribution_key_mismatch_count": 0,
        "therapeutic_area_equal_fraction_mismatch_count": 0,
        "therapeutic_area_weight_sum_mismatch_count": 0,
        "therapeutic_area_fractional_total_mismatch_count": 0,
    }
    pair_labels = [str(pair["label"]) for (pair) in (configuration["release_pairs"])]
    area_config = configuration["therapeutic_areas"]
    status_values = [
        str(area_config["stable_area_status"]),
        str(area_config["unstable_area_status"]),
        str(area_config["not_comparable_area_status"]),
        str(area_config["no_explicit_mapping_area_status"]),
    ]
    coverage: dict[str, dict[str, Any]] = {
        pair: {
            "primary_panel_count": 0,
            "stable_area_eligible_panel_count": 0,
            "outside_stable_area_panel_count": 0,
            "fractional_weight_total": Decimal(0),
            "therapeutic_area_count": 0,
            "status_counts": {status: 0 for (status) in (status_values)},
            "expected_status_counts": {status: 0 for (status) in (status_values)},
        }
        for (pair) in (pair_labels)
    }
    all_panel_keys: set[tuple[str, str]] = set()
    eligible_panel_keys: set[tuple[str, str]] = set()
    contributed_panel_keys: set[tuple[str, str]] = set()
    observed_status_keys: dict[str, set[tuple[str, str]]] = {
        status: set() for (status) in (status_values)
    }
    expected_status_keys: dict[str, set[tuple[str, str]]] = {
        status: set() for (status) in (status_values)
    }
    stable_status = str(area_config["stable_area_status"])
    no_mapping_status = str(area_config["no_explicit_mapping_area_status"])
    tolerance = finite_decimal(area_config["fractional_weight_tolerance"])
    for (row) in (disease_rows):
        pair = str(row["release_pair"])
        disease = str(row["canonical_disease_id"])
        panel_key = (pair, disease)
        all_panel_keys.add(panel_key)
        coverage[pair]["primary_panel_count"] += 1
        observed_status = str(row["therapeutic_area_membership_status"])
        if (observed_status in observed_status_keys):
            observed_status_keys[observed_status].add(panel_key)
            coverage[pair]["status_counts"][observed_status] += 1
        key = (str(row["old_release"]), disease)
        expected_status = str(classifications.get(key, no_mapping_status))
        if (expected_status in expected_status_keys):
            expected_status_keys[expected_status].add(panel_key)
            coverage[pair]["expected_status_counts"][expected_status] += 1
        areas = list(memberships.get(key, []))
        if (not areas):
            coverage[pair]["outside_stable_area_panel_count"] += 1
            continue
        eligible_panel_keys.add(panel_key)
        contributed_panel_keys.add(panel_key)
        coverage[pair]["stable_area_eligible_panel_count"] += 1
        weights = [finite_decimal(area["fractional_weight"]) for (area) in (areas)]
        expected_weight = Decimal(1) / Decimal(len(areas))
        diagnostics["therapeutic_area_equal_fraction_mismatch_count"] += sum(
            (abs(weight - expected_weight) > tolerance for (weight) in (weights))
        )
        weight_sum = sum(weights, Decimal(0))
        coverage[pair]["fractional_weight_total"] += weight_sum
        if (abs(weight_sum - Decimal(1)) > tolerance):
            diagnostics["therapeutic_area_weight_sum_mismatch_count"] += 1
        for (area) in (areas):
            grouped[
                pair,
                str(area["canonical_therapeutic_area_id"]),
                str(area["canonical_therapeutic_area_label"]),
            ].append((row, finite_decimal(area["fractional_weight"]), disease))
    observed_union: set[tuple[str, str]] = set()
    observed_intersection_count = 0
    for (status) in (status_values):
        observed_intersection_count += len(
            observed_union & observed_status_keys[status]
        )
        observed_union |= observed_status_keys[status]
    diagnostics["therapeutic_area_panel_partition_mismatch_count"] = (
        len(observed_union ^ all_panel_keys) + observed_intersection_count
    )
    diagnostics["therapeutic_area_status_key_mismatch_count"] = sum(
        (
            len(observed_status_keys[status] ^ expected_status_keys[status])
            for (status) in (status_values)
        )
    )
    diagnostics["therapeutic_area_status_count_mismatch_count"] = sum(
        (
            abs(
                int(coverage[pair]["status_counts"][status])
                - int(coverage[pair]["expected_status_counts"][status])
            )
            for (pair) in (pair_labels)
            for (status) in (status_values)
        )
    )
    diagnostics["therapeutic_area_contribution_key_mismatch_count"] = len(
        eligible_panel_keys ^ contributed_panel_keys
    ) + len(eligible_panel_keys ^ expected_status_keys[stable_status])
    diagnostics["therapeutic_area_fractional_total_mismatch_count"] = sum(
        (
            abs(
                coverage[pair]["fractional_weight_total"]
                - Decimal(coverage[pair]["stable_area_eligible_panel_count"])
            )
            > tolerance
            for (pair) in (pair_labels)
        )
    )
    output: list[dict[str, Any]] = []
    low_overlap_threshold = finite_float(
        configuration["signal_assessment"]["native_top_10_jaccard_below"]
    )
    minimum_panels = int(
        configuration["therapeutic_areas"]["well_populated_minimum_panels"]
    )
    for ((release_pair, area_id, area_label), values) in (sorted(grouped.items())):
        rows = [item[0] for (item) in (values)]
        weights = [item[1] for (item) in (values)]
        diseases = {item[2] for (item) in (values)}
        native_available = [
            index
            for ((index, row)) in (enumerate(rows))
            if (row.get("native_top_10_jaccard") not in {None, ""})
        ]
        native_values = [
            finite_float(rows[index]["native_top_10_jaccard"])
            for (index) in (native_available)
        ]
        native_weights = [weights[index] for (index) in (native_available)]
        fixed_changed = [
            parse_boolean(row["fixed_top_ranked_set_changed"]) for (row) in (rows)
        ]
        fixed_top_10_changed = [
            finite_float(row["fixed_top_10_jaccard"]) < 1.0 for (row) in (rows)
        ]
        native_top_10_changed = [value < 1.0 for (value) in (native_values)]
        native_low_overlap = [
            value < low_overlap_threshold for (value) in (native_values)
        ]
        native_top_changed_values = [
            parse_boolean(rows[index]["native_top_ranked_set_changed"])
            for (index) in (native_available)
        ]
        entrant_driven = [
            native_top_changed_values[position]
            and str(rows[index]["native_top_ranked_change_reason"])
            in {"ENTRANT", "MIXED"}
            and (int(rows[index]["native_top_ranked_entrant_count"]) > 0)
            for ((position, index)) in (enumerate(native_available))
        ]
        mapping_affected = [
            native_top_changed_values[position]
            and (
                str(rows[index]["native_top_ranked_change_reason"])
                in {"MAPPING_AFFECTED", "MIXED"}
                or int(rows[index]["native_top_ranked_mapping_affected_count"]) > 0
            )
            for ((position, index)) in (enumerate(native_available))
        ]
        tau_values = [finite_float(row["kendall_tau_b"]) for (row) in (rows)]
        displacement_values = [
            finite_float(row["p90_normalised_rank_displacement"]) for (row) in (rows)
        ]
        output.append(
            {
                "release_pair": release_pair,
                "therapeutic_area_id": area_id,
                "therapeutic_area_label": area_label,
                "eligible_disease_panel_count": len(diseases),
                "fractional_panel_count": sum(weights),
                "native_top_10_estimable_panel_count": len(native_available),
                "kendall_tau_b_q25": weighted_quantile(tau_values, weights, 0.25),
                "median_kendall_tau_b": weighted_quantile(tau_values, weights, 0.5),
                "kendall_tau_b_q75": weighted_quantile(tau_values, weights, 0.75),
                "median_panel_p90_displacement": weighted_quantile(
                    displacement_values, weights, 0.5
                ),
                "q90_panel_p90_displacement": weighted_quantile(
                    displacement_values, weights, 0.9
                ),
                "fixed_top_10_change_fraction": weighted_fraction(
                    fixed_top_10_changed, weights
                ),
                "native_top_10_change_fraction": weighted_fraction(
                    native_top_10_changed, native_weights
                ),
                "fixed_top_ranked_set_changed_count": sum(fixed_changed),
                "fixed_top_ranked_set_changed_fraction": weighted_fraction(
                    fixed_changed, weights
                ),
                "native_top_10_low_overlap_count": sum(native_low_overlap),
                "native_top_10_low_overlap_fraction": weighted_fraction(
                    native_low_overlap, native_weights
                ),
                "native_top_ranked_set_change_fraction": weighted_fraction(
                    native_top_changed_values, native_weights
                ),
                "release_native_entrant_driven_change_fraction": weighted_fraction(
                    entrant_driven, native_weights
                ),
                "mapping_affected_change_fraction": weighted_fraction(
                    mapping_affected, native_weights
                ),
                "area_comparison_status": "WELL_POPULATED"
                if (len(diseases) >= minimum_panels)
                else configuration["therapeutic_areas"]["small_area_status"].upper(),
                "area_rank": None,
            }
        )
    for (pair) in (pair_labels):
        coverage[pair]["therapeutic_area_count"] = sum(
            (str(row["release_pair"]) == pair for (row) in (output))
        )
    return (output, diagnostics, coverage)
