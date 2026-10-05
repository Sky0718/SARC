from __future__ import annotations
from pathlib import Path
from typing import Any
from typing import Mapping
from typing import Sequence
from .score_revision_types import ScoreRevisionBuildError

def parse_boolean(value: Any) -> bool:
    if (isinstance(value, bool)):
        return value
    if (isinstance(value, int) and value in {0, 1}):
        return bool(value)
    normalised = str(value).strip().lower()
    if (normalised in {"true", "1", "yes"}):
        return True
    if (normalised in {"false", "0", "no"}):
        return False
    raise ValueError(f"Invalid Boolean value: {value}")

def configure_duckdb_resources(
    connection: Any, temporary_root: Path, configuration: Mapping[str, Any]
) -> Path:
    resources = configuration["resources"]
    temp_directory = (
        temporary_root / str(resources["duckdb_temp_subdirectory"])
    ).resolve()
    if (temp_directory.parent != temporary_root.resolve()):
        raise ScoreRevisionBuildError(
            "DuckDB temporary directory escaped the Score revision build directory"
        )
    temp_directory.mkdir(parents = True, exist_ok = True)
    connection.execute(f"SET threads = {int(resources['threads'])}")
    connection.execute(f"SET memory_limit = {sql_literal(resources['memory_limit'])}")
    connection.execute(f"SET temp_directory = {sql_literal(temp_directory.as_posix())}")
    connection.execute(
        f"SET max_temp_directory_size = {sql_literal(resources['maximum_temp_directory_size'])}"
    )
    connection.execute("SET preserve_insertion_order = false")
    return temp_directory

def sql_literal(value: Any) -> str:
    return "'" + str(value).replace("'", "''") + "'"

def sql_identifier(value: str) -> str:
    return '"' + value.replace('"', '""') + '"'

def score_tier_case(expression: str, tiers: Sequence[Mapping[str, Any]]) -> str:
    clauses: list[str] = []
    for (tier) in (tiers):
        lower_operator = ">=" if (parse_boolean(tier["lower_inclusive"])) else ">"
        upper_operator = "<=" if (parse_boolean(tier["upper_inclusive"])) else "<"
        clauses.append(
            f"WHEN {expression} {lower_operator} {float(tier['lower'])} AND {expression} {upper_operator} {float(tier['upper'])} THEN {sql_literal(tier['label'])}"
        )
    return "CASE " + " ".join(clauses) + " ELSE 'unclassified' END"

def table_columns(connection: Any, relation_sql: str) -> set[str]:
    rows = connection.execute(f"DESCRIBE SELECT * FROM {relation_sql}").fetchall()
    return {str(row[0]) for (row) in (rows)}

def require_columns(
    connection: Any, relation_sql: str, expected: Sequence[str], label: str
) -> None:
    observed = table_columns(connection, relation_sql)
    missing = sorted(set(expected) - observed)
    if (missing):
        raise ScoreRevisionBuildError(f"{label} columns are absent: {missing}")

def create_input_views(
    connection: Any, paths: Mapping[str, Path], configuration: Mapping[str, Any]
) -> None:
    sources = {
        "support_construction_three_release": (
            "three_release_score_support",
            "read_parquet",
        ),
        "support_construction_transitions": (
            "association_state_transitions",
            "read_parquet",
        ),
        "support_construction_entrant_exit": ("entrant_exit_table", "read_parquet"),
        "support_construction_areas": ("therapeutic_area_mapping", "read_csv_auto"),
        "support_construction_native_members": (
            "release_native_panel_members",
            "read_parquet",
        ),
        "support_construction_evidence_domain": (
            "evidence_domain_support",
            "read_parquet",
        ),
    }
    for (view_name, (input_name, reader)) in (sources.items()):
        path_literal = sql_literal(paths[input_name].as_posix())
        reader_call = (
            f"{reader}({path_literal}, header=true)"
            if (reader == "read_csv_auto")
            else f"{reader}({path_literal})"
        )
        expected = configuration["input_contract"][input_name]
        require_columns(connection, reader_call, expected, input_name)
        connection.execute(
            f"CREATE VIEW {sql_identifier(view_name)} AS SELECT * FROM {reader_call}"
        )

def create_stable_area_table(connection: Any, configuration: Mapping[str, Any]) -> None:
    stability = sql_literal(configuration["therapeutic_areas"]["stability_value"])
    release_count = len(
        {pair["old_release"] for (pair) in (configuration["release_pairs"])}
        | {configuration["release_pairs"][-1]["new_release"]}
    )
    connection.execute(
        f"\n        CREATE TABLE stable_areas AS\n        SELECT\n            canonical_disease_id,\n            canonical_therapeutic_area_id AS therapeutic_area_id,\n            MAX(canonical_therapeutic_area_label) AS therapeutic_area_label,\n            MAX(CAST(fractional_weight AS DECIMAL(38,18))) AS fractional_weight\n        FROM support_construction_areas\n        WHERE UPPER(CAST(therapeutic_area_stability AS VARCHAR)) = {stability}\n          AND CAST(primary_area_contrast_eligible AS BOOLEAN)\n        GROUP BY canonical_disease_id, canonical_therapeutic_area_id\n        HAVING COUNT(DISTINCT release) = {release_count}\n           AND ABS(MAX(CAST(fractional_weight AS DOUBLE)) - MIN(CAST(fractional_weight AS DOUBLE))) <= 1e-12\n        "
    )

def create_pair_revision_table(
    connection: Any, configuration: Mapping[str, Any]
) -> None:
    primary_queries: list[str] = []
    for (pair) in (configuration["release_pairs"]):
        primary_queries.append(
            f"\n            SELECT\n                {sql_literal(pair['label'])} AS release_pair,\n                {sql_literal(configuration['supports']['primary'])} AS support_stratum,\n                s.canonical_disease_id,\n                s.canonical_target_id,\n                CAST(s.{sql_identifier(pair['old_score_field'])} AS DOUBLE) AS score_old,\n                CAST(s.{sql_identifier(pair['new_score_field'])} AS DOUBLE) AS score_new,\n                CAST(s.{sql_identifier(pair['old_evidence_count_field'])} AS BIGINT) AS evidence_count_old,\n                CAST(s.{sql_identifier(pair['new_evidence_count_field'])} AS BIGINT) AS evidence_count_new,\n                CAST(s.stable_entity_eligible AS BOOLEAN) AS stable_entity_eligible,\n                CAST(t.raw_identifier_mapping_affected AS BOOLEAN) AS raw_identifier_mapping_affected,\n                CAST(t.transition_class AS VARCHAR) AS transition_class\n            FROM support_construction_three_release s\n            INNER JOIN support_construction_transitions t\n              ON t.release_pair = {sql_literal(pair['label'])}\n             AND t.canonical_disease_id = s.canonical_disease_id\n             AND t.canonical_target_id = s.canonical_target_id\n            "
        )
    pair_specific = f"\n        SELECT\n            CAST(release_pair AS VARCHAR) AS release_pair,\n            {sql_literal(configuration['supports']['pair_specific_sensitivity'])} AS support_stratum,\n            canonical_disease_id,\n            canonical_target_id,\n            CAST(old_score AS DOUBLE) AS score_old,\n            CAST(new_score AS DOUBLE) AS score_new,\n            CAST(old_evidence_count AS BIGINT) AS evidence_count_old,\n            CAST(new_evidence_count AS BIGINT) AS evidence_count_new,\n            CAST(stable_entity_eligible AS BOOLEAN) AS stable_entity_eligible,\n            CAST(raw_identifier_mapping_affected AS BOOLEAN) AS raw_identifier_mapping_affected,\n            CAST(transition_class AS VARCHAR) AS transition_class\n        FROM support_construction_transitions\n        WHERE CAST(old_primary_eligible AS BOOLEAN)\n          AND CAST(new_primary_eligible AS BOOLEAN)\n          AND UPPER(CAST(old_state AS VARCHAR)) = 'PRESENT_POSITIVE'\n          AND UPPER(CAST(new_state AS VARCHAR)) = 'PRESENT_POSITIVE'\n    "
    combined = " UNION ALL ".join(primary_queries + [pair_specific])
    old_tier_case = score_tier_case("score_old", configuration["score_tiers"])
    new_tier_case = score_tier_case("score_new", configuration["score_tiers"])
    tolerance = float(configuration["scores"]["descriptive_tolerance"])
    calibration = configuration["implementation_scale_diagnostic"]
    calibration_pair = sql_literal(calibration["release_pair"])
    calibration_factor = float(calibration["old_score_multiplier"])
    connection.execute(
        f"\n        CREATE TABLE pair_score_revisions AS\n        WITH combined AS ({combined}),\n        calculated AS (\n            SELECT\n                *,\n                score_new - score_old AS signed_change,\n                ABS(score_new - score_old) AS absolute_change,\n                score_old = score_new AS exact_unchanged,\n                ABS(score_new - score_old) <= {tolerance} AS near_unchanged,\n                CASE\n                    WHEN evidence_count_old IS NULL OR evidence_count_new IS NULL THEN 'missing'\n                    WHEN evidence_count_old = evidence_count_new THEN 'unchanged'\n                    ELSE 'changed'\n                END AS evidence_count_status,\n                CASE\n                    WHEN raw_identifier_mapping_affected THEN 'raw_identifier_mapping_affected'\n                    ELSE 'canonical_identifier_stable'\n                END AS mapping_context\n            FROM combined\n        ), diagnostic_base AS (\n            SELECT\n                *,\n                release_pair = {calibration_pair} AS calibration_applicable,\n                CASE WHEN release_pair = {calibration_pair} THEN {calibration_factor} ELSE NULL END AS calibration_factor,\n                CASE WHEN release_pair = {calibration_pair} THEN score_old * {calibration_factor} ELSE NULL END AS calibrated_score_old\n            FROM calculated\n        ), diagnostic_changes AS (\n            SELECT\n                *,\n                CASE WHEN calibration_applicable THEN score_new - calibrated_score_old ELSE NULL END AS implementation_adjusted_signed_change,\n                CASE WHEN calibration_applicable THEN ABS(score_new - calibrated_score_old) ELSE NULL END AS implementation_adjusted_absolute_change\n            FROM diagnostic_base\n        )\n        SELECT\n            *,\n            {old_tier_case} AS old_score_tier,\n            {new_tier_case} AS new_score_tier\n        FROM diagnostic_changes\n        "
    )

def create_entrant_exit_context_table(
    connection: Any, configuration: Mapping[str, Any]
) -> None:
    context = configuration["entrant_exit_context"]
    tier_case = score_tier_case("applicable_score", configuration["score_tiers"])
    connection.execute(
        f"\n        CREATE TABLE entrant_exit_context AS\n        WITH classified AS (\n            SELECT\n                e.*,\n                old_native.average_native_rank AS computed_old_native_rank_context,\n                new_native.average_native_rank AS computed_new_native_rank_context,\n                old_native.tie_aware_top_rank_member AS computed_old_top_rank_context,\n                new_native.tie_aware_top_rank_member AS computed_new_top_rank_context,\n                CASE\n                    WHEN CAST(e.raw_identifier_mapping_affected AS BOOLEAN)\n                      OR UPPER(CAST(e.transition_class AS VARCHAR)) LIKE '%MAPPING%'\n                        THEN {sql_literal(context['mapping_affected_class'])}\n                    WHEN UPPER(CAST(e.old_state AS VARCHAR)) = 'PRESENT_ZERO'\n                      OR UPPER(CAST(e.new_state AS VARCHAR)) = 'PRESENT_ZERO'\n                        THEN {sql_literal(context['present_zero_context_class'])}\n                    WHEN UPPER(CAST(e.old_state AS VARCHAR)) = 'PRESENT_NULL'\n                      OR UPPER(CAST(e.new_state AS VARCHAR)) = 'PRESENT_NULL'\n                        THEN {sql_literal(context['present_null_context_class'])}\n                    WHEN UPPER(CAST(e.old_state AS VARCHAR)) = 'SCHEMA_INCOMPARABLE'\n                      OR UPPER(CAST(e.new_state AS VARCHAR)) = 'SCHEMA_INCOMPARABLE'\n                        THEN {sql_literal(context['schema_incomparable_context_class'])}\n                    WHEN UPPER(CAST(e.old_state AS VARCHAR)) = 'ABSENT_PAIR'\n                     AND UPPER(CAST(e.new_state AS VARCHAR)) = 'PRESENT_POSITIVE'\n                        THEN {sql_literal(context['entrant_class'])}\n                    WHEN UPPER(CAST(e.old_state AS VARCHAR)) = 'PRESENT_POSITIVE'\n                     AND UPPER(CAST(e.new_state AS VARCHAR)) = 'ABSENT_PAIR'\n                        THEN {sql_literal(context['exit_class'])}\n                    ELSE {sql_literal(context['remaining_other_context_class'])}\n                END AS context_class\n            FROM support_construction_entrant_exit e\n            LEFT JOIN native_rank_context old_native\n              ON old_native.release = e.old_release\n             AND old_native.canonical_disease_id = e.canonical_disease_id\n             AND old_native.canonical_target_id = e.canonical_target_id\n            LEFT JOIN native_rank_context new_native\n              ON new_native.release = e.new_release\n             AND new_native.canonical_disease_id = e.canonical_disease_id\n             AND new_native.canonical_target_id = e.canonical_target_id\n        ), applicable AS (\n            SELECT\n                *,\n                CASE\n                    WHEN context_class = {sql_literal(context['entrant_class'])} THEN CAST(new_score AS DOUBLE)\n                    WHEN context_class = {sql_literal(context['exit_class'])} THEN CAST(old_score AS DOUBLE)\n                    ELSE COALESCE(CAST(new_score AS DOUBLE), CAST(old_score AS DOUBLE))\n                END AS applicable_score,\n                CASE\n                    WHEN context_class = {sql_literal(context['entrant_class'])} THEN computed_new_native_rank_context\n                    WHEN context_class = {sql_literal(context['exit_class'])} THEN computed_old_native_rank_context\n                    ELSE NULL\n                END AS applicable_native_rank,\n                CASE\n                    WHEN context_class = {sql_literal(context['entrant_class'])} THEN computed_new_top_rank_context\n                    WHEN context_class = {sql_literal(context['exit_class'])} THEN computed_old_top_rank_context\n                    ELSE NULL\n                END AS applicable_top_rank_context\n            FROM classified\n        )\n        SELECT *,\n            CASE WHEN applicable_score IS NULL THEN 'score_not_defined' ELSE {tier_case} END AS score_tier\n        FROM applicable\n        "
    )

def create_native_rank_context_table(
    connection: Any, configuration: Mapping[str, Any]
) -> None:
    rank_cutoff = int(configuration["entrant_exit_context"]["high_rank_cutoff"])
    connection.execute(
        f"\n        CREATE TABLE native_rank_context AS\n        WITH ranked AS (\n            SELECT\n                release,\n                canonical_disease_id,\n                canonical_target_id,\n                CAST(score AS DOUBLE) AS score,\n                RANK() OVER (\n                    PARTITION BY release, canonical_disease_id\n                    ORDER BY CAST(score AS DOUBLE) DESC\n                ) AS first_tied_rank,\n                COUNT(*) OVER (\n                    PARTITION BY release, canonical_disease_id, CAST(score AS DOUBLE)\n                ) AS tied_target_count\n            FROM support_construction_native_members\n        )\n        SELECT\n            release,\n            canonical_disease_id,\n            canonical_target_id,\n            score,\n            first_tied_rank,\n            tied_target_count,\n            first_tied_rank + (tied_target_count - 1) / 2.0 AS average_native_rank,\n            first_tied_rank <= {rank_cutoff} AS tie_aware_top_rank_member\n        FROM ranked\n        "
    )

def selected_endpoint_values_sql(configuration: Mapping[str, Any]) -> str:
    values: list[str] = []
    for (endpoint) in (configuration["evidence_domain_context"]["selected_endpoints"]):
        values.append(
            "("
            + ",".join(
                [
                    sql_literal(endpoint["endpoint_level"]),
                    sql_literal(endpoint["endpoint_id"]),
                    sql_literal(endpoint["role"]),
                ]
            )
            + ")"
        )
    return ",".join(values)

def endpoint_filter_sql(configuration: Mapping[str, Any]) -> str:
    aliases = configuration["evidence_domain_context"]["endpoint_level_aliases"]
    clauses: list[str] = []
    for (endpoint) in (configuration["evidence_domain_context"]["selected_endpoints"]):
        accepted_levels = aliases[endpoint["endpoint_level"]]
        level_sql = ",".join((sql_literal(value) for (value) in (accepted_levels)))
        clauses.append(
            f"(endpoint_id = {sql_literal(endpoint['endpoint_id'])} AND LOWER(CAST(endpoint_level AS VARCHAR)) IN ({level_sql}))"
        )
    return " OR ".join(clauses)

def create_evidence_domain_context_tables(
    connection: Any, configuration: Mapping[str, Any]
) -> None:
    endpoint_values = selected_endpoint_values_sql(configuration)
    endpoint_filter = endpoint_filter_sql(configuration)
    pair_queries: list[str] = []
    for (pair) in (configuration["release_pairs"]):
        pair_queries.append(
            f"\n            SELECT\n                {sql_literal(pair['label'])} AS release_pair,\n                d.endpoint_level,\n                d.endpoint_id,\n                d.endpoint_role,\n                e.canonical_disease_id,\n                e.canonical_target_id,\n                MAX(CASE WHEN e.release = {sql_literal(pair['old_release'])} THEN e.score END) AS endpoint_score_old,\n                MAX(CASE WHEN e.release = {sql_literal(pair['new_release'])} THEN e.score END) AS endpoint_score_new\n            FROM selected_endpoint_definitions d\n            INNER JOIN selected_evidence_rows e\n              ON e.endpoint_id = d.endpoint_id\n             AND e.endpoint_role = d.endpoint_role\n            GROUP BY d.endpoint_level, d.endpoint_id, d.endpoint_role, e.canonical_disease_id, e.canonical_target_id\n            "
        )
    connection.execute(
        f"\n        CREATE TABLE selected_endpoint_definitions(endpoint_level VARCHAR, endpoint_id VARCHAR, endpoint_role VARCHAR);\n        INSERT INTO selected_endpoint_definitions VALUES {endpoint_values};\n        CREATE TABLE selected_evidence_rows AS\n        SELECT\n            CASE\n                WHEN endpoint_id = 'genetic_association' THEN 'selected_evidence_domain'\n                WHEN endpoint_id = 'gwas_credible_sets' THEN 'selected_source_decomposition'\n                ELSE 'unselected'\n            END AS endpoint_role,\n            LOWER(CAST(endpoint_level AS VARCHAR)) AS endpoint_level,\n            CAST(endpoint_id AS VARCHAR) AS endpoint_id,\n            CAST(release AS VARCHAR) AS release,\n            canonical_disease_id,\n            canonical_target_id,\n            CAST(score AS DOUBLE) AS score,\n            UPPER(CAST(association_state AS VARCHAR)) AS association_state,\n            CAST(domain_support_eligible AS BOOLEAN) AS domain_support_eligible\n        FROM support_construction_evidence_domain\n        WHERE ({endpoint_filter})\n          AND CAST(domain_support_eligible AS BOOLEAN);\n        CREATE TABLE endpoint_pair_scores AS\n        {' UNION ALL '.join(pair_queries)};\n        CREATE TABLE evidence_domain_matches AS\n        SELECT\n            p.release_pair,\n            p.support_stratum,\n            p.canonical_disease_id,\n            p.canonical_target_id,\n            p.absolute_change,\n            s.endpoint_level,\n            s.endpoint_id,\n            s.endpoint_role,\n            s.endpoint_score_old,\n            s.endpoint_score_new,\n            s.endpoint_score_new - s.endpoint_score_old AS endpoint_signed_change,\n            ABS(s.endpoint_score_new - s.endpoint_score_old) AS endpoint_absolute_change\n        FROM pair_score_revisions p\n        INNER JOIN endpoint_pair_scores s\n          ON s.release_pair = p.release_pair\n         AND s.canonical_disease_id = p.canonical_disease_id\n         AND s.canonical_target_id = p.canonical_target_id\n        WHERE s.endpoint_score_old IS NOT NULL\n          AND s.endpoint_score_new IS NOT NULL\n        "
    )
    threshold = float(
        configuration["evidence_domain_context"]["large_change_threshold"]
    )
    interpretation = sql_literal(
        configuration["evidence_domain_context"]["interpretation"]
    )
    connection.execute(
        f"\n        CREATE TABLE evidence_domain_tail_context AS\n        WITH overall_totals AS (\n            SELECT\n                release_pair,\n                support_stratum,\n                COUNT(*) AS all_pair_count,\n                SUM(CAST(absolute_change >= {threshold} AS BIGINT)) AS large_change_pair_count,\n                SUM(CAST(absolute_change < {threshold} AS BIGINT)) AS non_large_change_pair_count\n            FROM pair_score_revisions\n            GROUP BY release_pair, support_stratum\n        ), endpoint_totals AS (\n            SELECT\n                release_pair,\n                support_stratum,\n                endpoint_level,\n                endpoint_id,\n                endpoint_role,\n                COUNT(*) AS endpoint_persistent_pair_count,\n                SUM(CAST(absolute_change >= {threshold} AS BIGINT)) AS endpoint_persistent_large_change_count,\n                SUM(CAST(absolute_change < {threshold} AS BIGINT)) AS endpoint_persistent_non_large_change_count,\n                QUANTILE_CONT(CASE WHEN absolute_change >= {threshold} THEN endpoint_absolute_change END, 0.5) AS endpoint_change_median_in_overall_large_tail,\n                QUANTILE_CONT(CASE WHEN absolute_change >= {threshold} THEN endpoint_absolute_change END, 0.9) AS endpoint_change_q90_in_overall_large_tail\n            FROM evidence_domain_matches\n            GROUP BY release_pair, support_stratum, endpoint_level, endpoint_id, endpoint_role\n        ), groups AS (\n            SELECT\n                o.*,\n                d.endpoint_level,\n                d.endpoint_id,\n                d.endpoint_role\n            FROM overall_totals o\n            CROSS JOIN selected_endpoint_definitions d\n        ), disease_counts AS (\n            SELECT\n                release_pair,\n                support_stratum,\n                endpoint_id,\n                canonical_disease_id,\n                COUNT(*) AS contribution_count\n            FROM evidence_domain_matches\n            WHERE absolute_change >= {threshold}\n            GROUP BY release_pair, support_stratum, endpoint_id, canonical_disease_id\n        ), dominant_disease AS (\n            SELECT * EXCLUDE (order_index)\n            FROM (\n                SELECT *, ROW_NUMBER() OVER (\n                    PARTITION BY release_pair, support_stratum, endpoint_id\n                    ORDER BY contribution_count DESC, canonical_disease_id\n                ) AS order_index\n                FROM disease_counts\n            )\n            WHERE order_index = 1\n        ), target_counts AS (\n            SELECT\n                release_pair,\n                support_stratum,\n                endpoint_id,\n                canonical_target_id,\n                COUNT(*) AS contribution_count\n            FROM evidence_domain_matches\n            WHERE absolute_change >= {threshold}\n            GROUP BY release_pair, support_stratum, endpoint_id, canonical_target_id\n        ), dominant_target AS (\n            SELECT * EXCLUDE (order_index)\n            FROM (\n                SELECT *, ROW_NUMBER() OVER (\n                    PARTITION BY release_pair, support_stratum, endpoint_id\n                    ORDER BY contribution_count DESC, canonical_target_id\n                ) AS order_index\n                FROM target_counts\n            )\n            WHERE order_index = 1\n        ), area_counts AS (\n            SELECT\n                e.release_pair,\n                e.support_stratum,\n                e.endpoint_id,\n                a.therapeutic_area_id,\n                SUM(a.fractional_weight) AS contribution_weight\n            FROM evidence_domain_matches e\n            INNER JOIN stable_areas a USING (canonical_disease_id)\n            WHERE e.absolute_change >= {threshold}\n            GROUP BY e.release_pair, e.support_stratum, e.endpoint_id, a.therapeutic_area_id\n        ), area_totals AS (\n            SELECT release_pair, support_stratum, endpoint_id, SUM(contribution_weight) AS total_weight\n            FROM area_counts\n            GROUP BY release_pair, support_stratum, endpoint_id\n        ), dominant_area AS (\n            SELECT * EXCLUDE (order_index)\n            FROM (\n                SELECT *, ROW_NUMBER() OVER (\n                    PARTITION BY release_pair, support_stratum, endpoint_id\n                    ORDER BY contribution_weight DESC, therapeutic_area_id\n                ) AS order_index\n                FROM area_counts\n            )\n            WHERE order_index = 1\n        )\n        SELECT\n            g.release_pair,\n            g.support_stratum,\n            g.endpoint_level,\n            g.endpoint_id,\n            g.endpoint_role,\n            {threshold} AS large_change_threshold,\n            g.all_pair_count,\n            g.large_change_pair_count,\n            g.non_large_change_pair_count,\n            COALESCE(t.endpoint_persistent_pair_count, 0) AS endpoint_persistent_pair_count,\n            COALESCE(t.endpoint_persistent_large_change_count, 0) AS endpoint_persistent_large_change_count,\n            COALESCE(t.endpoint_persistent_non_large_change_count, 0) AS endpoint_persistent_non_large_change_count,\n            COALESCE(t.endpoint_persistent_large_change_count, 0) / NULLIF(CAST(g.large_change_pair_count AS DOUBLE), 0) AS endpoint_prevalence_in_large_tail,\n            COALESCE(t.endpoint_persistent_non_large_change_count, 0) / NULLIF(CAST(g.non_large_change_pair_count AS DOUBLE), 0) AS endpoint_prevalence_outside_large_tail,\n            endpoint_prevalence_in_large_tail - endpoint_prevalence_outside_large_tail AS endpoint_prevalence_difference,\n            endpoint_prevalence_in_large_tail / NULLIF(endpoint_prevalence_outside_large_tail, 0) AS endpoint_prevalence_ratio,\n            t.endpoint_change_median_in_overall_large_tail,\n            t.endpoint_change_q90_in_overall_large_tail,\n            d.canonical_disease_id AS dominant_disease_id,\n            d.contribution_count AS dominant_disease_count,\n            d.contribution_count / NULLIF(CAST(t.endpoint_persistent_large_change_count AS DOUBLE), 0) AS maximum_disease_share,\n            q.canonical_target_id AS dominant_target_id,\n            q.contribution_count AS dominant_target_count,\n            q.contribution_count / NULLIF(CAST(t.endpoint_persistent_large_change_count AS DOUBLE), 0) AS maximum_target_share,\n            a.therapeutic_area_id AS dominant_therapeutic_area_id,\n            a.contribution_weight AS dominant_therapeutic_area_weight,\n            a.contribution_weight / NULLIF(ar.total_weight, 0) AS maximum_therapeutic_area_share,\n            {interpretation} AS interpretation\n        FROM groups g\n        LEFT JOIN endpoint_totals t USING (release_pair, support_stratum, endpoint_level, endpoint_id, endpoint_role)\n        LEFT JOIN dominant_disease d USING (release_pair, support_stratum, endpoint_id)\n        LEFT JOIN dominant_target q USING (release_pair, support_stratum, endpoint_id)\n        LEFT JOIN dominant_area a USING (release_pair, support_stratum, endpoint_id)\n        LEFT JOIN area_totals ar USING (release_pair, support_stratum, endpoint_id)\n        "
    )

def threshold_summary_sql(
    cut_points: Sequence[float], expression: str = "absolute_change"
) -> str:
    fields: list[str] = []
    for (cut_point) in (cut_points):
        suffix = f"{cut_point:.2f}".replace(".", "_")
        fields.append(
            f"SUM(CAST({expression} >= {float(cut_point)} AS BIGINT)) / CAST(COUNT(*) AS DOUBLE) AS absolute_change_at_least_{suffix}_fraction"
        )
    return ",\n".join(fields)

def create_score_summary_table(
    connection: Any, configuration: Mapping[str, Any]
) -> None:
    tolerance = float(configuration["scores"]["descriptive_tolerance"])
    threshold_fields = threshold_summary_sql(
        configuration["scores"]["absolute_change_cut_points"]
    )
    status_values = ",".join(
        (
            f"({sql_literal(value)})"
            for (value) in (
                configuration["implementation_scale_diagnostic"][
                    "evidence_count_statuses"
                ]
            )
        )
    )
    mapping_values = ",".join(
        (f"({sql_literal(value)})" for (value) in (configuration["mapping_contexts"]))
    )
    release_pair_values = ",".join(
        (
            f"({sql_literal(pair['label'])})"
            for (pair) in (configuration["release_pairs"])
        )
    )
    support_values = ",".join(
        (f"({sql_literal(value)})" for (value) in (configuration["supports"].values()))
    )
    connection.execute(
        f"\n        CREATE TABLE score_summary AS\n        WITH expanded AS (\n            SELECT\n                *,\n                'OVERALL' AS summary_scope,\n                'all_mapping_contexts' AS mapping_context_group,\n                'all_evidence_count_statuses' AS evidence_count_status_group\n            FROM pair_score_revisions\n            UNION ALL\n            SELECT\n                *,\n                'MAPPING_CONTEXT' AS summary_scope,\n                mapping_context AS mapping_context_group,\n                'all_evidence_count_statuses' AS evidence_count_status_group\n            FROM pair_score_revisions\n            UNION ALL\n            SELECT\n                *,\n                'EVIDENCE_COUNT_STATUS' AS summary_scope,\n                'all_mapping_contexts' AS mapping_context_group,\n                evidence_count_status AS evidence_count_status_group\n            FROM pair_score_revisions\n        )\n        SELECT\n            release_pair,\n            support_stratum,\n            summary_scope,\n            mapping_context_group,\n            evidence_count_status_group,\n            COUNT(*) AS record_count,\n            COUNT(DISTINCT canonical_disease_id) AS disease_count,\n            COUNT(DISTINCT canonical_target_id) AS target_count,\n            SUM(CAST(exact_unchanged AS BIGINT)) / CAST(COUNT(*) AS DOUBLE) AS exact_unchanged_fraction,\n            SUM(CAST(absolute_change <= {tolerance} AS BIGINT)) / CAST(COUNT(*) AS DOUBLE) AS near_unchanged_fraction,\n            QUANTILE_CONT(signed_change, 0.5) AS signed_median,\n            QUANTILE_CONT(absolute_change, 0.5) AS absolute_median,\n            QUANTILE_CONT(absolute_change, 0.75) AS absolute_q75,\n            QUANTILE_CONT(absolute_change, 0.9) AS absolute_q90,\n            QUANTILE_CONT(absolute_change, 0.95) AS absolute_q95,\n            QUANTILE_CONT(absolute_change, 0.99) AS absolute_q99,\n            MAX(absolute_change) AS maximum_absolute_change,\n            {threshold_fields},\n            QUANTILE_CONT(score_old, 0.25) AS old_score_q25,\n            QUANTILE_CONT(score_old, 0.5) AS old_score_median,\n            QUANTILE_CONT(score_old, 0.75) AS old_score_q75,\n            QUANTILE_CONT(score_new, 0.25) AS new_score_q25,\n            QUANTILE_CONT(score_new, 0.5) AS new_score_median,\n            QUANTILE_CONT(score_new, 0.75) AS new_score_q75,\n            SUM(CAST(signed_change > 0 AS BIGINT)) AS positive_change_count,\n            SUM(CAST(signed_change < 0 AS BIGINT)) AS negative_change_count,\n            SUM(CAST(signed_change = 0 AS BIGINT)) AS zero_change_count\n        FROM expanded\n        GROUP BY release_pair, support_stratum, summary_scope, mapping_context_group, evidence_count_status_group\n        "
    )
    connection.execute(
        f"\n        INSERT INTO score_summary (\n            release_pair,\n            support_stratum,\n            summary_scope,\n            mapping_context_group,\n            evidence_count_status_group,\n            record_count,\n            disease_count,\n            target_count,\n            positive_change_count,\n            negative_change_count,\n            zero_change_count\n        )\n        SELECT\n            g.release_pair,\n            g.support_stratum,\n            'EVIDENCE_COUNT_STATUS',\n            'all_mapping_contexts',\n            g.evidence_count_status_group,\n            0,\n            0,\n            0,\n            0,\n            0,\n            0\n        FROM (\n            SELECT DISTINCT p.release_pair, p.support_stratum, s.evidence_count_status_group\n            FROM pair_score_revisions p\n            CROSS JOIN (VALUES {status_values}) s(evidence_count_status_group)\n        ) g\n        LEFT JOIN score_summary x\n          ON x.release_pair = g.release_pair\n         AND x.support_stratum = g.support_stratum\n         AND x.summary_scope = 'EVIDENCE_COUNT_STATUS'\n         AND x.evidence_count_status_group = g.evidence_count_status_group\n        WHERE x.release_pair IS NULL\n        "
    )
    connection.execute(
        f"\n        INSERT INTO score_summary (\n            release_pair,\n            support_stratum,\n            summary_scope,\n            mapping_context_group,\n            evidence_count_status_group,\n            record_count,\n            disease_count,\n            target_count,\n            positive_change_count,\n            negative_change_count,\n            zero_change_count\n        )\n        SELECT\n            g.release_pair,\n            g.support_stratum,\n            'MAPPING_CONTEXT',\n            g.mapping_context_group,\n            'all_evidence_count_statuses',\n            0,\n            0,\n            0,\n            0,\n            0,\n            0\n        FROM (\n            SELECT r.release_pair, s.support_stratum, m.mapping_context_group\n            FROM (VALUES {release_pair_values}) r(release_pair)\n            CROSS JOIN (VALUES {support_values}) s(support_stratum)\n            CROSS JOIN (VALUES {mapping_values}) m(mapping_context_group)\n        ) g\n        LEFT JOIN score_summary x\n          ON x.release_pair = g.release_pair\n         AND x.support_stratum = g.support_stratum\n         AND x.summary_scope = 'MAPPING_CONTEXT'\n         AND x.mapping_context_group = g.mapping_context_group\n        WHERE x.release_pair IS NULL\n        "
    )

def create_implementation_scale_diagnostic_table(
    connection: Any, configuration: Mapping[str, Any]
) -> None:
    diagnostic = configuration["implementation_scale_diagnostic"]
    release_pair = sql_literal(diagnostic["release_pair"])
    threshold_fields_raw = threshold_summary_sql(
        configuration["scores"]["absolute_change_cut_points"], "absolute_change"
    )
    threshold_fields_adjusted = threshold_summary_sql(
        configuration["scores"]["absolute_change_cut_points"],
        "implementation_adjusted_absolute_change",
    ).replace(" AS absolute_change_", " AS adjusted_absolute_change_")
    status_values = ",".join(
        (
            f"({sql_literal(value)})"
            for (value) in (diagnostic["evidence_count_statuses"])
        )
    )
    connection.execute(
        f"\n        CREATE TABLE implementation_scale_diagnostic AS\n        WITH selected AS (\n            SELECT * FROM pair_score_revisions\n            WHERE release_pair = {release_pair}\n              AND calibration_applicable\n        ), expanded AS (\n            SELECT\n                *,\n                'OVERALL' AS summary_scope,\n                'all_mapping_contexts' AS mapping_context_group,\n                'all_evidence_count_statuses' AS evidence_count_status_group\n            FROM selected\n            UNION ALL\n            SELECT\n                *,\n                'MAPPING_CONTEXT' AS summary_scope,\n                mapping_context AS mapping_context_group,\n                'all_evidence_count_statuses' AS evidence_count_status_group\n            FROM selected\n            UNION ALL\n            SELECT\n                *,\n                'EVIDENCE_COUNT_STATUS' AS summary_scope,\n                'all_mapping_contexts' AS mapping_context_group,\n                evidence_count_status AS evidence_count_status_group\n            FROM selected\n        )\n        SELECT\n            release_pair,\n            support_stratum,\n            summary_scope,\n            mapping_context_group,\n            evidence_count_status_group,\n            MAX(calibration_factor) AS old_score_multiplier,\n            {sql_literal(diagnostic['calibrated_old_score_rule'])} AS calibrated_old_score_rule,\n            {sql_literal(diagnostic['interpretation'])} AS interpretation,\n            COUNT(*) AS record_count,\n            QUANTILE_CONT(signed_change, 0.5) AS raw_signed_median,\n            QUANTILE_CONT(absolute_change, 0.5) AS raw_absolute_median,\n            QUANTILE_CONT(absolute_change, 0.9) AS raw_absolute_q90,\n            QUANTILE_CONT(absolute_change, 0.99) AS raw_absolute_q99,\n            QUANTILE_CONT(implementation_adjusted_signed_change, 0.5) AS adjusted_signed_median,\n            QUANTILE_CONT(implementation_adjusted_absolute_change, 0.5) AS adjusted_absolute_median,\n            QUANTILE_CONT(implementation_adjusted_absolute_change, 0.9) AS adjusted_absolute_q90,\n            QUANTILE_CONT(implementation_adjusted_absolute_change, 0.99) AS adjusted_absolute_q99,\n            MAX(calibrated_score_old) AS maximum_calibrated_old_score,\n            {threshold_fields_raw},\n            {threshold_fields_adjusted}\n        FROM expanded\n        GROUP BY release_pair, support_stratum, summary_scope, mapping_context_group, evidence_count_status_group\n        "
    )
    connection.execute(
        f"\n        INSERT INTO implementation_scale_diagnostic (\n            release_pair,\n            support_stratum,\n            summary_scope,\n            mapping_context_group,\n            evidence_count_status_group,\n            old_score_multiplier,\n            calibrated_old_score_rule,\n            interpretation,\n            record_count\n        )\n        SELECT\n            g.release_pair,\n            g.support_stratum,\n            'EVIDENCE_COUNT_STATUS',\n            'all_mapping_contexts',\n            g.evidence_count_status_group,\n            {float(diagnostic['old_score_multiplier'])},\n            {sql_literal(diagnostic['calibrated_old_score_rule'])},\n            {sql_literal(diagnostic['interpretation'])},\n            0\n        FROM (\n            SELECT DISTINCT p.release_pair, p.support_stratum, s.evidence_count_status_group\n            FROM pair_score_revisions p\n            CROSS JOIN (VALUES {status_values}) s(evidence_count_status_group)\n            WHERE p.release_pair = {release_pair}\n        ) g\n        LEFT JOIN implementation_scale_diagnostic x\n          ON x.release_pair = g.release_pair\n         AND x.support_stratum = g.support_stratum\n         AND x.summary_scope = 'EVIDENCE_COUNT_STATUS'\n         AND x.evidence_count_status_group = g.evidence_count_status_group\n        WHERE x.release_pair IS NULL\n        "
    )

def create_score_tier_summary_table(
    connection: Any, configuration: Mapping[str, Any]
) -> None:
    context = configuration["entrant_exit_context"]
    threshold_fields = threshold_summary_sql(
        configuration["scores"]["absolute_change_cut_points"]
    )
    connection.execute(
        f"\n        CREATE TABLE score_tier_summary AS\n        WITH transition_counts AS (\n            SELECT\n                release_pair,\n                score_tier,\n                SUM(CAST(context_class = {sql_literal(context['entrant_class'])} AS BIGINT)) AS entrant_count,\n                SUM(CAST(context_class = {sql_literal(context['exit_class'])} AS BIGINT)) AS exit_count,\n                SUM(CAST(context_class = {sql_literal(context['mapping_affected_class'])} AS BIGINT)) AS mapping_affected_count\n            FROM entrant_exit_context\n            GROUP BY release_pair, score_tier\n        ), persistent AS (\n            SELECT\n                release_pair,\n                support_stratum,\n                old_score_tier,\n                COUNT(*) AS pair_count,\n                COUNT(DISTINCT canonical_disease_id) AS disease_count,\n                COUNT(DISTINCT canonical_target_id) AS target_count,\n                QUANTILE_CONT(absolute_change, 0.5) AS median_absolute_change,\n                QUANTILE_CONT(absolute_change, 0.9) AS q90_absolute_change,\n                QUANTILE_CONT(signed_change, 0.5) AS median_signed_change,\n                SUM(CAST(signed_change > 0 AS BIGINT)) AS positive_change_count,\n                SUM(CAST(signed_change < 0 AS BIGINT)) AS negative_change_count,\n                {threshold_fields}\n            FROM pair_score_revisions\n            GROUP BY release_pair, support_stratum, old_score_tier\n        )\n        SELECT\n            p.*,\n            COALESCE(t.entrant_count, 0) AS release_native_entrant_count_in_tier,\n            COALESCE(t.exit_count, 0) AS release_native_exit_count_in_tier,\n            COALESCE(t.mapping_affected_count, 0) AS mapping_affected_transition_count_in_tier\n        FROM persistent p\n        LEFT JOIN transition_counts t\n          ON t.release_pair = p.release_pair\n         AND t.score_tier = p.old_score_tier\n        "
    )

def create_score_tier_area_distribution_table(connection: Any) -> None:
    connection.execute(
        "\n        CREATE TABLE score_tier_area_distribution AS\n        WITH assignments AS (\n            SELECT\n                p.release_pair,\n                p.support_stratum,\n                p.old_score_tier,\n                p.canonical_disease_id,\n                p.canonical_target_id,\n                a.therapeutic_area_id,\n                a.therapeutic_area_label,\n                CAST(a.fractional_weight AS DECIMAL(38,18)) AS fractional_weight\n            FROM pair_score_revisions p\n            INNER JOIN stable_areas a USING (canonical_disease_id)\n        ), totals AS (\n            SELECT\n                release_pair,\n                support_stratum,\n                old_score_tier,\n                SUM(fractional_weight) AS tier_fractional_pair_count\n            FROM assignments\n            GROUP BY release_pair, support_stratum, old_score_tier\n        )\n        SELECT\n            a.release_pair,\n            a.support_stratum,\n            a.old_score_tier,\n            a.therapeutic_area_id,\n            MAX(a.therapeutic_area_label) AS therapeutic_area_label,\n            COUNT(*) AS pair_area_assignment_count,\n            SUM(a.fractional_weight) AS fractional_pair_count,\n            MAX(t.tier_fractional_pair_count) AS tier_fractional_pair_count,\n            COUNT(DISTINCT a.canonical_disease_id) AS disease_count,\n            COUNT(DISTINCT a.canonical_target_id) AS target_count,\n            SUM(a.fractional_weight) / NULLIF(MAX(t.tier_fractional_pair_count), 0) AS within_tier_fraction\n        FROM assignments a\n        INNER JOIN totals t USING (release_pair, support_stratum, old_score_tier)\n        GROUP BY a.release_pair, a.support_stratum, a.old_score_tier, a.therapeutic_area_id\n        "
    )

def create_therapeutic_area_summary_table(
    connection: Any, configuration: Mapping[str, Any]
) -> None:
    large_threshold = float(configuration["scores"]["large_change_threshold"])
    minimum_diseases = int(
        configuration["therapeutic_areas"]["minimum_unique_diseases_for_comparison"]
    )
    connection.execute(
        f"\n        CREATE TABLE therapeutic_area_summary AS\n        WITH expanded AS (\n            SELECT\n                p.*,\n                a.therapeutic_area_id,\n                a.therapeutic_area_label,\n                a.fractional_weight\n            FROM pair_score_revisions p\n            INNER JOIN stable_areas a USING (canonical_disease_id)\n        ), ordered AS (\n            SELECT\n                *,\n                SUM(fractional_weight) OVER (\n                    PARTITION BY release_pair, support_stratum, therapeutic_area_id\n                    ORDER BY absolute_change, canonical_disease_id, canonical_target_id\n                    ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW\n                ) AS cumulative_weight,\n                SUM(fractional_weight) OVER (\n                    PARTITION BY release_pair, support_stratum, therapeutic_area_id\n                ) AS total_weight\n            FROM expanded\n        ), base AS (\n            SELECT\n                release_pair,\n                support_stratum,\n                therapeutic_area_id,\n                MAX(therapeutic_area_label) AS therapeutic_area_label,\n                COUNT(*) AS pair_area_assignment_count,\n                COUNT(DISTINCT canonical_disease_id || '|' || canonical_target_id) AS unique_pair_count,\n                SUM(fractional_weight) AS fractional_pair_count,\n                COUNT(DISTINCT canonical_disease_id) AS disease_count,\n                COUNT(DISTINCT canonical_target_id) AS target_count,\n                MIN(CASE WHEN cumulative_weight >= total_weight * 0.5 THEN absolute_change END) AS weighted_median_absolute_change,\n                MIN(CASE WHEN cumulative_weight >= total_weight * 0.9 THEN absolute_change END) AS weighted_q90_absolute_change,\n                SUM(CASE WHEN absolute_change >= {large_threshold} THEN fractional_weight ELSE CAST(0 AS DECIMAL(38,18)) END) AS large_change_fractional_pair_count,\n                SUM(CASE WHEN absolute_change >= {large_threshold} THEN fractional_weight ELSE CAST(0 AS DECIMAL(38,18)) END) / NULLIF(SUM(fractional_weight), 0) AS weighted_large_change_fraction,\n                COUNT(DISTINCT CASE WHEN absolute_change >= {large_threshold} THEN canonical_disease_id END) AS large_change_disease_count,\n                COUNT(DISTINCT CASE WHEN absolute_change >= {large_threshold} THEN canonical_target_id END) AS large_change_target_count,\n                COUNT(DISTINCT canonical_disease_id) >= {minimum_diseases} AS primary_area_comparison_eligible,\n                CASE\n                    WHEN COUNT(DISTINCT canonical_disease_id) >= {minimum_diseases} THEN 'comparison_eligible'\n                    ELSE 'descriptive_only'\n                END AS reporting_status\n            FROM ordered\n            GROUP BY release_pair, support_stratum, therapeutic_area_id\n        ), disease_counts AS (\n            SELECT\n                release_pair,\n                support_stratum,\n                therapeutic_area_id,\n                canonical_disease_id,\n                SUM(fractional_weight) AS contribution_weight\n            FROM expanded\n            WHERE absolute_change >= {large_threshold}\n            GROUP BY release_pair, support_stratum, therapeutic_area_id, canonical_disease_id\n        ), dominant_disease AS (\n            SELECT * EXCLUDE (order_index)\n            FROM (\n                SELECT *, ROW_NUMBER() OVER (\n                    PARTITION BY release_pair, support_stratum, therapeutic_area_id\n                    ORDER BY contribution_weight DESC, canonical_disease_id\n                ) AS order_index\n                FROM disease_counts\n            )\n            WHERE order_index = 1\n        ), target_counts AS (\n            SELECT\n                release_pair,\n                support_stratum,\n                therapeutic_area_id,\n                canonical_target_id,\n                SUM(fractional_weight) AS contribution_weight\n            FROM expanded\n            WHERE absolute_change >= {large_threshold}\n            GROUP BY release_pair, support_stratum, therapeutic_area_id, canonical_target_id\n        ), dominant_target AS (\n            SELECT * EXCLUDE (order_index)\n            FROM (\n                SELECT *, ROW_NUMBER() OVER (\n                    PARTITION BY release_pair, support_stratum, therapeutic_area_id\n                    ORDER BY contribution_weight DESC, canonical_target_id\n                ) AS order_index\n                FROM target_counts\n            )\n            WHERE order_index = 1\n        )\n        SELECT\n            b.*,\n            d.canonical_disease_id AS dominant_large_change_disease_id,\n            d.contribution_weight AS dominant_large_change_disease_weight,\n            d.contribution_weight / NULLIF(b.large_change_fractional_pair_count, 0) AS maximum_disease_share,\n            t.canonical_target_id AS dominant_large_change_target_id,\n            t.contribution_weight AS dominant_large_change_target_weight,\n            t.contribution_weight / NULLIF(b.large_change_fractional_pair_count, 0) AS maximum_target_share\n        FROM base b\n        LEFT JOIN dominant_disease d USING (release_pair, support_stratum, therapeutic_area_id)\n        LEFT JOIN dominant_target t USING (release_pair, support_stratum, therapeutic_area_id)\n        "
    )

def create_concentration_tables(
    connection: Any, configuration: Mapping[str, Any]
) -> None:
    threshold = float(configuration["concentration"]["large_change_threshold"])
    connection.execute(
        f"\n        CREATE TABLE concentration_groups AS\n        SELECT\n            release_pair,\n            support_stratum,\n            COUNT(*) AS all_pair_count,\n            QUANTILE_CONT(absolute_change, 0.9) AS full_q90_absolute_change,\n            SUM(CAST(absolute_change >= {threshold} AS BIGINT)) AS large_change_pair_count,\n            COUNT(DISTINCT CASE WHEN absolute_change >= {threshold} THEN canonical_disease_id END) AS large_change_disease_count,\n            COUNT(DISTINCT CASE WHEN absolute_change >= {threshold} THEN canonical_target_id END) AS large_change_target_count\n        FROM pair_score_revisions\n        GROUP BY release_pair, support_stratum\n        "
    )
    connection.execute(
        f"\n        CREATE TABLE dominant_disease AS\n        WITH counts AS (\n            SELECT release_pair, support_stratum, canonical_disease_id, COUNT(*) AS contribution_count\n            FROM pair_score_revisions\n            WHERE absolute_change >= {threshold}\n            GROUP BY release_pair, support_stratum, canonical_disease_id\n        )\n        SELECT * EXCLUDE (order_index)\n        FROM (\n            SELECT *, ROW_NUMBER() OVER (\n                PARTITION BY release_pair, support_stratum\n                ORDER BY contribution_count DESC, canonical_disease_id\n            ) AS order_index\n            FROM counts\n        )\n        WHERE order_index = 1\n        "
    )
    connection.execute(
        f"\n        CREATE TABLE dominant_target AS\n        WITH counts AS (\n            SELECT release_pair, support_stratum, canonical_target_id, COUNT(*) AS contribution_count\n            FROM pair_score_revisions\n            WHERE absolute_change >= {threshold}\n            GROUP BY release_pair, support_stratum, canonical_target_id\n        )\n        SELECT * EXCLUDE (order_index)\n        FROM (\n            SELECT *, ROW_NUMBER() OVER (\n                PARTITION BY release_pair, support_stratum\n                ORDER BY contribution_count DESC, canonical_target_id\n            ) AS order_index\n            FROM counts\n        )\n        WHERE order_index = 1\n        "
    )
    connection.execute(
        f"\n        CREATE TABLE dominant_area AS\n        WITH counts AS (\n            SELECT\n                p.release_pair,\n                p.support_stratum,\n                a.therapeutic_area_id,\n                SUM(a.fractional_weight) AS contribution_weight\n            FROM pair_score_revisions p\n            INNER JOIN stable_areas a USING (canonical_disease_id)\n            WHERE p.absolute_change >= {threshold}\n            GROUP BY p.release_pair, p.support_stratum, a.therapeutic_area_id\n        )\n        SELECT * EXCLUDE (order_index)\n        FROM (\n            SELECT *, ROW_NUMBER() OVER (\n                PARTITION BY release_pair, support_stratum\n                ORDER BY contribution_weight DESC, therapeutic_area_id\n            ) AS order_index\n            FROM counts\n        )\n        WHERE order_index = 1\n        "
    )
    connection.execute(
        f"\n        CREATE TABLE concentration_diagnostics AS\n        WITH area_totals AS (\n            SELECT\n                p.release_pair,\n                p.support_stratum,\n                SUM(a.fractional_weight) AS stable_area_large_change_weight\n            FROM pair_score_revisions p\n            INNER JOIN stable_areas a USING (canonical_disease_id)\n            WHERE p.absolute_change >= {threshold}\n            GROUP BY p.release_pair, p.support_stratum\n        ), leave_disease AS (\n            SELECT\n                p.release_pair,\n                p.support_stratum,\n                COUNT(*) AS remaining_pair_count,\n                QUANTILE_CONT(p.absolute_change, 0.9) AS remaining_q90\n            FROM pair_score_revisions p\n            INNER JOIN dominant_disease d USING (release_pair, support_stratum)\n            WHERE p.canonical_disease_id <> d.canonical_disease_id\n            GROUP BY p.release_pair, p.support_stratum\n        ), leave_target AS (\n            SELECT\n                p.release_pair,\n                p.support_stratum,\n                COUNT(*) AS remaining_pair_count,\n                QUANTILE_CONT(p.absolute_change, 0.9) AS remaining_q90\n            FROM pair_score_revisions p\n            INNER JOIN dominant_target t USING (release_pair, support_stratum)\n            WHERE p.canonical_target_id <> t.canonical_target_id\n            GROUP BY p.release_pair, p.support_stratum\n        )\n        SELECT\n            g.release_pair,\n            g.support_stratum,\n            {threshold} AS large_change_threshold,\n            g.all_pair_count,\n            g.large_change_pair_count,\n            g.large_change_disease_count,\n            g.large_change_target_count,\n            d.canonical_disease_id AS dominant_disease_id,\n            d.contribution_count AS dominant_disease_large_change_count,\n            d.contribution_count / NULLIF(CAST(g.large_change_pair_count AS DOUBLE), 0) AS maximum_disease_share,\n            t.canonical_target_id AS dominant_target_id,\n            t.contribution_count AS dominant_target_large_change_count,\n            t.contribution_count / NULLIF(CAST(g.large_change_pair_count AS DOUBLE), 0) AS maximum_target_share,\n            a.therapeutic_area_id AS dominant_therapeutic_area_id,\n            a.contribution_weight AS dominant_therapeutic_area_weight,\n            a.contribution_weight / NULLIF(ar.stable_area_large_change_weight, 0) AS maximum_therapeutic_area_share,\n            ar.stable_area_large_change_weight,\n            g.full_q90_absolute_change,\n            ld.remaining_pair_count AS leave_dominant_disease_out_pair_count,\n            ld.remaining_q90 AS leave_dominant_disease_out_q90_absolute_change,\n            ld.remaining_q90 - g.full_q90_absolute_change AS leave_dominant_disease_out_q90_difference,\n            lt.remaining_pair_count AS leave_dominant_target_out_pair_count,\n            lt.remaining_q90 AS leave_dominant_target_out_q90_absolute_change,\n            lt.remaining_q90 - g.full_q90_absolute_change AS leave_dominant_target_out_q90_difference\n        FROM concentration_groups g\n        LEFT JOIN dominant_disease d USING (release_pair, support_stratum)\n        LEFT JOIN dominant_target t USING (release_pair, support_stratum)\n        LEFT JOIN dominant_area a USING (release_pair, support_stratum)\n        LEFT JOIN area_totals ar USING (release_pair, support_stratum)\n        LEFT JOIN leave_disease ld USING (release_pair, support_stratum)\n        LEFT JOIN leave_target lt USING (release_pair, support_stratum)\n        "
    )

def create_trajectory_table(connection: Any, configuration: Mapping[str, Any]) -> None:
    tolerance = float(configuration["trajectories"]["stable_tolerance"])
    tiers = configuration["score_tiers"]
    tier_case = score_tier_case("score_25_12", tiers)
    category_values = ",".join(
        (
            f"({sql_literal(value)})"
            for (value) in (configuration["trajectories"]["categories"])
        )
    )
    tier_group_values = ",".join(
        (
            f"({sql_literal(value)})"
            for (value) in (
                [*[str(tier["label"]) for (tier) in (tiers)], "all_score_tiers"]
            )
        )
    )
    connection.execute(
        f"\n        CREATE TABLE trajectory_rows AS\n        WITH changes AS (\n            SELECT\n                canonical_disease_id,\n                canonical_target_id,\n                CAST(score_25_12 AS DOUBLE) AS score_25_12,\n                CAST(score_26_03 AS DOUBLE) AS score_26_03,\n                CAST(score_26_06 AS DOUBLE) AS score_26_06,\n                CAST(score_26_03 AS DOUBLE) - CAST(score_25_12 AS DOUBLE) AS first_change,\n                CAST(score_26_06 AS DOUBLE) - CAST(score_26_03 AS DOUBLE) AS second_change,\n                CAST(stable_entity_eligible AS BOOLEAN) AS stable_entity_eligible\n            FROM support_construction_three_release\n        ), states AS (\n            SELECT\n                *,\n                CASE WHEN ABS(first_change) <= {tolerance} THEN 'stable' WHEN first_change > 0 THEN 'increase' ELSE 'decrease' END AS first_state,\n                CASE WHEN ABS(second_change) <= {tolerance} THEN 'stable' WHEN second_change > 0 THEN 'increase' ELSE 'decrease' END AS second_state\n            FROM changes\n        )\n        SELECT\n            *,\n            CASE\n                WHEN first_state = 'stable' AND second_state = 'stable' THEN 'stable-stable'\n                WHEN first_state = 'stable' OR second_state = 'stable' THEN 'one-stable-one-changing'\n                ELSE first_state || '-' || second_state\n            END AS trajectory,\n            {tier_case} AS starting_score_tier\n        FROM states\n        "
    )
    connection.execute(
        f"\n        CREATE TABLE trajectory_summary AS\n        WITH expanded AS (\n            SELECT *, starting_score_tier AS score_tier_group FROM trajectory_rows\n            UNION ALL\n            SELECT *, 'all_score_tiers' AS score_tier_group FROM trajectory_rows\n        ), categories AS (\n            SELECT trajectory FROM (VALUES {category_values}) c(trajectory)\n        ), tier_groups AS (\n            SELECT score_tier_group FROM (VALUES {tier_group_values}) t(score_tier_group)\n        ), summary_grid AS (\n            SELECT c.trajectory, t.score_tier_group\n            FROM categories c\n            CROSS JOIN tier_groups t\n        ), aggregated AS (\n            SELECT\n                trajectory,\n                score_tier_group,\n                COUNT(*) AS pair_count,\n                COUNT(DISTINCT canonical_disease_id) AS disease_count,\n                COUNT(DISTINCT canonical_target_id) AS target_count,\n                QUANTILE_CONT(ABS(first_change), 0.5) AS first_change_absolute_median,\n                QUANTILE_CONT(ABS(second_change), 0.5) AS second_change_absolute_median\n            FROM expanded\n            GROUP BY trajectory, score_tier_group\n        ), totals AS (\n            SELECT score_tier_group, COUNT(*) AS tier_total\n            FROM expanded\n            GROUP BY score_tier_group\n        ), overall AS (\n            SELECT COUNT(*) AS overall_total FROM trajectory_rows\n        )\n        SELECT\n            '25.12_to_26.06' AS release_span,\n            g.trajectory,\n            g.score_tier_group,\n            COALESCE(a.pair_count, 0) AS pair_count,\n            COALESCE(a.disease_count, 0) AS disease_count,\n            COALESCE(a.target_count, 0) AS target_count,\n            CASE\n                WHEN COALESCE(t.tier_total, 0) = 0 THEN 0.0\n                ELSE COALESCE(a.pair_count, 0) / CAST(t.tier_total AS DOUBLE)\n            END AS within_tier_fraction,\n            CASE\n                WHEN o.overall_total = 0 THEN 0.0\n                ELSE COALESCE(a.pair_count, 0) / CAST(o.overall_total AS DOUBLE)\n            END AS overall_fraction,\n            a.first_change_absolute_median,\n            a.second_change_absolute_median\n        FROM summary_grid g\n        LEFT JOIN aggregated a USING (trajectory, score_tier_group)\n        LEFT JOIN totals t USING (score_tier_group)\n        CROSS JOIN overall o\n        "
    )

def create_entrant_exit_summary_table(
    connection: Any, configuration: Mapping[str, Any]
) -> None:
    missing_status = sql_literal(
        configuration["entrant_exit_context"]["missing_rank_context_status"]
    )
    context_values = ",".join(
        (
            f"({sql_literal(configuration['entrant_exit_context'][key])})"
            for (key) in (
                [
                    "entrant_class",
                    "exit_class",
                    "mapping_affected_class",
                    "present_zero_context_class",
                    "present_null_context_class",
                    "schema_incomparable_context_class",
                    "remaining_other_context_class",
                ]
            )
        )
    )
    connection.execute(
        f"\n        CREATE TABLE entrant_exit_summary AS\n        WITH overall AS (\n            SELECT\n                'overall' AS summary_scope,\n                release_pair,\n                context_class,\n                CAST(NULL AS VARCHAR) AS therapeutic_area_id,\n                COUNT(*) AS record_count,\n                CAST(COUNT(*) AS DOUBLE) AS fractional_record_count,\n                COUNT(DISTINCT canonical_disease_id) AS disease_count,\n                COUNT(DISTINCT canonical_target_id) AS target_count,\n                COUNT(applicable_score) AS defined_score_count,\n                QUANTILE_CONT(applicable_score, 0.5) AS applicable_score_median,\n                QUANTILE_CONT(applicable_score, 0.9) AS applicable_score_q90,\n                COUNT(old_score) AS old_defined_score_count,\n                QUANTILE_CONT(old_score, 0.5) AS old_score_median,\n                QUANTILE_CONT(old_score, 0.9) AS old_score_q90,\n                COUNT(new_score) AS new_defined_score_count,\n                QUANTILE_CONT(new_score, 0.5) AS new_score_median,\n                QUANTILE_CONT(new_score, 0.9) AS new_score_q90,\n                COUNT(applicable_native_rank) AS native_rank_context_count,\n                COUNT(*) - COUNT(applicable_native_rank) AS missing_native_rank_context_count,\n                CAST(SUM(CAST(applicable_top_rank_context AS BIGINT)) AS DOUBLE) AS top_rank_context_count,\n                SUM(CAST(applicable_top_rank_context AS BIGINT)) / NULLIF(CAST(COUNT(applicable_top_rank_context) AS DOUBLE), 0) AS top_rank_context_fraction,\n                CASE WHEN COUNT(*) - COUNT(applicable_native_rank) > 0 THEN {missing_status} ELSE 'complete' END AS rank_context_status\n            FROM entrant_exit_context\n            GROUP BY release_pair, context_class\n        ), area AS (\n            SELECT\n                'therapeutic_area' AS summary_scope,\n                e.release_pair,\n                e.context_class,\n                a.therapeutic_area_id,\n                COUNT(*) AS record_count,\n                SUM(a.fractional_weight) AS fractional_record_count,\n                COUNT(DISTINCT e.canonical_disease_id) AS disease_count,\n                COUNT(DISTINCT e.canonical_target_id) AS target_count,\n                COUNT(e.applicable_score) AS defined_score_count,\n                QUANTILE_CONT(e.applicable_score, 0.5) AS applicable_score_median,\n                QUANTILE_CONT(e.applicable_score, 0.9) AS applicable_score_q90,\n                COUNT(e.old_score) AS old_defined_score_count,\n                QUANTILE_CONT(e.old_score, 0.5) AS old_score_median,\n                QUANTILE_CONT(e.old_score, 0.9) AS old_score_q90,\n                COUNT(e.new_score) AS new_defined_score_count,\n                QUANTILE_CONT(e.new_score, 0.5) AS new_score_median,\n                QUANTILE_CONT(e.new_score, 0.9) AS new_score_q90,\n                COUNT(e.applicable_native_rank) AS native_rank_context_count,\n                COUNT(*) - COUNT(e.applicable_native_rank) AS missing_native_rank_context_count,\n                SUM(CASE WHEN e.applicable_top_rank_context THEN a.fractional_weight ELSE CAST(0 AS DECIMAL(38,18)) END) AS top_rank_context_count,\n                SUM(CASE WHEN e.applicable_top_rank_context THEN a.fractional_weight ELSE CAST(0 AS DECIMAL(38,18)) END) / NULLIF(SUM(CASE WHEN e.applicable_top_rank_context IS NOT NULL THEN a.fractional_weight ELSE CAST(0 AS DECIMAL(38,18)) END), 0) AS top_rank_context_fraction,\n                CASE WHEN COUNT(*) - COUNT(e.applicable_native_rank) > 0 THEN {missing_status} ELSE 'complete' END AS rank_context_status\n            FROM entrant_exit_context e\n            INNER JOIN stable_areas a USING (canonical_disease_id)\n            GROUP BY e.release_pair, e.context_class, a.therapeutic_area_id\n        )\n        SELECT * FROM overall\n        UNION ALL\n        SELECT * FROM area\n        "
    )
    connection.execute(
        f"\n        INSERT INTO entrant_exit_summary (\n            summary_scope,\n            release_pair,\n            context_class,\n            therapeutic_area_id,\n            record_count,\n            fractional_record_count,\n            disease_count,\n            target_count,\n            defined_score_count,\n            old_defined_score_count,\n            new_defined_score_count,\n            native_rank_context_count,\n            missing_native_rank_context_count,\n            top_rank_context_count,\n            rank_context_status\n        )\n        SELECT\n            'overall',\n            g.release_pair,\n            g.context_class,\n            NULL,\n            0,\n            0.0,\n            0,\n            0,\n            0,\n            0,\n            0,\n            0,\n            0,\n            0.0,\n            'complete'\n        FROM (\n            SELECT DISTINCT p.release_pair, c.context_class\n            FROM pair_score_revisions p\n            CROSS JOIN (VALUES {context_values}) c(context_class)\n        ) g\n        LEFT JOIN entrant_exit_summary s\n          ON s.summary_scope = 'overall'\n         AND s.release_pair = g.release_pair\n         AND s.context_class = g.context_class\n        WHERE s.release_pair IS NULL\n        "
    )

def export_table(
    connection: Any, table: str, order_by: str, destination: Path, parquet: bool
) -> None:
    destination.parent.mkdir(parents = True, exist_ok = True)
    destination_sql = sql_literal(destination.as_posix())
    query = f"SELECT * FROM {sql_identifier(table)} ORDER BY {order_by}"
    if (parquet):
        connection.execute(
            f"COPY ({query}) TO {destination_sql} (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000)"
        )
    else:
        connection.execute(
            f"COPY ({query}) TO {destination_sql} (FORMAT CSV, HEADER TRUE, DELIMITER ',')"
        )

def export_specs() -> dict[str, tuple[str, str, bool]]:
    return {
        "pair_score_revisions": (
            "pair_score_revisions",
            "release_pair, support_stratum, canonical_disease_id, canonical_target_id",
            True,
        ),
        "score_summary": (
            "score_summary",
            "release_pair, support_stratum, summary_scope, mapping_context_group, evidence_count_status_group",
            False,
        ),
        "implementation_scale_diagnostic": (
            "implementation_scale_diagnostic",
            "release_pair, support_stratum, summary_scope, mapping_context_group, evidence_count_status_group",
            False,
        ),
        "score_tier_summary": (
            "score_tier_summary",
            "release_pair, support_stratum, old_score_tier",
            False,
        ),
        "score_tier_area_distribution": (
            "score_tier_area_distribution",
            "release_pair, support_stratum, old_score_tier, therapeutic_area_id",
            False,
        ),
        "therapeutic_area_summary": (
            "therapeutic_area_summary",
            "release_pair, support_stratum, therapeutic_area_id",
            False,
        ),
        "concentration_diagnostics": (
            "concentration_diagnostics",
            "release_pair, support_stratum",
            False,
        ),
        "entrant_exit_summary": (
            "entrant_exit_summary",
            "release_pair, summary_scope, context_class, therapeutic_area_id NULLS FIRST",
            False,
        ),
        "evidence_domain_tail_context": (
            "evidence_domain_tail_context",
            "release_pair, support_stratum, endpoint_role, endpoint_id",
            False,
        ),
        "trajectory_table": (
            "trajectory_summary",
            "release_span, score_tier_group, trajectory",
            False,
        ),
    }

def export_stage(
    connection: Any, stage: Path, configuration: Mapping[str, Any]
) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    for (output_name, (table, order_by, parquet)) in (export_specs().items()):
        configured_name = Path(configuration["outputs"][output_name]).name
        destination = stage / configured_name
        export_table(connection, table, order_by, destination, parquet)
        paths[output_name] = destination
    return paths

def create_analysis_tables(connection: Any, configuration: Mapping[str, Any]) -> None:
    create_stable_area_table(connection, configuration)
    create_pair_revision_table(connection, configuration)
    create_native_rank_context_table(connection, configuration)
    create_entrant_exit_context_table(connection, configuration)
    create_evidence_domain_context_tables(connection, configuration)
    create_score_summary_table(connection, configuration)
    create_implementation_scale_diagnostic_table(connection, configuration)
    create_score_tier_summary_table(connection, configuration)
    create_score_tier_area_distribution_table(connection)
    create_therapeutic_area_summary_table(connection, configuration)
    create_concentration_tables(connection, configuration)
    create_trajectory_table(connection, configuration)
    create_entrant_exit_summary_table(connection, configuration)
