from __future__ import annotations
import csv
from pathlib import Path
from typing import Any
from typing import Iterator
from typing import Mapping
from typing import Sequence
from ot_release_sensitivity.pipeline.ranking import finite_float
from ot_release_sensitivity.pipeline.ranking import jaccard_overlap
from ot_release_sensitivity.pipeline.ranking import kendall_tau_b
from ot_release_sensitivity.pipeline.ranking import linear_quantile
from ot_release_sensitivity.pipeline.ranking import normalised_rank_displacements
from ot_release_sensitivity.pipeline.ranking import parse_boolean
from ot_release_sensitivity.pipeline.ranking import spearman_correlation
from ot_release_sensitivity.pipeline.ranking import tie_aware_top_ranked_set
from ot_release_sensitivity.pipeline.ranking import tie_aware_top_set
from .sensitivity_types import CrossLayerSensitivityBuildError
from .sensitivity_metrics import (
    joined_identifiers,
    split_identifiers,
    proportion,
    calculate_fixed_panel_metrics,
    classify_persistence,
    release_note_alignment_class,
)

def sql_literal(value: Any) -> str:
    return "'" + str(value).replace("'", "''") + "'"

def sql_identifier(value: str) -> str:
    return '"' + str(value).replace('"', '""') + '"'

def resolve_project_path(project_root: Path, value: str) -> Path:
    path = (project_root / Path(value)).resolve()
    if (not path.is_relative_to(project_root.resolve())):
        raise CrossLayerSensitivityBuildError(
            f"Project-relative path escapes the project root: {value}"
        )
    return path

def relation_columns(connection: Any, relation: str) -> set[str]:
    return {
        str(row[1])
        for (row) in (
            connection.execute(f"PRAGMA table_info({sql_literal(relation)})").fetchall()
        )
    }

def create_release_change_matrix_table(
    connection: Any, path: Path, configuration: Mapping[str, Any]
) -> None:
    contract = configuration["input_reading"]["release_change_matrix"]
    expected_header = list(contract["expected_header"])
    with path.open("r", encoding = "utf-8", newline = "") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise CrossLayerSensitivityBuildError(
                "Release-change matrix is empty"
            ) from exc
        if (header != expected_header):
            raise CrossLayerSensitivityBuildError(
                "Release-change matrix header does not match the frozen contract"
            )
        rows: list[list[str]] = []
        for (source_row) in (reader):
            if (len(source_row) != len(expected_header)):
                raise CrossLayerSensitivityBuildError(
                    "Release-change matrix row width is outside the frozen reader contract"
                )
            rows.append(source_row)
    if (not rows):
        raise CrossLayerSensitivityBuildError(
            "Release-change matrix contains no records"
        )
    fields = ", ".join(
        (f"{sql_identifier(field)} VARCHAR" for (field) in (expected_header))
    )
    parameters = ", ".join(("?" for (_) in (expected_header)))
    connection.execute(f"CREATE TABLE release_change_matrix ({fields})")
    connection.executemany(
        f"INSERT INTO release_change_matrix VALUES ({parameters})", rows
    )

def create_input_views(
    connection: Any, paths: Mapping[str, Path], configuration: Mapping[str, Any]
) -> None:
    for (name, path) in (paths.items()):
        if (name == "release_change_matrix"):
            create_release_change_matrix_table(connection, path, configuration)
        else:
            reader = (
                "read_parquet"
                if (path.suffix.casefold() == ".parquet")
                else "read_csv_auto"
            )
            options = (
                "" if (reader == "read_parquet") else ", header=true, sample_size=-1"
            )
            connection.execute(
                f"CREATE VIEW {sql_identifier(name)} AS SELECT * FROM {reader}({sql_literal(path.as_posix())}{options})"
            )
        missing = sorted(
            set(configuration["input_contract"][name])
            - relation_columns(connection, name)
        )
        if (missing):
            raise CrossLayerSensitivityBuildError(
                f"Input {name} lacks required fields: {missing}"
            )

def query_dicts(connection: Any, query: str) -> list[dict[str, Any]]:
    cursor = connection.execute(query)
    fields = [str(item[0]) for (item) in (cursor.description)]
    return [dict(zip(fields, row)) for (row) in (cursor.fetchall())]

def query_grouped_dicts(
    connection: Any, query: str, key_fields: Sequence[str], batch_size: int = 10000
) -> Iterator[tuple[tuple[Any, ...], list[dict[str, Any]]]]:
    cursor = connection.execute(query)
    fields = [str(item[0]) for (item) in (cursor.description)]
    current_key: tuple[Any, ...] | None = None
    current_rows: list[dict[str, Any]] = []
    while (True):
        batch = cursor.fetchmany(batch_size)
        if (not batch):
            break
        for (raw) in (batch):
            row = dict(zip(fields, raw))
            key = tuple((row[field] for (field) in (key_fields)))
            if (current_key is not None and key != current_key):
                yield (current_key, current_rows)
                current_rows = []
            current_key = key
            current_rows.append(row)
    if (current_key is not None):
        yield (current_key, current_rows)

def create_stable_area_table(connection: Any, configuration: Mapping[str, Any]) -> None:
    reference = sql_literal(configuration["therapeutic_areas"]["reference_release"])
    stability = sql_literal(configuration["therapeutic_areas"]["stability_value"])
    connection.execute(
        f"\n        CREATE TEMP TABLE stable_areas AS\n        SELECT\n            canonical_disease_id,\n            canonical_therapeutic_area_id AS therapeutic_area_id,\n            MAX(canonical_therapeutic_area_label) AS therapeutic_area_label,\n            MAX(CAST(fractional_weight AS DOUBLE)) AS fractional_weight\n        FROM therapeutic_area_mapping\n        WHERE release = {reference}\n          AND UPPER(CAST(therapeutic_area_stability AS VARCHAR)) = {stability}\n          AND TRY_CAST(primary_area_contrast_eligible AS BOOLEAN)\n        GROUP BY canonical_disease_id, canonical_therapeutic_area_id\n        "
    )

def create_fixed_revision_tables(
    connection: Any, configuration: Mapping[str, Any]
) -> None:
    support = sql_literal(configuration["ranking"]["support_name"])
    endpoint_level = sql_literal(
        str(configuration["evidence_domain"]["endpoint_level"]).casefold()
    )
    endpoint_id = sql_literal(
        str(configuration["evidence_domain"]["endpoint_id"]).casefold()
    )
    reference = sql_literal(configuration["therapeutic_areas"]["reference_release"])
    connection.execute(
        "\n        CREATE TEMP TABLE primary_panel_keys AS\n        SELECT *\n        FROM disease_ranking_panels\n        QUALIFY COUNT(*) OVER (PARTITION BY release_pair, canonical_disease_id) = 1\n        "
    )
    duplicate_primary = connection.execute(
        "\n        SELECT COUNT(*)\n        FROM (\n            SELECT release_pair, canonical_disease_id\n            FROM disease_ranking_panels\n            GROUP BY release_pair, canonical_disease_id\n            HAVING COUNT(*) <> 1\n        )\n        "
    ).fetchone()[0]
    if (int(duplicate_primary) != 0):
        raise CrossLayerSensitivityBuildError(
            "Ranking sensitivity disease-ranking registry contains duplicate panel keys"
        )
    connection.execute(
        f"\n        CREATE TEMP TABLE stable_eligible_panels AS\n        SELECT\n            r.release_pair,\n            r.canonical_disease_id,\n            r.canonical_label,\n            r.canonical_therapeutic_area_ids,\n            CAST(r.stable_fixed_target_count AS BIGINT) AS stable_fixed_target_count,\n            CAST(r.fixed_target_count AS BIGINT) AS primary_fixed_target_count\n        FROM fixed_ranking_panel_registry r\n        JOIN primary_panel_keys p USING (release_pair, canonical_disease_id)\n        WHERE r.support_name = {support}\n          AND TRY_CAST(r.primary_panel_eligible AS BOOLEAN)\n          AND TRY_CAST(r.stable_panel_eligible AS BOOLEAN)\n        "
    )
    connection.execute(
        f"\n        CREATE TEMP TABLE primary_pair_revisions AS\n        SELECT\n            m.release_pair,\n            m.old_release,\n            m.new_release,\n            m.canonical_disease_id,\n            m.canonical_target_id,\n            CAST(m.old_score AS DOUBLE) AS old_score,\n            CAST(m.new_score AS DOUBLE) AS new_score,\n            CAST(m.new_score AS DOUBLE) - CAST(m.old_score AS DOUBLE) AS signed_change,\n            ABS(CAST(m.new_score AS DOUBLE) - CAST(m.old_score AS DOUBLE)) AS absolute_change\n        FROM fixed_ranking_panel_members m\n        JOIN primary_panel_keys p USING (release_pair, canonical_disease_id)\n        WHERE m.support_name = {support}\n        "
    )
    connection.execute(
        f"\n        CREATE TEMP TABLE stable_revision_base AS\n        SELECT\n            m.release_pair,\n            m.old_release,\n            m.new_release,\n            m.canonical_disease_id,\n            p.canonical_label,\n            p.canonical_therapeutic_area_ids,\n            m.canonical_target_id,\n            CAST(m.old_score AS DOUBLE) AS old_score,\n            CAST(m.new_score AS DOUBLE) AS new_score,\n            TRUE AS stable_entity_eligible\n        FROM fixed_ranking_panel_members m\n        JOIN stable_eligible_panels p USING (release_pair, canonical_disease_id)\n        WHERE m.support_name = {support}\n          AND TRY_CAST(m.stable_entity_member AS BOOLEAN)\n        "
    )
    connection.execute(
        "\n        CREATE TABLE stable_entity_pair_revisions AS\n        WITH ranked AS (\n            SELECT\n                *,\n                COUNT(*) OVER (PARTITION BY release_pair, canonical_disease_id) AS fixed_target_count,\n                RANK() OVER (\n                    PARTITION BY release_pair, canonical_disease_id\n                    ORDER BY old_score DESC\n                ) + (COUNT(*) OVER (\n                    PARTITION BY release_pair, canonical_disease_id, old_score\n                ) - 1) / 2.0 AS old_rank,\n                RANK() OVER (\n                    PARTITION BY release_pair, canonical_disease_id\n                    ORDER BY new_score DESC\n                ) + (COUNT(*) OVER (\n                    PARTITION BY release_pair, canonical_disease_id, new_score\n                ) - 1) / 2.0 AS new_rank\n            FROM stable_revision_base\n        )\n        SELECT\n            release_pair,\n            old_release,\n            new_release,\n            canonical_disease_id,\n            canonical_label,\n            canonical_therapeutic_area_ids,\n            canonical_target_id,\n            old_score,\n            new_score,\n            new_score - old_score AS signed_change,\n            ABS(new_score - old_score) AS absolute_change,\n            old_rank,\n            new_rank,\n            ABS(new_rank - old_rank) / (fixed_target_count - 1) AS normalised_absolute_rank_displacement,\n            fixed_target_count,\n            stable_entity_eligible,\n            'STABLE_ENTITY_FIXED' AS support_denominator\n        FROM ranked\n        "
    )
    connection.execute(
        f"\n        CREATE TEMP TABLE disease_labels AS\n        SELECT canonical_disease_id, MAX(canonical_label) AS canonical_label\n        FROM disease_crosswalk\n        WHERE release = {reference} AND canonical_disease_id IS NOT NULL AND canonical_disease_id <> ''\n        GROUP BY canonical_disease_id\n        "
    )
    connection.execute(
        f"\n        CREATE TEMP TABLE evidence_eligible_panels AS\n        SELECT\n            r.release_pair,\n            r.canonical_disease_id,\n            COALESCE(l.canonical_label, '') AS canonical_label,\n            r.canonical_therapeutic_area_ids,\n            CAST(r.fixed_target_count AS BIGINT) AS fixed_target_count\n        FROM evidence_domain_panel_registry r\n        LEFT JOIN disease_labels l USING (canonical_disease_id)\n        WHERE LOWER(CAST(r.endpoint_level AS VARCHAR)) = {endpoint_level}\n          AND LOWER(CAST(r.endpoint_id AS VARCHAR)) = {endpoint_id}\n          AND TRY_CAST(r.panel_eligible AS BOOLEAN)\n        "
    )
    pair_queries: list[str] = []
    for (pair) in (configuration["release_pairs"]):
        pair_label = sql_literal(pair["label"])
        old_release = sql_literal(pair["old_release"])
        new_release = sql_literal(pair["new_release"])
        pair_queries.append(
            f"\n            SELECT\n                {pair_label} AS release_pair,\n                {old_release} AS old_release,\n                {new_release} AS new_release,\n                o.canonical_disease_id,\n                p.canonical_label,\n                p.canonical_therapeutic_area_ids,\n                o.canonical_target_id,\n                CAST(o.score AS DOUBLE) AS old_score,\n                CAST(n.score AS DOUBLE) AS new_score,\n                TRY_CAST(o.stable_entity_eligible AS BOOLEAN)\n                  AND TRY_CAST(n.stable_entity_eligible AS BOOLEAN) AS stable_entity_eligible,\n                o.source_object_relative_path AS old_source_object_relative_path,\n                n.source_object_relative_path AS new_source_object_relative_path\n            FROM evidence_domain_support o\n            JOIN evidence_domain_support n\n              ON LOWER(CAST(o.endpoint_level AS VARCHAR)) = LOWER(CAST(n.endpoint_level AS VARCHAR))\n             AND LOWER(CAST(o.endpoint_id AS VARCHAR)) = LOWER(CAST(n.endpoint_id AS VARCHAR))\n             AND o.canonical_disease_id = n.canonical_disease_id\n             AND o.canonical_target_id = n.canonical_target_id\n            JOIN evidence_eligible_panels p\n              ON p.release_pair = {pair_label}\n             AND p.canonical_disease_id = o.canonical_disease_id\n            WHERE LOWER(CAST(o.endpoint_level AS VARCHAR)) = {endpoint_level}\n              AND LOWER(CAST(o.endpoint_id AS VARCHAR)) = {endpoint_id}\n              AND o.release = {old_release}\n              AND n.release = {new_release}\n              AND TRY_CAST(o.domain_support_eligible AS BOOLEAN)\n              AND TRY_CAST(n.domain_support_eligible AS BOOLEAN)\n            "
        )
    connection.execute(
        "CREATE TEMP TABLE evidence_revision_base AS "
        + " UNION ALL ".join(pair_queries)
    )
    connection.execute(
        "\n        CREATE TABLE evidence_domain_pair_revisions AS\n        WITH ranked AS (\n            SELECT\n                *,\n                COUNT(*) OVER (PARTITION BY release_pair, canonical_disease_id) AS fixed_target_count,\n                RANK() OVER (\n                    PARTITION BY release_pair, canonical_disease_id\n                    ORDER BY old_score DESC\n                ) + (COUNT(*) OVER (\n                    PARTITION BY release_pair, canonical_disease_id, old_score\n                ) - 1) / 2.0 AS old_rank,\n                RANK() OVER (\n                    PARTITION BY release_pair, canonical_disease_id\n                    ORDER BY new_score DESC\n                ) + (COUNT(*) OVER (\n                    PARTITION BY release_pair, canonical_disease_id, new_score\n                ) - 1) / 2.0 AS new_rank\n            FROM evidence_revision_base\n        )\n        SELECT\n            'DATATYPE' AS endpoint_level,\n            'genetic_association' AS endpoint_id,\n            release_pair,\n            old_release,\n            new_release,\n            canonical_disease_id,\n            canonical_label,\n            canonical_therapeutic_area_ids,\n            canonical_target_id,\n            old_score,\n            new_score,\n            new_score - old_score AS signed_change,\n            ABS(new_score - old_score) AS absolute_change,\n            old_rank,\n            new_rank,\n            ABS(new_rank - old_rank) / (fixed_target_count - 1) AS normalised_absolute_rank_displacement,\n            fixed_target_count,\n            stable_entity_eligible,\n            old_source_object_relative_path,\n            new_source_object_relative_path,\n            'GENETIC_ASSOCIATION_FIXED' AS support_denominator\n        FROM ranked\n        "
    )

def score_summary_rows(connection: Any, table: str, layer: str) -> list[dict[str, Any]]:
    return query_dicts(
        connection,
        f"\n        SELECT\n            release_pair,\n            {sql_literal(layer)} AS analysis_layer,\n            {sql_literal(layer)} AS support_denominator,\n            COUNT(*) AS eligible_target_pair_count,\n            COUNT(DISTINCT canonical_disease_id) AS eligible_disease_panel_count,\n            COUNT(DISTINCT canonical_target_id) AS distinct_target_count,\n            MEDIAN(signed_change) AS median_signed_change,\n            MEDIAN(absolute_change) AS median_absolute_change,\n            QUANTILE_CONT(absolute_change, 0.90) AS q90_absolute_change,\n            QUANTILE_CONT(absolute_change, 0.95) AS q95_absolute_change,\n            MAX(absolute_change) AS maximum_absolute_change,\n            AVG(CASE WHEN absolute_change = 0 THEN 1.0 ELSE 0.0 END) AS exact_unchanged_fraction,\n            AVG(CASE WHEN absolute_change >= 0.01 THEN 1.0 ELSE 0.0 END) AS absolute_change_at_least_0_01_fraction,\n            AVG(CASE WHEN absolute_change >= 0.05 THEN 1.0 ELSE 0.0 END) AS absolute_change_at_least_0_05_fraction,\n            AVG(CASE WHEN absolute_change >= 0.10 THEN 1.0 ELSE 0.0 END) AS absolute_change_at_least_0_10_fraction\n        FROM {sql_identifier(table)}\n        GROUP BY release_pair\n        ORDER BY release_pair\n        ",
    )

def build_fixed_metric_rows(
    connection: Any,
    table: str,
    configuration: Mapping[str, Any],
    layer: str,
    primary_by_key: Mapping[tuple[str, str], Mapping[str, Any]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    query = f"\n        SELECT\n            release_pair,\n            old_release,\n            new_release,\n            canonical_disease_id,\n            canonical_label,\n            canonical_therapeutic_area_ids,\n            canonical_target_id,\n            old_score,\n            new_score,\n            stable_entity_eligible\n        FROM {sql_identifier(table)}\n        ORDER BY release_pair, canonical_disease_id, canonical_target_id\n    "
    for (key, panel_rows) in (query_grouped_dicts(
        connection, query, ["release_pair", "canonical_disease_id"]
    )):
        metrics = calculate_fixed_panel_metrics(panel_rows, configuration, layer)
        primary = primary_by_key.get((str(key[0]), str(key[1])))
        old_domain = split_identifiers(metrics["old_top_ranked_target_ids"])
        new_domain = split_identifiers(metrics["new_top_ranked_target_ids"])
        old_overall = (
            split_identifiers(primary.get("old_fixed_top_ranked_target_ids"))
            if (primary)
            else set()
        )
        new_overall = (
            split_identifiers(primary.get("new_fixed_top_ranked_target_ids"))
            if (primary)
            else set()
        )
        domain_change = old_domain ^ new_domain
        overall_change = old_overall ^ new_overall
        metrics.update(
            {
                "overall_primary_panel_available": primary is not None,
                "overall_top_ranked_set_changed": None
                if (primary is None)
                else parse_boolean(primary["fixed_top_ranked_set_changed"]),
                "changed_top_target_overlap_jaccard_with_overall": None
                if (primary is None or not domain_change or (not overall_change))
                else len(domain_change & overall_change)
                / len(domain_change | overall_change),
                "corroborates_overall_top_change": None
                if (primary is None)
                else bool(
                    parse_boolean(primary["fixed_top_ranked_set_changed"])
                    and metrics["top_ranked_set_changed"]
                    and domain_change & overall_change
                ),
                "descriptive_not_causal": True,
            }
        )
        rows.append(metrics)
    return rows

def native_change_reason(
    old_set: set[str] | None,
    new_set: set[str] | None,
    old_members: set[str],
    new_members: set[str],
) -> str:
    if (old_set is None or new_set is None):
        return "UNRESOLVED"
    changed = old_set ^ new_set
    if (not changed):
        return "UNCHANGED"
    reasons: set[str] = set()
    for (target) in (changed):
        if (target in old_members and target in new_members):
            reasons.add("PERSISTENT_PAIR_RERANKING")
        elif (target in new_members):
            reasons.add("ENTRANT")
        elif (target in old_members):
            reasons.add("EXIT")
        else:
            reasons.add("UNRESOLVED")
    return next(iter(reasons)) if (len(reasons) == 1) else "MIXED"

def calculate_stable_native_panel(
    rows: Sequence[Mapping[str, Any]], configuration: Mapping[str, Any]
) -> dict[str, Any]:
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
    old_members = set(old_scores)
    new_members = set(new_scores)
    minimum = finite_float(configuration["ranking"]["minimum_positive_top_set_score"])
    result: dict[str, Any] = {
        "analysis_layer": "STABLE_ENTITY_RELEASE_NATIVE",
        "release_pair": rows[0]["release_pair"],
        "old_release": rows[0]["old_release"],
        "new_release": rows[0]["new_release"],
        "canonical_disease_id": rows[0]["canonical_disease_id"],
        "canonical_label": rows[0]["canonical_label"],
        "canonical_therapeutic_area_ids": rows[0]["canonical_therapeutic_area_ids"],
        "old_native_target_count": len(old_members),
        "new_native_target_count": len(new_members),
        "persistent_member_count": len(old_members & new_members),
        "entrant_member_count": len(new_members - old_members),
        "exit_member_count": len(old_members - new_members),
        "support_denominator": "STABLE_ENTITY_RELEASE_NATIVE",
    }
    for (cutoff) in (configuration["ranking"]["top_set_cutoffs"]):
        old_set = tie_aware_top_set(old_scores, int(cutoff), minimum)
        new_set = tie_aware_top_set(new_scores, int(cutoff), minimum)
        result[f"native_top_{cutoff}_jaccard"] = jaccard_overlap(old_set, new_set)
        result[f"old_native_top_{cutoff}_effective_size"] = (
            None if (old_set is None) else len(old_set)
        )
        result[f"new_native_top_{cutoff}_effective_size"] = (
            None if (new_set is None) else len(new_set)
        )
        result[f"native_top_{cutoff}_changed"] = (
            None if (old_set is None or new_set is None) else old_set != new_set
        )
        result[f"native_top_{cutoff}_change_reason"] = native_change_reason(
            old_set, new_set, old_members, new_members
        )
    old_top = tie_aware_top_ranked_set(old_scores, minimum)
    new_top = tie_aware_top_ranked_set(new_scores, minimum)
    result.update(
        {
            "old_native_top_ranked_target_ids": joined_identifiers(old_top),
            "new_native_top_ranked_target_ids": joined_identifiers(new_top),
            "native_top_ranked_jaccard": jaccard_overlap(old_top, new_top),
            "native_top_ranked_set_changed": None
            if (old_top is None or new_top is None)
            else old_top != new_top,
            "native_top_ranked_change_reason": native_change_reason(
                old_top, new_top, old_members, new_members
            ),
            "descriptive_not_causal": True,
        }
    )
    return result

def build_stable_native_rows(
    connection: Any, configuration: Mapping[str, Any]
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for (pair) in (configuration["release_pairs"]):
        label = sql_literal(pair["label"])
        old_release = sql_literal(pair["old_release"])
        new_release = sql_literal(pair["new_release"])
        query = f"\n            WITH old_members AS (\n                SELECT canonical_disease_id, canonical_target_id, CAST(score AS DOUBLE) AS old_score\n                FROM release_native_panel_members\n                WHERE release = {old_release} AND TRY_CAST(stable_entity_member AS BOOLEAN)\n            ),\n            new_members AS (\n                SELECT canonical_disease_id, canonical_target_id, CAST(score AS DOUBLE) AS new_score\n                FROM release_native_panel_members\n                WHERE release = {new_release} AND TRY_CAST(stable_entity_member AS BOOLEAN)\n            ),\n            combined AS (\n                SELECT\n                    COALESCE(o.canonical_disease_id, n.canonical_disease_id) AS canonical_disease_id,\n                    COALESCE(o.canonical_target_id, n.canonical_target_id) AS canonical_target_id,\n                    o.old_score,\n                    n.new_score,\n                    o.canonical_target_id IS NOT NULL AS old_present,\n                    n.canonical_target_id IS NOT NULL AS new_present\n                FROM old_members o\n                FULL OUTER JOIN new_members n USING (canonical_disease_id, canonical_target_id)\n            )\n            SELECT\n                {label} AS release_pair,\n                {old_release} AS old_release,\n                {new_release} AS new_release,\n                c.canonical_disease_id,\n                p.canonical_label,\n                p.canonical_therapeutic_area_ids,\n                c.canonical_target_id,\n                c.old_score,\n                c.new_score,\n                c.old_present,\n                c.new_present\n            FROM combined c\n            JOIN stable_eligible_panels p\n              ON p.release_pair = {label}\n             AND p.canonical_disease_id = c.canonical_disease_id\n            ORDER BY c.canonical_disease_id, c.canonical_target_id\n        "
        for (key, rows) in (query_grouped_dicts(
            connection, query, ["canonical_disease_id"]
        )):
            if (not rows):
                continue
            output.append(calculate_stable_native_panel(rows, configuration))
    return output

def discover_diagnostic_inputs(project_root, configuration):
    result = {}
    dataset = configuration["diagnostic_datasources"]["dataset_name"]
    releases = sorted(
        {pair["old_release"] for (pair) in (configuration["release_pairs"])}
        | {configuration["release_pairs"][-1]["new_release"]}
    )
    for (release) in (releases):
        directory = Path(project_root) / "data" / "raw" / release / dataset
        paths = sorted(directory.rglob("*.parquet")) if (directory.is_dir()) else []
        if (not paths):
            raise FileNotFoundError(f"Missing datasource partitions: {directory}")
        result[release] = [
            (path.relative_to(project_root).as_posix(), path) for (path) in (paths)
        ]
    return result

def parquet_list_sql(paths: Sequence[Path]) -> str:
    return "[" + ",".join((sql_literal(path.as_posix()) for (path) in (paths))) + "]"

def create_changed_target_table(
    connection: Any, primary_rows: Sequence[Mapping[str, Any]]
) -> None:
    connection.execute(
        "\n        CREATE TEMP TABLE changed_primary_targets (\n            release_pair VARCHAR,\n            canonical_disease_id VARCHAR,\n            canonical_target_id VARCHAR,\n            old_top_ranked_member BOOLEAN,\n            new_top_ranked_member BOOLEAN\n        )\n        "
    )
    rows: list[tuple[Any, ...]] = []
    for (panel) in (primary_rows):
        if (not parse_boolean(panel["fixed_top_ranked_set_changed"])):
            continue
        old_set = split_identifiers(panel["old_fixed_top_ranked_target_ids"])
        new_set = split_identifiers(panel["new_fixed_top_ranked_target_ids"])
        for (target) in (sorted(old_set | new_set)):
            rows.append(
                (
                    panel["release_pair"],
                    panel["canonical_disease_id"],
                    target,
                    target in old_set,
                    target in new_set,
                )
            )
    if (rows):
        connection.executemany(
            "INSERT INTO changed_primary_targets VALUES (?, ?, ?, ?, ?)", rows
        )

def diagnostic_projection_query(
    project_root: Path,
    release: str,
    paths: Sequence[tuple[str, Path]],
    configuration: Mapping[str, Any],
) -> str:
    field_contract = configuration["diagnostic_datasources"]["fields"][release]
    relation = f"read_parquet({parquet_list_sql([path for ((relative, path)) in (paths)])}, filename=true, union_by_name=true)"
    endpoint = sql_identifier(field_contract["endpoint_field"])
    filters = []
    if ("filter_field" in field_contract):
        filters.append(
            f"CAST(r.{sql_identifier(field_contract['filter_field'])} AS VARCHAR) = {sql_literal(field_contract['filter_value'])}"
        )
    family_cases: list[str] = []
    endpoint_values: list[str] = []
    for (family) in (configuration["diagnostic_datasources"]["families"]):
        endpoint_value = str(family["release_endpoints"][release])
        endpoint_values.append(endpoint_value)
        family_cases.append(
            f"WHEN CAST(r.{endpoint} AS VARCHAR) = {sql_literal(endpoint_value)} THEN {sql_literal(family['family'])}"
        )
    filters.append(
        f"CAST(r.{endpoint} AS VARCHAR) IN ({','.join((sql_literal(value) for (value) in (sorted(set(endpoint_values)))))})"
    )
    source_cases = [
        f"WHEN REPLACE(CAST(r.filename AS VARCHAR), {sql_literal(chr(92))}, '/') = {sql_literal(path.as_posix())} THEN {sql_literal(relative)}"
        for ((relative, path)) in (paths)
    ]
    source_path_expression = "CASE " + " ".join(source_cases) + " ELSE NULL END"
    return f"\n        SELECT\n            {sql_literal(release)} AS release,\n            CASE {' '.join(family_cases)} ELSE 'UNSELECTED' END AS diagnostic_family,\n            CAST(r.{endpoint} AS VARCHAR) AS endpoint_id,\n            d.canonical_disease_id,\n            t.canonical_target_id,\n            TRY_CAST(r.{sql_identifier(field_contract['score'])} AS DOUBLE) AS score,\n            TRY_CAST(r.{sql_identifier(field_contract['evidence_count'])} AS BIGINT) AS evidence_count,\n            {source_path_expression} AS source_object_relative_path\n        FROM {relation} r\n        JOIN disease_crosswalk d\n          ON d.release = {sql_literal(release)}\n         AND d.source_disease_id = CAST(r.{sql_identifier(field_contract['disease_id'])} AS VARCHAR)\n         AND TRY_CAST(d.primary_support_eligible AS BOOLEAN)\n        JOIN target_crosswalk t\n          ON t.release = {sql_literal(release)}\n         AND t.source_target_id = CAST(r.{sql_identifier(field_contract['target_id'])} AS VARCHAR)\n         AND TRY_CAST(t.primary_support_eligible AS BOOLEAN)\n        JOIN (\n            SELECT DISTINCT canonical_disease_id, canonical_target_id\n            FROM changed_primary_targets\n        ) c USING (canonical_disease_id, canonical_target_id)\n        WHERE {' AND '.join(filters)}\n    "

def create_diagnostic_tables(
    connection: Any,
    project_root: Path,
    source_paths: Mapping[str, Sequence[tuple[str, Path]]],
    configuration: Mapping[str, Any],
) -> dict[str, int]:
    projections = [
        diagnostic_projection_query(project_root, release, paths, configuration)
        for ((release, paths)) in (sorted(source_paths.items()))
    ]
    connection.execute(
        "CREATE TEMP TABLE diagnostic_support_unchecked AS "
        + " UNION ALL ".join(projections)
    )
    duplicate_count = connection.execute(
        "\n        SELECT COALESCE(SUM(n - 1), 0)\n        FROM (\n            SELECT release, diagnostic_family, canonical_disease_id, canonical_target_id, COUNT(*) AS n\n            FROM diagnostic_support_unchecked\n            GROUP BY release, diagnostic_family, canonical_disease_id, canonical_target_id\n            HAVING COUNT(*) > 1\n        )\n        "
    ).fetchone()[0]
    if (int(duplicate_count) != 0):
        raise CrossLayerSensitivityBuildError(
            "Selected diagnostic datasource support contains duplicate canonical keys"
        )
    connection.execute(
        "CREATE TEMP TABLE diagnostic_support AS SELECT * FROM diagnostic_support_unchecked"
    )
    family_queries: list[str] = []
    families = {
        str(item["family"]): item
        for (item) in (configuration["diagnostic_datasources"]["families"])
    }
    threshold = finite_float(
        configuration["scores"]["large_diagnostic_change_threshold"]
    )
    for (pair) in (configuration["release_pairs"]):
        for (family_name, family) in (sorted(families.items())):
            comparable = str(pair["label"]) in set(family["comparable_release_pairs"])
            comparison_status = (
                "COMPARABLE"
                if (comparable)
                else str(
                    configuration["diagnostic_datasources"]["definition_break_status"]
                )
            )
            old_endpoint = str(family["release_endpoints"][pair["old_release"]])
            new_endpoint = str(family["release_endpoints"][pair["new_release"]])
            family_queries.append(
                f"\n                SELECT\n                    c.release_pair,\n                    {sql_literal(pair['old_release'])} AS old_release,\n                    {sql_literal(pair['new_release'])} AS new_release,\n                    c.canonical_disease_id,\n                    c.canonical_target_id,\n                    c.old_top_ranked_member,\n                    c.new_top_ranked_member,\n                    CASE\n                        WHEN c.old_top_ranked_member AND c.new_top_ranked_member THEN 'PERSISTENT_TOP_RANKED'\n                        WHEN c.old_top_ranked_member THEN 'OLD_TOP_RANKED_ONLY'\n                        WHEN c.new_top_ranked_member THEN 'NEW_TOP_RANKED_ONLY'\n                        ELSE 'UNRESOLVED'\n                    END AS primary_top_ranked_membership,\n                    {sql_literal(family_name)} AS diagnostic_family,\n                    {sql_literal(old_endpoint)} AS old_endpoint_id,\n                    {sql_literal(new_endpoint)} AS new_endpoint_id,\n                    {sql_literal(comparison_status)} AS comparison_status,\n                    o.score AS old_score,\n                    n.score AS new_score,\n                    o.evidence_count AS old_evidence_count,\n                    n.evidence_count AS new_evidence_count,\n                    CASE\n                        WHEN o.canonical_target_id IS NOT NULL AND n.canonical_target_id IS NOT NULL THEN 'BOTH_PRESENT'\n                        WHEN o.canonical_target_id IS NOT NULL THEN 'OLD_ONLY'\n                        WHEN n.canonical_target_id IS NOT NULL THEN 'NEW_ONLY'\n                        ELSE 'ABSENT_BOTH'\n                    END AS diagnostic_presence,\n                    CASE WHEN {str(comparable).upper()} AND o.score IS NOT NULL AND n.score IS NOT NULL THEN n.score - o.score END AS signed_change,\n                    CASE WHEN {str(comparable).upper()} AND o.score IS NOT NULL AND n.score IS NOT NULL THEN ABS(n.score - o.score) END AS absolute_change,\n                    CASE WHEN {str(comparable).upper()} AND o.score IS NOT NULL AND n.score IS NOT NULL THEN ABS(n.score - o.score) >= {threshold} END AS large_diagnostic_change,\n                    o.source_object_relative_path AS old_source_object_relative_path,\n                    n.source_object_relative_path AS new_source_object_relative_path,\n                    TRUE AS diagnostic_only_not_causal\n                FROM changed_primary_targets c\n                LEFT JOIN diagnostic_support o\n                  ON o.release = {sql_literal(pair['old_release'])}\n                 AND o.diagnostic_family = {sql_literal(family_name)}\n                 AND o.canonical_disease_id = c.canonical_disease_id\n                 AND o.canonical_target_id = c.canonical_target_id\n                LEFT JOIN diagnostic_support n\n                  ON n.release = {sql_literal(pair['new_release'])}\n                 AND n.diagnostic_family = {sql_literal(family_name)}\n                 AND n.canonical_disease_id = c.canonical_disease_id\n                 AND n.canonical_target_id = c.canonical_target_id\n                WHERE c.release_pair = {sql_literal(pair['label'])}\n                "
            )
    if (family_queries):
        connection.execute(
            "CREATE TABLE diagnostic_target_context AS "
            + " UNION ALL ".join(family_queries)
        )
    else:
        connection.execute(
            "\n            CREATE TABLE diagnostic_target_context AS\n            SELECT NULL::VARCHAR AS release_pair WHERE FALSE\n            "
        )
    return {"diagnostic_duplicate_key_count": int(duplicate_count)}

def build_diagnostic_summary_rows(connection: Any) -> list[dict[str, Any]]:
    return query_dicts(
        connection,
        "\n        SELECT\n            release_pair,\n            old_release,\n            new_release,\n            canonical_disease_id,\n            diagnostic_family,\n            MAX(comparison_status) AS comparison_status,\n            COUNT(*) AS changed_top_target_count,\n            COUNT(*) FILTER (WHERE diagnostic_presence <> 'ABSENT_BOTH') AS target_count_with_diagnostic_support,\n            COUNT(*) FILTER (WHERE diagnostic_presence = 'BOTH_PRESENT') AS target_count_with_both_release_support,\n            COUNT(*) FILTER (WHERE diagnostic_presence IN ('OLD_ONLY', 'NEW_ONLY')) AS target_count_with_presence_change,\n            COUNT(*) FILTER (WHERE large_diagnostic_change) AS target_count_with_large_comparable_change,\n            CASE\n                WHEN COUNT(*) FILTER (WHERE comparison_status = 'COMPARABLE' AND absolute_change IS NOT NULL) = 0 THEN NULL\n                ELSE COUNT(*) FILTER (WHERE large_diagnostic_change)::DOUBLE\n                  / COUNT(*) FILTER (WHERE comparison_status = 'COMPARABLE' AND absolute_change IS NOT NULL)\n            END AS large_comparable_change_fraction,\n            CASE\n                WHEN COUNT(*) FILTER (WHERE diagnostic_presence <> 'ABSENT_BOTH') > 0 THEN 'WITH_DIAGNOSTIC_EVIDENCE'\n                ELSE 'WITHOUT_DIAGNOSTIC_EVIDENCE'\n            END AS evidence_group,\n            BOOL_OR(COALESCE(large_diagnostic_change, FALSE)) AS any_large_comparable_change,\n            TRUE AS diagnostic_only_not_causal\n        FROM diagnostic_target_context\n        GROUP BY release_pair, old_release, new_release, canonical_disease_id, diagnostic_family\n        ORDER BY release_pair, canonical_disease_id, diagnostic_family\n        ",
    )

def build_release_note_rows(
    connection: Any, configuration: Mapping[str, Any]
) -> list[dict[str, Any]]:
    matrix_rows = query_dicts(
        connection, "SELECT * FROM release_change_matrix ORDER BY release, category"
    )
    pair_by_release = {
        str(pair["new_release"]): str(pair["label"])
        for (pair) in (configuration["release_pairs"])
    }
    output: list[dict[str, Any]] = []
    for (row) in (matrix_rows):
        release = str(row["release"])
        alignment = release_note_alignment_class(
            str(row["category"]), str(row["evidence_status"]), configuration
        )
        output.append(
            {
                "release": release,
                "release_pair_context": pair_by_release.get(
                    release, "BASELINE_CONTEXT"
                ),
                "category": row["category"],
                "documented_change": row["documented_change"],
                "analytical_relevance": row["analytical_relevance"],
                "source_evidence_status": row["evidence_status"],
                "alignment_class": alignment,
                "observed_evidence": {
                    "DIRECTLY OBSERVABLE IN RELEASE OUTPUT": "UPSTREAM_RELEASE_OUTPUT_OR_SCHEMA_AUDIT",
                    "POSSIBLY RELEVANT": "MEASURED_CONTEXT_WITHOUT_MECHANISM_ASSIGNMENT",
                    "INFRASTRUCTURE ONLY": "OUTSIDE_FROZEN_ANALYTICAL_ENDPOINT",
                    "NOT TESTED": "NO_RELEASE_OUTPUT_TEST_ASSIGNED",
                }[alignment],
                "mechanism_attribution": configuration["release_note_alignment"][
                    "mechanism_attribution"
                ],
                "source_url": row["source_url"],
                "source_snapshot": row["source_snapshot"],
                "causality_boundary": row["causality_boundary"],
                "context_only_not_causal": True,
            }
        )
    return output

def stable_area_memberships(connection: Any) -> dict[str, set[str]]:
    rows = query_dicts(
        connection,
        "SELECT canonical_disease_id, therapeutic_area_id FROM stable_areas ORDER BY canonical_disease_id, therapeutic_area_id",
    )
    result: dict[str, set[str]] = {}
    for (row) in (rows):
        result.setdefault(str(row["canonical_disease_id"]), set()).add(
            str(row["therapeutic_area_id"])
        )
    return result

def build_sensitivity_comparison_rows(
    primary_rows: Sequence[Mapping[str, Any]],
    stable_rows: Sequence[Mapping[str, Any]],
    domain_rows: Sequence[Mapping[str, Any]],
    stable_native_rows: Sequence[Mapping[str, Any]],
    primary_score_rows: Sequence[Mapping[str, Any]],
    stable_score_rows: Sequence[Mapping[str, Any]],
    domain_score_rows: Sequence[Mapping[str, Any]],
    area_memberships: Mapping[str, set[str]],
    configuration: Mapping[str, Any],
) -> list[dict[str, Any]]:
    score_maps = {
        "PRIMARY_OVERALL_FIXED": {
            str(row["release_pair"]): row for (row) in (primary_score_rows)
        },
        "STABLE_ENTITY_FIXED": {
            str(row["release_pair"]): row for (row) in (stable_score_rows)
        },
        "GENETIC_ASSOCIATION_FIXED": {
            str(row["release_pair"]): row for (row) in (domain_score_rows)
        },
    }
    stable_native_map = {
        (str(row["release_pair"]), str(row["canonical_disease_id"])): row
        for (row) in (stable_native_rows)
    }
    layers: list[tuple[str, Sequence[Mapping[str, Any]]]] = [
        ("PRIMARY_OVERALL_FIXED", primary_rows),
        ("STABLE_ENTITY_FIXED", stable_rows),
        ("GENETIC_ASSOCIATION_FIXED", domain_rows),
    ]
    output: list[dict[str, Any]] = []
    for (pair) in ([str(item["label"]) for (item) in (configuration["release_pairs"])]):
        primary_pair = [
            row for (row) in (primary_rows) if (str(row["release_pair"]) == pair)
        ]
        primary_keys = {
            (pair, str(row["canonical_disease_id"])) for (row) in (primary_pair)
        }
        for (layer, source) in (layers):
            rows = [row for (row) in (source) if (str(row["release_pair"]) == pair)]
            score = score_maps[layer].get(pair, {})
            if (layer == "PRIMARY_OVERALL_FIXED"):
                tau_values = [finite_float(row["kendall_tau_b"]) for (row) in (rows)]
                top10 = [
                    finite_float(row["fixed_top_10_jaccard"]) < 1.0
                    for (row) in (rows)
                    if (row.get("fixed_top_10_jaccard") not in {None, ""})
                ]
                top_ranked = [
                    parse_boolean(row["fixed_top_ranked_set_changed"])
                    for (row) in (rows)
                ]
                mapping_values = [
                    str(row["mapping_stability"]) != "NO_MAPPING_CONTEXT_FLAG"
                    for (row) in (rows)
                ]
                native_changed = [
                    row
                    for (row) in (rows)
                    if (parse_boolean(row.get("native_top_ranked_set_changed", False)))
                ]
                entrant_values = [
                    int(row.get("native_top_ranked_entrant_count") or 0) > 0
                    or int(row.get("native_top_ranked_exit_count") or 0) > 0
                    or str(row.get("native_top_ranked_change_reason", ""))
                    in {"ENTRANT", "EXIT", "MIXED"}
                    for (row) in (native_changed)
                ]
                mapping_definition = "Ranking sensitivity panel-level mapping context"
                entrant_definition = "Changed release-native top-ranked panels"
            else:
                tau_values = [finite_float(row["kendall_tau_b"]) for (row) in (rows)]
                top10 = [
                    parse_boolean(row["top_10_changed"])
                    for (row) in (rows)
                    if (row.get("top_10_changed") not in {None, ""})
                ]
                top_ranked = [
                    parse_boolean(row["top_ranked_set_changed"])
                    for (row) in (rows)
                    if (row.get("top_ranked_set_changed") not in {None, ""})
                ]
                if (layer == "STABLE_ENTITY_FIXED"):
                    stable_keys = {
                        (pair, str(row["canonical_disease_id"])) for (row) in (rows)
                    }
                    mapping_values = [
                        key not in stable_keys for (key) in (sorted(primary_keys))
                    ]
                    native_changed_rows = [
                        stable_native_map[key]
                        for (key) in (sorted(stable_keys))
                        if (
                            key in stable_native_map
                            and parse_boolean(
                                stable_native_map[key].get(
                                    "native_top_ranked_set_changed", False
                                )
                            )
                        )
                    ]
                    entrant_values = [
                        str(row["native_top_ranked_change_reason"])
                        in {"ENTRANT", "EXIT", "MIXED"}
                        for (row) in (native_changed_rows)
                    ]
                    mapping_definition = (
                        "Primary panels removed by the strict entity-stability control"
                    )
                    entrant_definition = (
                        "Changed stable-entity release-native top-ranked panels"
                    )
                else:
                    primary_overlap = [
                        row
                        for (row) in (rows)
                        if ((pair, str(row["canonical_disease_id"])) in primary_keys)
                    ]
                    primary_index = {
                        (pair, str(row["canonical_disease_id"])): row
                        for (row) in (primary_pair)
                    }
                    mapping_values = [
                        str(
                            primary_index[pair, str(row["canonical_disease_id"])][
                                "mapping_stability"
                            ]
                        )
                        != "NO_MAPPING_CONTEXT_FLAG"
                        for (row) in (primary_overlap)
                    ]
                    entrant_values = []
                    mapping_definition = "Ranking sensitivity mapping context among evidence-domain panels overlapping primary support"
                    entrant_definition = (
                        "Not estimable on fixed evidence-domain support"
                    )
            top10_n, top10_d, top10_value = proportion(top10)
            top_n, top_d, top_value = proportion(top_ranked)
            mapping_n, mapping_d, mapping_value = proportion(mapping_values)
            entrant_n, entrant_d, entrant_value = proportion(entrant_values)
            diseases = {str(row["canonical_disease_id"]) for (row) in (rows)}
            areas = (
                set().union(
                    *(area_memberships.get(disease, set()) for (disease) in (diseases))
                )
                if (diseases)
                else set()
            )
            output.append(
                {
                    "release_pair": pair,
                    "analysis_layer": layer,
                    "support_denominator": layer,
                    "denominator_definition": {
                        "PRIMARY_OVERALL_FIXED": "Ranking sensitivity primary pair-specific persistent overall-score panels",
                        "STABLE_ENTITY_FIXED": "Independently rebuilt fixed panels restricted to frozen stable diseases and targets",
                        "GENETIC_ASSOCIATION_FIXED": "Independent fixed support for the release-provided genetic-association score",
                    }[layer],
                    "eligible_disease_panel_count": len(rows),
                    "eligible_target_pair_count": int(
                        score.get("eligible_target_pair_count", 0) or 0
                    ),
                    "median_absolute_score_change": score.get("median_absolute_change"),
                    "q90_absolute_score_change": score.get("q90_absolute_change"),
                    "median_kendall_tau_b": linear_quantile(tau_values, 0.5),
                    "top_10_change_numerator": top10_n,
                    "top_10_change_denominator": top10_d,
                    "top_10_change_fraction": top10_value,
                    "top_ranked_change_numerator": top_n,
                    "top_ranked_change_denominator": top_d,
                    "top_ranked_change_fraction": top_value,
                    "mapping_affected_numerator": mapping_n,
                    "mapping_affected_denominator": mapping_d,
                    "mapping_affected_fraction": mapping_value,
                    "mapping_affected_definition": mapping_definition,
                    "entrant_driven_numerator": entrant_n if (entrant_d) else None,
                    "entrant_driven_denominator": entrant_d if (entrant_d) else None,
                    "entrant_driven_fraction": entrant_value,
                    "entrant_driven_definition": entrant_definition,
                    "therapeutic_area_coverage": len(areas),
                    "distinct_support_denominator": True,
                }
            )
    return output

def build_persistence_rows(
    primary_rows: Sequence[Mapping[str, Any]],
    stable_rows: Sequence[Mapping[str, Any]],
    domain_rows: Sequence[Mapping[str, Any]],
    diagnostic_rows: Sequence[Mapping[str, Any]],
    stable_registry_rows: Sequence[Mapping[str, Any]],
    area_memberships: Mapping[str, set[str]],
) -> list[dict[str, Any]]:
    primary = {
        (str(row["release_pair"]), str(row["canonical_disease_id"])): row
        for (row) in (primary_rows)
    }
    stable = {
        (str(row["release_pair"]), str(row["canonical_disease_id"])): row
        for (row) in (stable_rows)
    }
    domain = {
        (str(row["release_pair"]), str(row["canonical_disease_id"])): row
        for (row) in (domain_rows)
    }
    registry = {
        (str(row["release_pair"]), str(row["canonical_disease_id"])): row
        for (row) in (stable_registry_rows)
    }
    diagnostic_keys = {
        (str(row["release_pair"]), str(row["canonical_disease_id"]))
        for (row) in (diagnostic_rows)
        if (int(row.get("target_count_with_diagnostic_support") or 0) > 0)
    }
    output: list[dict[str, Any]] = []
    for (key) in (sorted(set(primary) | set(stable) | set(domain))):
        primary_row = primary.get(key)
        stable_row = stable.get(key)
        domain_row = domain.get(key)
        primary_changed = (
            None
            if (primary_row is None)
            else parse_boolean(primary_row["fixed_top_ranked_set_changed"])
        )
        stable_changed = (
            None
            if (stable_row is None)
            else parse_boolean(stable_row["top_ranked_set_changed"])
        )
        domain_changed = (
            None
            if (domain_row is None)
            else parse_boolean(domain_row["top_ranked_set_changed"])
        )
        if (not any(
            (
                bool(value)
                for (value) in ([primary_changed, stable_changed, domain_changed])
            )
        )):
            continue
        source = primary_row or stable_row or domain_row
        registry_row = registry.get(key, {})
        primary_count = int(
            registry_row.get("fixed_target_count")
            or (primary_row.get("fixed_target_count") if (primary_row) else 0)
            or 0
        )
        stable_count = int(
            registry_row.get("stable_fixed_target_count")
            or (stable_row.get("fixed_target_count") if (stable_row) else 0)
            or 0
        )
        removed_fraction = (
            None
            if (primary_count == 0)
            else (primary_count - stable_count) / primary_count
        )
        classification = classify_persistence(
            primary_changed,
            stable_row is not None,
            stable_changed,
            domain_row is not None,
            domain_changed,
        )
        output.append(
            {
                "release_pair": key[0],
                "canonical_disease_id": key[1],
                "canonical_label": source.get("canonical_label", ""),
                "canonical_therapeutic_area_ids": source.get(
                    "canonical_therapeutic_area_ids", ""
                ),
                "stable_therapeutic_area_ids": joined_identifiers(
                    area_memberships.get(key[1], set())
                ),
                "overall_primary_eligible": primary_row is not None,
                "stable_entity_eligible": stable_row is not None,
                "evidence_domain_eligible": domain_row is not None,
                "overall_top_ranked_set_changed": primary_changed,
                "stable_entity_top_ranked_set_changed": stable_changed,
                "evidence_domain_top_ranked_set_changed": domain_changed,
                "overall_kendall_tau_b": None
                if (primary_row is None)
                else primary_row["kendall_tau_b"],
                "stable_entity_kendall_tau_b": None
                if (stable_row is None)
                else stable_row["kendall_tau_b"],
                "evidence_domain_kendall_tau_b": None
                if (domain_row is None)
                else domain_row["kendall_tau_b"],
                "evidence_domain_corroborates_overall": None
                if (domain_row is None)
                else domain_row["corroborates_overall_top_change"],
                "mapping_context_present": False
                if (primary_row is None)
                else str(primary_row["mapping_stability"]) != "NO_MAPPING_CONTEXT_FLAG",
                "primary_fixed_target_count": primary_count or None,
                "stable_fixed_target_count": stable_count or None,
                "entity_stability_control_removed_fraction": removed_fraction,
                "diagnostic_datasource_context_available": key in diagnostic_keys,
                "persistence_class": classification,
                "descriptive_not_causal": True,
            }
        )
    return output
