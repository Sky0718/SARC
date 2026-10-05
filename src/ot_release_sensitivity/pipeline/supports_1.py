from __future__ import annotations
from types import SimpleNamespace
from pathlib import Path
import re
import shutil
import duckdb
import pyarrow.parquet as pq
from ot_release_sensitivity.pipeline.harmonisation import discover_inputs
from ot_release_sensitivity.pipeline.harmonisation import projection_sql
from ot_release_sensitivity.pipeline.harmonisation import quote_literal
from ot_release_sensitivity.pipeline.harmonisation import required_assertion_sql

def validate_memory_limit(value):
    if (not re.fullmatch("[1-9][0-9]*(KB|MB|GB|TB)", str(value).upper())):
        raise ValueError("Invalid DuckDB memory limit")
    return str(value).upper()

def sql_path(path):
    return quote_literal(str(Path(path).resolve()).replace("\\", "/"))

def parquet_relation(path):
    return f"read_parquet({sql_path(path)})"

def csv_relation(path):
    return f"read_csv({sql_path(path)}, header = true, all_varchar = true)"

def copy_query(connection, query, destination, file_format):
    path = Path(destination)
    path.parent.mkdir(parents = True, exist_ok = True)
    partial = path.with_suffix(path.suffix + ".partial")
    partial.unlink(missing_ok = True)
    if (file_format == "PARQUET"):
        connection.execute(
            f"COPY ({query}) TO {sql_path(partial)} (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000)"
        )
    else:
        connection.execute(
            f"""COPY ({query}) TO {sql_path(partial)} (FORMAT CSV, HEADER TRUE, DELIMITER ',', QUOTE '"', ESCAPE '"')"""
        )
    partial.replace(path)

def cleanup_partial_outputs(configuration, *, root):
    for (relative) in (configuration["outputs"].values()):
        path = root / relative
        path.with_suffix(path.suffix + ".partial").unlink(missing_ok = True)

def cleanup_stage_outputs(temp_directory, *, root):
    root = Path(temp_directory).resolve()
    interim_root = (root / "data" / "interim").resolve()
    if (not root.is_relative_to(interim_root)):
        raise RuntimeError(
            "Support construction stage cleanup must remain within data/interim"
        )
    for (path) in (root.glob("support_construction_*.parquet*")):
        if (path.is_file() and path.resolve().is_relative_to(root)):
            path.unlink(missing_ok = True)

def configure_resources(connection, configuration, args, *, root):
    minimum_free = int(configuration["resources"]["minimum_free_gib"]) * 1024**3
    available = shutil.disk_usage(root).free
    if (available < minimum_free):
        raise RuntimeError(
            f"Insufficient free disk space for Support construction: {available} bytes available, {minimum_free} required"
        )
    temp_directory = (
        root / configuration["resources"]["duckdb_temp_directory"]
    ).resolve()
    interim_root = (root / "data" / "interim").resolve()
    if (not temp_directory.is_relative_to(interim_root)):
        raise RuntimeError("DuckDB temporary directory must remain within data/interim")
    temp_directory.mkdir(parents = True, exist_ok = True)
    temp_sql = str(temp_directory).replace("\\", "/")
    connection.execute(f"SET temp_directory = {quote_literal(temp_sql)}")
    connection.execute(f"SET threads = {int(args.threads)}")
    connection.execute(
        f"SET memory_limit = {quote_literal(validate_memory_limit(args.memory_limit))}"
    )
    connection.execute("SET preserve_insertion_order = false")
    return temp_directory

def partition_expression(disease_column, target_column, partition_count):
    return f"CAST(hash({disease_column}, {target_column}) % {int(partition_count)} AS INTEGER)"

def canonical_association_query(
    configuration, inputs, release, build_partition, partition_count
):
    parts = []
    dataset = "association_overall_direct"
    contract = configuration["dataset_fields"][dataset][release]
    parts.append(projection_sql(inputs[release, dataset], release, dataset, contract))
    raw = " UNION ALL ".join(parts)
    lower, upper = configuration["association_support"]["score_bounds"]
    partition = partition_expression(
        "COALESCE(canonical_disease_id, source_disease_id, '')",
        "COALESCE(canonical_target_id, source_target_id, '')",
        partition_count,
    )
    return f"\nWITH raw AS ({raw}),\nmapped AS (\n    SELECT\n        r.*,\n        d.canonical_disease_id,\n        t.canonical_target_id,\n        d.identifier_mapping_status AS disease_mapping_status,\n        d.mapping_status AS disease_crosswalk_status,\n        d.exclusion_reason AS disease_mapping_detail,\n        t.mapping_status AS target_mapping_status,\n        t.exclusion_reason AS target_mapping_detail,\n        d.canonical_therapeutic_area_ids,\n        COALESCE(d.primary_support_eligible = 'true', false) AS disease_primary_mapping_eligible,\n        COALESCE(t.primary_support_eligible = 'true', false) AS target_primary_mapping_eligible,\n        COALESCE(d.stable_entity_eligible = 'true', false) AND COALESCE(t.stable_entity_eligible = 'true', false) AS stable_entity_eligible,\n        COALESCE(d.primary_support_eligible = 'true', false) AND COALESCE(t.primary_support_eligible = 'true', false) AND COALESCE(d.canonical_disease_id, '') <> '' AND COALESCE(t.canonical_target_id, '') <> '' AS mapping_eligible\n    FROM raw r\n    LEFT JOIN support_construction_disease_crosswalk d\n      ON r.release = d.release AND r.source_disease_id = d.source_disease_id\n    LEFT JOIN support_construction_target_crosswalk t\n      ON r.release = t.release AND r.source_target_id = t.source_target_id\n),\npartitioned AS (\n    SELECT *, {partition} AS build_partition\n    FROM mapped\n),\nclassified AS (\n    SELECT *,\n        CASE\n            WHEN score IS NULL THEN 'PRESENT_NULL'\n            WHEN NOT isfinite(score) OR score < {lower} OR score > {upper} THEN 'SCHEMA_INCOMPARABLE'\n            WHEN score = 0 THEN 'PRESENT_ZERO'\n            ELSE 'PRESENT_POSITIVE'\n        END AS association_state\n    FROM partitioned\n    WHERE build_partition = {int(build_partition)}\n)\nSELECT\n    release,\n    source_disease_id,\n    source_target_id,\n    canonical_disease_id,\n    canonical_target_id,\n    score,\n    evidence_count,\n    disease_mapping_status,\n    disease_crosswalk_status,\n    disease_mapping_detail,\n    target_mapping_status,\n    target_mapping_detail,\n    association_state,\n    disease_primary_mapping_eligible,\n    target_primary_mapping_eligible,\n    mapping_eligible,\n    mapping_eligible AND association_state = 'PRESENT_POSITIVE' AND evidence_count IS NOT NULL AND evidence_count >= 1 AS primary_association_eligible,\n    stable_entity_eligible,\n    canonical_therapeutic_area_ids,\n    build_partition,\n    score_source_dataset,\n    source_object_relative_path,\n    release || '|' || source_disease_id || '|' || source_target_id AS source_row_locator,\n    CASE\n        WHEN disease_primary_mapping_eligible IS NOT TRUE THEN COALESCE(NULLIF(disease_crosswalk_status, ''), 'DISEASE_MAPPING_INELIGIBLE')\n        WHEN target_primary_mapping_eligible IS NOT TRUE THEN COALESCE(NULLIF(target_mapping_status, ''), 'TARGET_IDENTIFIER_NOT_EXACT')\n        WHEN canonical_disease_id IS NULL OR canonical_disease_id = '' THEN 'DISEASE_MAPPING_INELIGIBLE'\n        WHEN canonical_target_id IS NULL OR canonical_target_id = '' THEN 'TARGET_IDENTIFIER_NOT_EXACT'\n        WHEN association_state = 'PRESENT_NULL' THEN 'PRESENT_NULL'\n        WHEN association_state = 'SCHEMA_INCOMPARABLE' AND score IS NOT NULL AND NOT isfinite(score) THEN 'SCORE_NONFINITE'\n        WHEN association_state = 'SCHEMA_INCOMPARABLE' THEN 'SCORE_OUT_OF_BOUNDS'\n        WHEN association_state = 'PRESENT_ZERO' THEN 'PRESENT_ZERO'\n        WHEN evidence_count IS NULL OR evidence_count < 1 THEN 'MISSING_OR_NONPOSITIVE_EVIDENCE_COUNT'\n        ELSE ''\n    END AS exclusion_reason\nFROM classified\nORDER BY source_disease_id, source_target_id, source_object_relative_path\n"

def merge_parquet_files(sources, destination, row_group_size):
    path = Path(destination)
    path.parent.mkdir(parents = True, exist_ok = True)
    partial = path.with_suffix(path.suffix + ".partial")
    partial.unlink(missing_ok = True)
    writer = None
    expected_schema = None
    try:
        for (source) in (sources):
            parquet_file = pq.ParquetFile(source)
            if (writer is None):
                expected_schema = parquet_file.schema_arrow
                writer = pq.ParquetWriter(partial, expected_schema, compression = "zstd")
            elif (not parquet_file.schema_arrow.equals(
                expected_schema, check_metadata = False
            )):
                raise RuntimeError("Partitioned Parquet stage schemas differ")
            for (batch) in (parquet_file.iter_batches(batch_size = row_group_size)):
                writer.write_batch(batch, row_group_size = row_group_size)
        if (writer is None):
            raise RuntimeError("No partitioned Parquet stages were available")
        writer.close()
        writer = None
        partial.replace(path)
    except Exception:
        if (writer is not None):
            try:
                writer.close()
            except Exception:
                pass
        partial.unlink(missing_ok = True)
        raise

def build_partitioned_parquet(
    connection,
    labelled_queries,
    destination,
    temp_directory,
    stem,
    row_group_size,
    *,
    root,
):
    stages = []
    temp_root = Path(temp_directory).resolve()
    interim_root = (root / "data" / "interim").resolve()
    if (not temp_root.is_relative_to(interim_root)):
        raise RuntimeError(
            "Support construction partition stages must remain within data/interim"
        )
    try:
        for (index, (label, query)) in (enumerate(labelled_queries)):
            safe_label = re.sub("[^A-Za-z0-9_]+", "_", label).strip("_")
            stage = temp_root / f"{stem}_{index:04d}_{safe_label}.parquet"
            if (not stage.resolve().is_relative_to(temp_root)):
                raise RuntimeError(
                    "Support construction stage path escaped the configured temporary directory"
                )
            stage.unlink(missing_ok = True)
            stage.with_suffix(stage.suffix + ".partial").unlink(missing_ok = True)
            stages.append(stage)
            copy_query(connection, query, stage, "PARQUET")
        merge_parquet_files(stages, destination, row_group_size)
    finally:
        for (stage) in (stages):
            stage.unlink(missing_ok = True)
            stage.with_suffix(stage.suffix + ".partial").unlink(missing_ok = True)

def state_transition_query(canonical_path, old_release, new_release, build_partition):
    relation = parquet_relation(canonical_path)
    pair = f"{old_release}_to_{new_release}"
    return f"\nWITH old_rows AS (\n    SELECT * FROM {relation}\n    WHERE release = {quote_literal(old_release)} AND build_partition = {int(build_partition)} AND mapping_eligible\n),\nnew_rows AS (\n    SELECT * FROM {relation}\n    WHERE release = {quote_literal(new_release)} AND build_partition = {int(build_partition)} AND mapping_eligible\n),\njoined AS (\n    SELECT\n        {quote_literal(pair)} AS release_pair,\n        {quote_literal(old_release)} AS old_release,\n        {quote_literal(new_release)} AS new_release,\n        {int(build_partition)}::INTEGER AS build_partition,\n        COALESCE(o.canonical_disease_id, n.canonical_disease_id) AS canonical_disease_id,\n        COALESCE(o.canonical_target_id, n.canonical_target_id) AS canonical_target_id,\n        o.source_disease_id AS old_source_disease_id,\n        n.source_disease_id AS new_source_disease_id,\n        o.source_target_id AS old_source_target_id,\n        n.source_target_id AS new_source_target_id,\n        o.score AS old_score,\n        n.score AS new_score,\n        o.evidence_count AS old_evidence_count,\n        n.evidence_count AS new_evidence_count,\n        COALESCE(o.association_state, 'ABSENT_PAIR') AS old_state,\n        COALESCE(n.association_state, 'ABSENT_PAIR') AS new_state,\n        COALESCE(o.primary_association_eligible, false) AS old_primary_eligible,\n        COALESCE(n.primary_association_eligible, false) AS new_primary_eligible,\n        o.source_disease_id IS NOT NULL AND n.source_disease_id IS NOT NULL AND (o.source_disease_id <> n.source_disease_id OR o.source_target_id <> n.source_target_id) AS raw_identifier_mapping_affected,\n        COALESCE(o.stable_entity_eligible, false) AND COALESCE(n.stable_entity_eligible, false) AS stable_entity_eligible\n    FROM old_rows o\n    FULL OUTER JOIN new_rows n\n      USING (canonical_disease_id, canonical_target_id)\n)\nSELECT *,\n    CASE\n        WHEN old_state = 'SCHEMA_INCOMPARABLE' OR new_state = 'SCHEMA_INCOMPARABLE' THEN 'SCHEMA_INCOMPARABLE'\n        WHEN old_primary_eligible AND new_primary_eligible AND raw_identifier_mapping_affected THEN 'PERSISTENT_MAPPING_AFFECTED'\n        WHEN old_primary_eligible AND new_primary_eligible THEN 'PERSISTENT_POSITIVE'\n        WHEN old_state = 'ABSENT_PAIR' AND new_primary_eligible THEN 'ENTRANT'\n        WHEN old_primary_eligible AND new_state = 'ABSENT_PAIR' THEN 'EXIT'\n        WHEN raw_identifier_mapping_affected THEN 'MAPPING_AFFECTED'\n        ELSE 'NON_PRIMARY_STATE_TRANSITION'\n    END AS transition_class\nFROM joined\nORDER BY canonical_disease_id, canonical_target_id\n"

def three_release_query(canonical_path, build_partition):
    relation = parquet_relation(canonical_path)
    return f"\nWITH r25 AS (SELECT * FROM {relation} WHERE release = '25.12' AND build_partition = {int(build_partition)} AND primary_association_eligible),\nr03 AS (SELECT * FROM {relation} WHERE release = '26.03' AND build_partition = {int(build_partition)} AND primary_association_eligible),\nr06 AS (SELECT * FROM {relation} WHERE release = '26.06' AND build_partition = {int(build_partition)} AND primary_association_eligible)\nSELECT\n    {int(build_partition)}::INTEGER AS build_partition,\n    r03.canonical_disease_id,\n    r03.canonical_target_id,\n    r25.score AS score_25_12,\n    r03.score AS score_26_03,\n    r06.score AS score_26_06,\n    r25.evidence_count AS evidence_count_25_12,\n    r03.evidence_count AS evidence_count_26_03,\n    r06.evidence_count AS evidence_count_26_06,\n    r25.stable_entity_eligible AND r03.stable_entity_eligible AND r06.stable_entity_eligible AS stable_entity_eligible,\n    r03.canonical_therapeutic_area_ids\nFROM r25\nJOIN r03 USING (canonical_disease_id, canonical_target_id)\nJOIN r06 USING (canonical_disease_id, canonical_target_id)\nORDER BY canonical_disease_id, canonical_target_id\n"

def fixed_members_query(
    canonical_path, triple_path, support_name, old_release, new_release, build_partition
):
    relation = parquet_relation(canonical_path)
    triple = parquet_relation(triple_path)
    pair = f"{old_release}_to_{new_release}"
    if (support_name == "PAIR_SPECIFIC_PERSISTENT"):
        return f"\nSELECT\n    'PAIR_SPECIFIC_PERSISTENT' AS support_name,\n    {quote_literal(pair)} AS release_pair,\n    {quote_literal(old_release)} AS old_release,\n    {quote_literal(new_release)} AS new_release,\n    {int(build_partition)}::INTEGER AS build_partition,\n    o.canonical_disease_id,\n    o.canonical_target_id,\n    o.score AS old_score,\n    n.score AS new_score,\n    o.evidence_count AS old_evidence_count,\n    n.evidence_count AS new_evidence_count,\n    o.stable_entity_eligible AND n.stable_entity_eligible AS stable_entity_member,\n    o.canonical_therapeutic_area_ids\nFROM {relation} o\nJOIN {relation} n USING (canonical_disease_id, canonical_target_id)\nWHERE o.release = {quote_literal(old_release)} AND n.release = {quote_literal(new_release)}\n  AND o.build_partition = {int(build_partition)} AND n.build_partition = {int(build_partition)}\n  AND o.primary_association_eligible AND n.primary_association_eligible\nORDER BY o.canonical_disease_id, o.canonical_target_id\n"
    if (support_name != "THREE_RELEASE_PERSISTENT"):
        raise ValueError(f"Unknown fixed support: {support_name}")
    old_suffix = old_release.replace(".", "_")
    new_suffix = new_release.replace(".", "_")
    return f"\nSELECT\n    'THREE_RELEASE_PERSISTENT' AS support_name,\n    {quote_literal(pair)} AS release_pair,\n    {quote_literal(old_release)} AS old_release,\n    {quote_literal(new_release)} AS new_release,\n    {int(build_partition)}::INTEGER AS build_partition,\n    canonical_disease_id,\n    canonical_target_id,\n    score_{old_suffix} AS old_score,\n    score_{new_suffix} AS new_score,\n    evidence_count_{old_suffix} AS old_evidence_count,\n    evidence_count_{new_suffix} AS new_evidence_count,\n    stable_entity_eligible AS stable_entity_member,\n    canonical_therapeutic_area_ids\nFROM {triple}\nWHERE build_partition = {int(build_partition)}\nORDER BY canonical_disease_id, canonical_target_id\n"

def fixed_registry_query(configuration, members_path, disease_crosswalk_path):
    members = parquet_relation(members_path)
    diseases = csv_relation(disease_crosswalk_path)
    lower = configuration["ranking_panels"]["lower_sensitivity_minimum_targets"]
    primary = configuration["ranking_panels"]["primary_minimum_targets"]
    higher = configuration["ranking_panels"]["higher_sensitivity_minimum_targets"]
    distinct = configuration["ranking_panels"]["minimum_distinct_scores_per_release"]
    return f"\nWITH labels AS (\n    SELECT canonical_disease_id, source_label AS canonical_label\n    FROM {diseases}\n    WHERE release = '26.03' AND canonical_disease_id <> ''\n),\npanels AS (\n    SELECT\n        support_name,\n        release_pair,\n        canonical_disease_id,\n        any_value(canonical_therapeutic_area_ids) AS canonical_therapeutic_area_ids,\n        COUNT(*) AS fixed_target_count,\n        COUNT(DISTINCT old_score) AS old_distinct_score_count,\n        COUNT(DISTINCT new_score) AS new_distinct_score_count,\n        COUNT(*) FILTER (WHERE stable_entity_member) AS stable_fixed_target_count,\n        COUNT(DISTINCT CASE WHEN stable_entity_member THEN old_score END) AS stable_old_distinct_score_count,\n        COUNT(DISTINCT CASE WHEN stable_entity_member THEN new_score END) AS stable_new_distinct_score_count\n    FROM {members}\n    GROUP BY support_name, release_pair, canonical_disease_id\n)\nSELECT\n    p.support_name,\n    p.release_pair,\n    p.canonical_disease_id,\n    l.canonical_label,\n    p.canonical_therapeutic_area_ids,\n    p.fixed_target_count,\n    p.old_distinct_score_count,\n    p.new_distinct_score_count,\n    p.fixed_target_count >= {lower} AND p.old_distinct_score_count >= {distinct} AND p.new_distinct_score_count >= {distinct} AS eligible_min_20,\n    p.fixed_target_count >= {primary} AND p.old_distinct_score_count >= {distinct} AND p.new_distinct_score_count >= {distinct} AS eligible_primary_30,\n    p.fixed_target_count >= {higher} AND p.old_distinct_score_count >= {distinct} AND p.new_distinct_score_count >= {distinct} AS eligible_min_50,\n    p.support_name = 'PAIR_SPECIFIC_PERSISTENT' AND p.fixed_target_count >= {primary} AND p.old_distinct_score_count >= {distinct} AND p.new_distinct_score_count >= {distinct} AS primary_panel_eligible,\n    p.stable_fixed_target_count,\n    p.stable_old_distinct_score_count,\n    p.stable_new_distinct_score_count,\n    p.stable_fixed_target_count >= {primary} AND p.stable_old_distinct_score_count >= {distinct} AND p.stable_new_distinct_score_count >= {distinct} AS stable_panel_eligible,\n    CASE\n        WHEN p.fixed_target_count < {primary} THEN 'FEWER_THAN_MINIMUM_FIXED_TARGETS'\n        WHEN p.old_distinct_score_count < {distinct} THEN 'DEGENERATE_OLD_SCORE_DISTRIBUTION'\n        WHEN p.new_distinct_score_count < {distinct} THEN 'DEGENERATE_NEW_SCORE_DISTRIBUTION'\n        ELSE ''\n    END AS exclusion_reason\nFROM panels p\nLEFT JOIN labels l USING (canonical_disease_id)\nORDER BY support_name, release_pair, canonical_disease_id\n"

def native_members_query(canonical_path, release, build_partition):
    relation = parquet_relation(canonical_path)
    return f"\nSELECT\n    release,\n    canonical_disease_id,\n    canonical_target_id,\n    score,\n    evidence_count,\n    build_partition,\n    stable_entity_eligible AS stable_entity_member,\n    canonical_therapeutic_area_ids\nFROM {relation}\nWHERE release = {quote_literal(release)} AND build_partition = {int(build_partition)} AND primary_association_eligible\nORDER BY canonical_disease_id, canonical_target_id\n"

def native_registry_query(configuration, members_path, disease_crosswalk_path):
    members = parquet_relation(members_path)
    diseases = csv_relation(disease_crosswalk_path)
    return f"\nWITH labels AS (\n    SELECT canonical_disease_id, source_label AS canonical_label\n    FROM {diseases}\n    WHERE release = '26.03' AND canonical_disease_id <> ''\n)\nSELECT\n    m.release,\n    m.canonical_disease_id,\n    l.canonical_label,\n    any_value(m.canonical_therapeutic_area_ids) AS canonical_therapeutic_area_ids,\n    COUNT(*) AS native_target_count,\n    COUNT(DISTINCT m.score) AS distinct_score_count,\n    COUNT(*) FILTER (WHERE m.stable_entity_member) AS stable_native_target_count,\n    COUNT(*) > 0 AS release_native_panel_eligible,\n    COUNT(*) >= 5 AS top_5_estimable,\n    COUNT(*) >= 10 AS top_10_estimable,\n    COUNT(*) >= 20 AS top_20_estimable\nFROM {members} m\nLEFT JOIN labels l USING (canonical_disease_id)\nGROUP BY m.release, m.canonical_disease_id, l.canonical_label\nORDER BY CASE m.release WHEN '25.12' THEN 1 WHEN '26.03' THEN 2 WHEN '26.06' THEN 3 END, m.canonical_disease_id\n"

def evidence_domain_query(
    configuration,
    inputs,
    endpoint_level,
    dataset,
    endpoint_id,
    release,
    build_partition,
    partition_count,
):
    contract = configuration["dataset_fields"][dataset][release]
    projected = projection_sql(
        inputs[release, dataset], release, dataset, contract, endpoint_id
    )
    raw = f"SELECT {quote_literal(endpoint_level)} AS endpoint_level, {quote_literal(endpoint_id)} AS endpoint_id, * FROM ({projected})"
    lower, upper = configuration["association_support"]["score_bounds"]
    partition = partition_expression(
        "COALESCE(d.canonical_disease_id, r.source_disease_id, '')",
        "COALESCE(t.canonical_target_id, r.source_target_id, '')",
        partition_count,
    )
    return f"\nWITH raw AS ({raw}),\nmapped AS (\nSELECT\n    r.endpoint_level,\n    r.endpoint_id,\n    r.release,\n    r.source_disease_id,\n    r.source_target_id,\n    d.canonical_disease_id,\n    t.canonical_target_id,\n    r.score,\n    r.evidence_count,\n    {partition} AS build_partition,\n    CASE\n        WHEN r.score IS NULL THEN 'PRESENT_NULL'\n        WHEN NOT isfinite(r.score) OR r.score < {lower} OR r.score > {upper} THEN 'SCHEMA_INCOMPARABLE'\n        WHEN r.score = 0 THEN 'PRESENT_ZERO'\n        ELSE 'PRESENT_POSITIVE'\n    END AS association_state,\n    d.primary_support_eligible = 'true' AND t.primary_support_eligible = 'true' AND r.score IS NOT NULL AND isfinite(r.score) AND r.score BETWEEN {lower} AND {upper} AS domain_support_eligible,\n    d.stable_entity_eligible = 'true' AND t.stable_entity_eligible = 'true' AS stable_entity_eligible,\n    d.canonical_therapeutic_area_ids,\n    r.score_source_dataset,\n    r.source_object_relative_path\nFROM raw r\nLEFT JOIN support_construction_disease_crosswalk d\n  ON r.release = d.release AND r.source_disease_id = d.source_disease_id\nLEFT JOIN support_construction_target_crosswalk t\n  ON r.release = t.release AND r.source_target_id = t.source_target_id\n)\nSELECT * FROM mapped\nWHERE build_partition = {int(build_partition)}\nORDER BY canonical_disease_id, canonical_target_id, source_disease_id, source_target_id, source_object_relative_path\n"

def evidence_registry_query(configuration, evidence_path):
    evidence = parquet_relation(evidence_path)
    primary = configuration["ranking_panels"]["primary_minimum_targets"]
    distinct = configuration["ranking_panels"]["minimum_distinct_scores_per_release"]
    parts = []
    for (old_release, new_release) in (configuration["release_pairs"]):
        pair = f"{old_release}_to_{new_release}"
        parts.append(
            f"\nSELECT\n    o.endpoint_level,\n    o.endpoint_id,\n    {quote_literal(pair)} AS release_pair,\n    o.canonical_disease_id,\n    COUNT(*) AS fixed_target_count,\n    COUNT(DISTINCT o.score) AS old_distinct_score_count,\n    COUNT(DISTINCT n.score) AS new_distinct_score_count,\n    COUNT(*) >= {primary} AND COUNT(DISTINCT o.score) >= {distinct} AND COUNT(DISTINCT n.score) >= {distinct} AS panel_eligible,\n    any_value(o.canonical_therapeutic_area_ids) AS canonical_therapeutic_area_ids\nFROM {evidence} o\nJOIN {evidence} n\n  ON o.endpoint_level = n.endpoint_level\n AND o.endpoint_id = n.endpoint_id\n AND o.canonical_disease_id = n.canonical_disease_id\n AND o.canonical_target_id = n.canonical_target_id\nWHERE o.release = {quote_literal(old_release)} AND n.release = {quote_literal(new_release)}\n  AND o.domain_support_eligible AND n.domain_support_eligible\nGROUP BY o.endpoint_level, o.endpoint_id, o.canonical_disease_id\n"
        )
    return (
        " UNION ALL ".join(parts)
        + " ORDER BY endpoint_level, endpoint_id, release_pair, canonical_disease_id"
    )

def entrant_exit_query(transition_path, release_pair, build_partition):
    transitions = parquet_relation(transition_path)
    return f"\nSELECT *,\n    CASE\n        WHEN transition_class = 'ENTRANT' THEN 'ENTRANT'\n        WHEN transition_class = 'EXIT' THEN 'EXIT'\n        WHEN transition_class IN ('PERSISTENT_MAPPING_AFFECTED', 'MAPPING_AFFECTED') THEN 'MAPPING_AFFECTED'\n        ELSE 'NON_PRIMARY_STATE'\n    END AS record_category,\n    old_primary_eligible AS old_release_native_member,\n    new_primary_eligible AS new_release_native_member,\n    'release_native_panel_members.parquet' AS rank_context_source\nFROM {transitions}\nWHERE release_pair = {quote_literal(release_pair)}\n  AND build_partition = {int(build_partition)}\n  AND transition_class <> 'PERSISTENT_POSITIVE'\nORDER BY canonical_disease_id, canonical_target_id\n"

def exclusion_registry_queries(
    configuration, canonical_path, fixed_registry_path, partition_count, *, root
):
    canonical = parquet_relation(canonical_path)
    diseases = csv_relation(root / configuration["outputs"]["disease_crosswalk"])
    targets = csv_relation(root / configuration["outputs"]["target_crosswalk"])
    fixed = csv_relation(fixed_registry_path)
    queries = []
    for (release) in (configuration["releases"]):
        queries.append(
            (
                f"disease_{release}",
                f"\nSELECT\n    'DISEASE_MAPPING' AS exclusion_scope,\n    'DISEASE_ENTITY' AS unit_type,\n    release,\n    '' AS release_pair,\n    -1::INTEGER AS build_partition,\n    release || '|' || source_disease_id AS unit_id,\n    CASE WHEN mapping_status = 'EXCLUDED_ENTITY_TYPE' THEN 'EXCLUDED_ENTITY_TYPE' ELSE identifier_mapping_status END AS primary_reason,\n    exclusion_reason AS detail\nFROM {diseases}\nWHERE release = {quote_literal(release)} AND primary_support_eligible <> 'true'\nORDER BY unit_id\n",
            )
        )
    for (release) in (configuration["releases"]):
        queries.append(
            (
                f"target_{release}",
                f"\nSELECT\n    'TARGET_MAPPING' AS exclusion_scope,\n    'TARGET_ENTITY' AS unit_type,\n    release,\n    '' AS release_pair,\n    -1::INTEGER AS build_partition,\n    release || '|' || source_target_id AS unit_id,\n    'TARGET_IDENTIFIER_NOT_EXACT' AS primary_reason,\n    exclusion_reason AS detail\nFROM {targets}\nWHERE release = {quote_literal(release)} AND primary_support_eligible <> 'true'\nORDER BY unit_id\n",
            )
        )
    for (release) in (configuration["releases"]):
        for (build_partition) in (range(partition_count)):
            queries.append(
                (
                    f"association_{release}_p{build_partition:02d}",
                    f"\nSELECT\n    'OVERALL_ASSOCIATION_SUPPORT' AS exclusion_scope,\n    'ASSOCIATION_ROW' AS unit_type,\n    release,\n    '' AS release_pair,\n    build_partition,\n    source_row_locator AS unit_id,\n    exclusion_reason AS primary_reason,\n    association_state AS detail\nFROM {canonical}\nWHERE release = {quote_literal(release)}\n  AND build_partition = {int(build_partition)}\n  AND primary_association_eligible IS NOT TRUE\nORDER BY unit_id\n",
                )
            )
    for (old_release, new_release) in (configuration["release_pairs"]):
        pair = f"{old_release}_to_{new_release}"
        queries.append(
            (
                f"fixed_{pair}",
                f"\nSELECT\n    'FIXED_PANEL_ELIGIBILITY' AS exclusion_scope,\n    'DISEASE_PANEL' AS unit_type,\n    '' AS release,\n    release_pair,\n    -1::INTEGER AS build_partition,\n    support_name || '|' || release_pair || '|' || canonical_disease_id AS unit_id,\n    exclusion_reason AS primary_reason,\n    'fixed_target_count=' || fixed_target_count AS detail\nFROM {fixed}\nWHERE release_pair = {quote_literal(pair)} AND support_name = 'PAIR_SPECIFIC_PERSISTENT' AND primary_panel_eligible <> 'true'\nORDER BY unit_id\n",
            )
        )
    return queries

def support_flow_query(
    configuration,
    canonical_path,
    transition_path,
    triple_path,
    fixed_registry_path,
    native_registry_path,
    evidence_registry_path,
    disease_crosswalk_path,
    target_crosswalk_path,
    therapeutic_area_path,
):
    canonical = parquet_relation(canonical_path)
    transitions = parquet_relation(transition_path)
    triple = parquet_relation(triple_path)
    fixed = csv_relation(fixed_registry_path)
    native = csv_relation(native_registry_path)
    evidence = csv_relation(evidence_registry_path)
    diseases = csv_relation(disease_crosswalk_path)
    targets = csv_relation(target_crosswalk_path)
    areas = csv_relation(therapeutic_area_path)
    return f"\nWITH association_release AS (\n    SELECT\n        release,\n        COUNT(*) AS raw_rows,\n        COUNT(DISTINCT (source_disease_id, source_target_id)) AS unique_source_pairs,\n        COUNT(*) FILTER (WHERE mapping_eligible) AS canonical_mapping_rows,\n        COUNT(*) FILTER (WHERE primary_association_eligible) AS primary_rows,\n        COUNT(DISTINCT source_disease_id) AS source_diseases,\n        COUNT(DISTINCT source_target_id) AS source_targets,\n        COUNT(DISTINCT canonical_disease_id) FILTER (WHERE mapping_eligible) AS canonical_diseases,\n        COUNT(DISTINCT canonical_target_id) FILTER (WHERE mapping_eligible) AS canonical_targets,\n        SUM(CAST(score AS DECIMAL(38,18))) AS raw_score_mass,\n        SUM(CAST(score AS DECIMAL(38,18))) FILTER (WHERE primary_association_eligible) AS primary_score_mass\n    FROM {canonical}\n    GROUP BY release\n),\ndisease_denominators AS (\n    SELECT release, COUNT(*) AS all_entities, COUNT(*) FILTER (WHERE release_native_eligible = 'true') AS eligible_entities\n    FROM {diseases}\n    GROUP BY release\n),\ntarget_denominators AS (\n    SELECT release, COUNT(*) AS all_entities\n    FROM {targets}\n    GROUP BY release\n),\ntransition_denominators AS (\n    SELECT release_pair, COUNT(*) AS union_pairs\n    FROM {transitions}\n    GROUP BY release_pair\n),\nfixed_denominators AS (\n    SELECT support_name, release_pair, COUNT(*) AS candidate_panels\n    FROM {fixed}\n    GROUP BY support_name, release_pair\n),\nnative_denominators AS (\n    SELECT release, COUNT(*) AS candidate_panels\n    FROM {native}\n    GROUP BY release\n),\nevidence_denominators AS (\n    SELECT endpoint_level, endpoint_id, release_pair, COUNT(*) AS candidate_panels\n    FROM {evidence}\n    GROUP BY endpoint_level, endpoint_id, release_pair\n),\nflow AS (\n    SELECT 'OVERALL_ASSOCIATION' AS scope, release, '' AS release_pair, 'RELEASE_NATIVE' AS support_name, 'ASSOCIATION_ROW' AS unit, 'RAW_ASSOCIATION_ROWS' AS stage, raw_rows AS count, NULL::BIGINT AS denominator, NULL::DOUBLE AS fraction, source_diseases AS disease_count, source_targets AS target_count, NULL::BIGINT AS therapeutic_area_count, raw_score_mass AS score_mass FROM association_release\n    UNION ALL\n    SELECT 'OVERALL_ASSOCIATION', release, '', 'RELEASE_NATIVE', 'SOURCE_DISEASE_TARGET_PAIR', 'UNIQUE_SOURCE_PAIRS', unique_source_pairs, raw_rows, unique_source_pairs::DOUBLE / NULLIF(raw_rows, 0), source_diseases, source_targets, NULL, raw_score_mass FROM association_release\n    UNION ALL\n    SELECT 'OVERALL_ASSOCIATION', release, '', 'CANONICAL', 'ASSOCIATION_ROW', 'CANONICAL_MAPPINGS_AVAILABLE', canonical_mapping_rows, raw_rows, canonical_mapping_rows::DOUBLE / NULLIF(raw_rows, 0), canonical_diseases, canonical_targets, NULL, NULL FROM association_release\n    UNION ALL\n    SELECT 'OVERALL_ASSOCIATION', release, '', 'RELEASE_NATIVE', 'ASSOCIATION_ROW', 'PRIMARY_ELIGIBLE', primary_rows, raw_rows, primary_rows::DOUBLE / NULLIF(raw_rows, 0), canonical_diseases, canonical_targets, NULL, primary_score_mass FROM association_release\n    UNION ALL\n    SELECT 'DISEASE_MAPPING', d.release, '', d.identifier_mapping_status, 'DISEASE_ENTITY', 'EXACT_OR_EXPLICIT_MAPPING', COUNT(*), den.eligible_entities, COUNT(*)::DOUBLE / NULLIF(den.eligible_entities, 0), COUNT(*), NULL, NULL, NULL\n    FROM {diseases} d JOIN disease_denominators den USING (release)\n    WHERE d.release_native_eligible = 'true' AND d.identifier_mapping_status IN ('EXACT_ID', 'EXPLICIT_REPLACEMENT', 'EXPLICIT_CROSS_REFERENCE', 'ONE_TO_ONE_CANONICAL')\n    GROUP BY d.release, d.identifier_mapping_status, den.eligible_entities\n    UNION ALL\n    SELECT 'DISEASE_MAPPING', d.release, '', d.identifier_mapping_status, 'DISEASE_ENTITY', 'AMBIGUOUS_OR_UNRESOLVED_MAPPING', COUNT(*), den.eligible_entities, COUNT(*)::DOUBLE / NULLIF(den.eligible_entities, 0), COUNT(*), NULL, NULL, NULL\n    FROM {diseases} d JOIN disease_denominators den USING (release)\n    WHERE d.release_native_eligible = 'true' AND d.identifier_mapping_status IN ('AMBIGUOUS_ONE_TO_MANY', 'AMBIGUOUS_MANY_TO_ONE', 'LABEL_ONLY_CANDIDATE', 'UNRESOLVED')\n    GROUP BY d.release, d.identifier_mapping_status, den.eligible_entities\n    UNION ALL\n    SELECT 'DISEASE_MAPPING', d.release, '', d.source_entity_type, 'DISEASE_ENTITY', 'EXCLUDED_ENTITY_TYPE', COUNT(*), den.all_entities, COUNT(*)::DOUBLE / NULLIF(den.all_entities, 0), COUNT(*), NULL, NULL, NULL\n    FROM {diseases} d JOIN disease_denominators den USING (release)\n    WHERE d.mapping_status = 'EXCLUDED_ENTITY_TYPE'\n    GROUP BY d.release, d.source_entity_type, den.all_entities\n    UNION ALL\n    SELECT 'TARGET_MAPPING', t.release, '', 'EXACT_ID', 'TARGET_ENTITY', 'EXACT_TARGET_MAPPINGS', COUNT(*), den.all_entities, COUNT(*)::DOUBLE / NULLIF(den.all_entities, 0), NULL, COUNT(*), NULL, NULL\n    FROM {targets} t JOIN target_denominators den USING (release)\n    WHERE t.mapping_status = 'EXACT_ID'\n    GROUP BY t.release, den.all_entities\n    UNION ALL\n    SELECT 'THERAPEUTIC_AREA', a.release, '', 'STABLE_MULTI_MEMBERSHIP', 'THERAPEUTIC_AREA_MEMBERSHIP', 'PRIMARY_STABLE_MEMBERSHIP', COUNT(*), NULL, NULL, COUNT(DISTINCT canonical_disease_id), NULL, COUNT(DISTINCT canonical_therapeutic_area_id), NULL\n    FROM {areas} a\n    WHERE primary_area_contrast_eligible = 'true'\n    GROUP BY a.release\n    UNION ALL\n    SELECT 'OVERALL_ASSOCIATION', '', t.release_pair, 'PAIR_SPECIFIC_PERSISTENT', 'CANONICAL_PAIR', 'PAIR_PERSISTENT', COUNT(*), den.union_pairs, COUNT(*)::DOUBLE / NULLIF(den.union_pairs, 0), COUNT(DISTINCT t.canonical_disease_id), COUNT(DISTINCT t.canonical_target_id), NULL, NULL\n    FROM {transitions} t JOIN transition_denominators den USING (release_pair)\n    WHERE t.transition_class IN ('PERSISTENT_POSITIVE', 'PERSISTENT_MAPPING_AFFECTED')\n    GROUP BY t.release_pair, den.union_pairs\n    UNION ALL\n    SELECT 'OVERALL_ASSOCIATION', '', t.release_pair, 'RELEASE_NATIVE', 'CANONICAL_PAIR', t.transition_class, COUNT(*), den.union_pairs, COUNT(*)::DOUBLE / NULLIF(den.union_pairs, 0), COUNT(DISTINCT t.canonical_disease_id), COUNT(DISTINCT t.canonical_target_id), NULL, NULL\n    FROM {transitions} t JOIN transition_denominators den USING (release_pair)\n    WHERE t.transition_class IN ('ENTRANT', 'EXIT', 'PERSISTENT_MAPPING_AFFECTED', 'MAPPING_AFFECTED')\n    GROUP BY t.release_pair, t.transition_class, den.union_pairs\n    UNION ALL\n    SELECT 'OVERALL_ASSOCIATION', '', 'three_release', 'THREE_RELEASE_PERSISTENT', 'CANONICAL_PAIR', 'THREE_RELEASE_PERSISTENT', COUNT(*), NULL, NULL, COUNT(DISTINCT canonical_disease_id), COUNT(DISTINCT canonical_target_id), NULL, NULL FROM {triple}\n    UNION ALL\n    SELECT 'FIXED_PANEL', '', f.release_pair, f.support_name, 'DISEASE_PANEL', 'ELIGIBLE_PRIMARY_30', COUNT(*), den.candidate_panels, COUNT(*)::DOUBLE / NULLIF(den.candidate_panels, 0), COUNT(DISTINCT f.canonical_disease_id), NULL, NULL, NULL\n    FROM {fixed} f JOIN fixed_denominators den USING (support_name, release_pair)\n    WHERE f.eligible_primary_30 = 'true'\n    GROUP BY f.release_pair, f.support_name, den.candidate_panels\n    UNION ALL\n    SELECT 'FIXED_PANEL', '', f.release_pair, f.support_name, 'DISEASE_PANEL', 'ELIGIBLE_STABLE_ENTITY', COUNT(*), den.candidate_panels, COUNT(*)::DOUBLE / NULLIF(den.candidate_panels, 0), COUNT(DISTINCT f.canonical_disease_id), NULL, NULL, NULL\n    FROM {fixed} f JOIN fixed_denominators den USING (support_name, release_pair)\n    WHERE f.stable_panel_eligible = 'true'\n    GROUP BY f.release_pair, f.support_name, den.candidate_panels\n    UNION ALL\n    SELECT 'RELEASE_NATIVE_PANEL', n.release, '', 'RELEASE_NATIVE', 'DISEASE_PANEL', 'ELIGIBLE_RELEASE_NATIVE', COUNT(*), den.candidate_panels, COUNT(*)::DOUBLE / NULLIF(den.candidate_panels, 0), COUNT(DISTINCT n.canonical_disease_id), NULL, NULL, NULL\n    FROM {native} n JOIN native_denominators den USING (release)\n    WHERE n.release_native_panel_eligible = 'true'\n    GROUP BY n.release, den.candidate_panels\n    UNION ALL\n    SELECT 'EVIDENCE_DOMAIN_PANEL', '', e.release_pair, e.endpoint_level || ':' || e.endpoint_id, 'DISEASE_PANEL', 'ELIGIBLE_PRIMARY_30', COUNT(*), den.candidate_panels, COUNT(*)::DOUBLE / NULLIF(den.candidate_panels, 0), COUNT(DISTINCT e.canonical_disease_id), NULL, NULL, NULL\n    FROM {evidence} e JOIN evidence_denominators den USING (endpoint_level, endpoint_id, release_pair)\n    WHERE e.panel_eligible = 'true'\n    GROUP BY e.release_pair, e.endpoint_level, e.endpoint_id, den.candidate_panels\n)\nSELECT scope, release, release_pair, support_name, unit, stage, count, denominator, fraction, disease_count, target_count, therapeutic_area_count, score_mass\nFROM flow\nORDER BY scope, release, release_pair, support_name, unit, stage\n"

def run(root, configuration, threads = 8, memory_limit = "8GB"):
    root = Path(root).resolve()
    args = SimpleNamespace(threads = threads, memory_limit = memory_limit)
    inputs = discover_inputs(root, configuration)
    outputs = {key: root / value for ((key, value)) in (configuration["outputs"].items())}
    cleanup_partial_outputs(configuration, root = root)
    connection = duckdb.connect()
    try:
        temp_directory = configure_resources(connection, configuration, args, root = root)
        cleanup_stage_outputs(temp_directory, root = root)
        partition_count = int(
            configuration["determinism"]["large_output_hash_partitions"]
        )
        row_group_size = int(configuration["determinism"]["parquet_row_group_size"])
        connection.execute(
            f"CREATE TEMP VIEW support_construction_disease_crosswalk AS SELECT * FROM {csv_relation(outputs['disease_crosswalk'])}"
        )
        connection.execute(
            f"CREATE TEMP VIEW support_construction_target_crosswalk AS SELECT * FROM {csv_relation(outputs['target_crosswalk'])}"
        )
        for (release) in (configuration["releases"]):
            contract = configuration["dataset_fields"]["association_overall_direct"][
                release
            ]
            invalid = connection.execute(
                required_assertion_sql(
                    inputs[release, "association_overall_direct"], contract
                )
            ).fetchone()[0]
            if (int(invalid or 0) != 0):
                raise RuntimeError(
                    f"Overall aggregation assertion failed for {release}: {invalid}"
                )
        canonical_queries = []
        for (release) in (configuration["releases"]):
            for (build_partition) in (range(partition_count)):
                label = f"{release}_p{build_partition:02d}"
                canonical_queries.append(
                    (
                        label,
                        canonical_association_query(
                            configuration,
                            inputs,
                            release,
                            build_partition,
                            partition_count,
                        ),
                    )
                )
        build_partitioned_parquet(
            connection,
            canonical_queries,
            outputs["canonical_associations"],
            temp_directory,
            "support_construction_canonical",
            row_group_size,
            root = root,
        )
        transition_queries = []
        for (old_release, new_release) in (configuration["release_pairs"]):
            pair = f"{old_release}_to_{new_release}"
            for (build_partition) in (range(partition_count)):
                label = f"{pair}_p{build_partition:02d}"
                transition_queries.append(
                    (
                        label,
                        state_transition_query(
                            outputs["canonical_associations"],
                            old_release,
                            new_release,
                            build_partition,
                        ),
                    )
                )
        build_partitioned_parquet(
            connection,
            transition_queries,
            outputs["association_state_transitions"],
            temp_directory,
            "support_construction_transitions",
            row_group_size,
            root = root,
        )
        triple_queries = []
        for (build_partition) in (range(partition_count)):
            triple_queries.append(
                (
                    f"p{build_partition:02d}",
                    three_release_query(
                        outputs["canonical_associations"], build_partition
                    ),
                )
            )
        build_partitioned_parquet(
            connection,
            triple_queries,
            outputs["three_release_score_support"],
            temp_directory,
            "support_construction_three_release",
            row_group_size,
            root = root,
        )
        fixed_queries = []
        for (support_name) in (("PAIR_SPECIFIC_PERSISTENT", "THREE_RELEASE_PERSISTENT")):
            for (old_release, new_release) in (configuration["release_pairs"]):
                pair = f"{old_release}_to_{new_release}"
                for (build_partition) in (range(partition_count)):
                    label = f"{support_name}_{pair}_p{build_partition:02d}"
                    fixed_queries.append(
                        (
                            label,
                            fixed_members_query(
                                outputs["canonical_associations"],
                                outputs["three_release_score_support"],
                                support_name,
                                old_release,
                                new_release,
                                build_partition,
                            ),
                        )
                    )
        build_partitioned_parquet(
            connection,
            fixed_queries,
            outputs["fixed_ranking_panel_members"],
            temp_directory,
            "support_construction_fixed_members",
            row_group_size,
            root = root,
        )
        copy_query(
            connection,
            fixed_registry_query(
                configuration,
                outputs["fixed_ranking_panel_members"],
                outputs["disease_crosswalk"],
            ),
            outputs["fixed_ranking_panel_registry"],
            "CSV",
        )
        native_queries = []
        for (release) in (configuration["releases"]):
            for (build_partition) in (range(partition_count)):
                label = f"{release}_p{build_partition:02d}"
                native_queries.append(
                    (
                        label,
                        native_members_query(
                            outputs["canonical_associations"], release, build_partition
                        ),
                    )
                )
        build_partitioned_parquet(
            connection,
            native_queries,
            outputs["release_native_panel_members"],
            temp_directory,
            "support_construction_native_members",
            row_group_size,
            root = root,
        )
        copy_query(
            connection,
            native_registry_query(
                configuration,
                outputs["release_native_panel_members"],
                outputs["disease_crosswalk"],
            ),
            outputs["release_native_panel_registry"],
            "CSV",
        )
        evidence_queries = []
        evidence_selections = (
            (
                "DATATYPE",
                "association_by_datatype_direct",
                configuration["evidence_domain"]["datatype_id"],
            ),
            (
                "DATASOURCE",
                "association_by_datasource_direct",
                configuration["evidence_domain"]["datasource_id"],
            ),
        )
        for (endpoint_level, dataset, endpoint_id) in (evidence_selections):
            for (release) in (configuration["releases"]):
                for (build_partition) in (range(partition_count)):
                    label = f"{endpoint_level}_{endpoint_id}_{release}_p{build_partition:02d}"
                    query = evidence_domain_query(
                        configuration,
                        inputs,
                        endpoint_level,
                        dataset,
                        endpoint_id,
                        release,
                        build_partition,
                        partition_count,
                    )
                    evidence_queries.append((label, query))
        build_partitioned_parquet(
            connection,
            evidence_queries,
            outputs["evidence_domain_support"],
            temp_directory,
            "support_construction_evidence",
            row_group_size,
            root = root,
        )
        copy_query(
            connection,
            evidence_registry_query(configuration, outputs["evidence_domain_support"]),
            outputs["evidence_domain_panel_registry"],
            "CSV",
        )
        entrant_queries = []
        for (old_release, new_release) in (configuration["release_pairs"]):
            pair = f"{old_release}_to_{new_release}"
            for (build_partition) in (range(partition_count)):
                label = f"{pair}_p{build_partition:02d}"
                entrant_queries.append(
                    (
                        label,
                        entrant_exit_query(
                            outputs["association_state_transitions"],
                            pair,
                            build_partition,
                        ),
                    )
                )
        build_partitioned_parquet(
            connection,
            entrant_queries,
            outputs["entrant_exit_table"],
            temp_directory,
            "support_construction_entrant_exit",
            row_group_size,
            root = root,
        )
        exclusion_queries = exclusion_registry_queries(
            configuration,
            outputs["canonical_associations"],
            outputs["fixed_ranking_panel_registry"],
            partition_count,
            root = root,
        )
        build_partitioned_parquet(
            connection,
            exclusion_queries,
            outputs["exclusion_registry"],
            temp_directory,
            "support_construction_exclusions",
            row_group_size,
            root = root,
        )
        copy_query(
            connection,
            support_flow_query(
                configuration,
                outputs["canonical_associations"],
                outputs["association_state_transitions"],
                outputs["three_release_score_support"],
                outputs["fixed_ranking_panel_registry"],
                outputs["release_native_panel_registry"],
                outputs["evidence_domain_panel_registry"],
                outputs["disease_crosswalk"],
                outputs["target_crosswalk"],
                outputs["therapeutic_area_mapping"],
            ),
            outputs["support_flow"],
            "CSV",
        )
    except Exception:
        cleanup_partial_outputs(configuration, root = root)
        if ("temp_directory" in locals()):
            cleanup_stage_outputs(temp_directory, root = root)
        raise
    finally:
        connection.close()
