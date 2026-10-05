from __future__ import annotations
from pathlib import Path
from typing import Any
from typing import Mapping
from typing import Sequence
from .ranking_types import DuckDbDataset, RankingSensitivityBuildError

def write_parquet_dataset(
    connection: Any,
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    fields: Sequence[str],
    types: Sequence[str],
    order: Sequence[str],
    table_suffix: str,
    batch_size: int = 5000,
) -> None:
    from .ranking_1 import sql_identifier, sql_literal

    if (len(fields) != len(types)):
        raise RankingSensitivityBuildError(
            "Parquet field and type declarations differ in length"
        )
    if (batch_size < 1):
        raise ValueError("Batch size must be positive")
    table = f"ranking_sensitivity_export_{table_suffix}"
    connection.execute(f"DROP TABLE IF EXISTS {sql_identifier(table)}")
    declaration = ", ".join(
        (
            f"{sql_identifier(field)} {data_type}"
            for ((field, data_type)) in (zip(fields, types))
        )
    )
    connection.execute(f"CREATE TABLE {sql_identifier(table)} ({declaration})")
    if (rows):
        parameters = ",".join(("?" for (_) in (fields)))
        statement = f"INSERT INTO {sql_identifier(table)} VALUES ({parameters})"
        for (start) in (range(0, len(rows), batch_size)):
            batch = rows[start : start + batch_size]
            connection.executemany(
                statement,
                [[row.get(field) for (field) in (fields)] for (row) in (batch)],
            )
    path.parent.mkdir(parents = True, exist_ok = True)
    order_sql = ", ".join((sql_identifier(field) for (field) in (order)))
    connection.execute(
        f"COPY (SELECT * FROM {sql_identifier(table)} ORDER BY {order_sql}) TO {sql_literal(path.as_posix())} (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000)"
    )
    connection.execute(f"DROP TABLE {sql_identifier(table)}")

def write_parquet_table_dataset(
    connection: Any, path: Path, dataset: DuckDbDataset, order: Sequence[str]
) -> None:
    from .ranking_1 import flush_duckdb_dataset, sql_identifier, sql_literal

    flush_duckdb_dataset(connection, dataset)
    path.parent.mkdir(parents = True, exist_ok = True)
    order_sql = ", ".join((sql_identifier(field) for (field) in (order)))
    connection.execute(
        f"COPY (SELECT * FROM {sql_identifier(dataset.table)} ORDER BY {order_sql}) TO {sql_literal(path.as_posix())} (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000)"
    )

def export_stage(
    connection: Any,
    stage: Path,
    datasets: Mapping[str, Any],
    configuration: Mapping[str, Any],
) -> dict[str, Path]:
    from .ranking_2 import output_specifications, write_csv_dataset

    specifications = output_specifications()
    paths: dict[str, Path] = {}
    for (name, specification) in (specifications.items()):
        connection.execute(
            f"SET threads = {int(configuration['resources']['export_threads'])}"
        )
        destination = stage / Path(configuration["outputs"][name]).name
        rows = datasets[name]
        if (isinstance(rows, DuckDbDataset)):
            if (not specification["parquet"]):
                raise RankingSensitivityBuildError(
                    "A DuckDB-backed dataset must use Parquet output"
                )
            write_parquet_table_dataset(
                connection, destination, rows, specification["order"]
            )
        elif (specification["parquet"]):
            write_parquet_dataset(
                connection,
                destination,
                rows,
                specification["fields"],
                specification["types"],
                specification["order"],
                f"{name}_{stage.name}",
            )
        else:
            write_csv_dataset(
                destination, rows, specification["fields"], specification["order"]
            )
        paths[name] = destination
    return paths

def build_ranking_sensitivity_datasets(
    connection: Any, configuration: Mapping[str, Any]
) -> dict[str, Any]:
    from .ranking_2 import (
        build_disease_registry,
        build_fixed_outputs,
        build_native_outputs,
        build_sensitivity_signal_assessment,
        build_therapeutic_area_summary,
        build_turnover_decomposition,
        output_specifications,
        stable_area_memberships,
    )
    from .ranking_1 import create_duckdb_dataset, flush_duckdb_datasets

    top_member_specification = output_specifications()["top_set_members"]
    top_member_dataset = create_duckdb_dataset(
        connection,
        "ranking_sensitivity_top_set_member_data",
        top_member_specification["fields"],
        top_member_specification["types"],
        int(configuration["resources"]["bulk_insert_batch_rows"]),
    )
    (
        fixed_metrics,
        displacement_dataset,
        fixed_top_rows,
        primary_fixed,
        fixed_diagnostics,
    ) = build_fixed_outputs(connection, configuration, top_member_dataset)
    (
        native_metrics,
        native_top_rows,
        primary_native,
        turnover_source,
        native_diagnostics,
    ) = build_native_outputs(
        connection, configuration, primary_fixed, top_member_dataset
    )
    disease_rows = build_disease_registry(primary_fixed, primary_native, configuration)
    memberships, area_classifications, area_source_diagnostics = (
        stable_area_memberships(connection, configuration)
    )
    for (row) in (disease_rows):
        mapping_key = (str(row["old_release"]), str(row["canonical_disease_id"]))
        row["therapeutic_area_membership_status"] = area_classifications.get(
            mapping_key,
            configuration["therapeutic_areas"]["no_explicit_mapping_area_status"],
        )
    therapeutic_rows, therapeutic_diagnostics, therapeutic_coverage = (
        build_therapeutic_area_summary(
            disease_rows, memberships, area_classifications, configuration
        )
    )
    turnover_rows = build_turnover_decomposition(turnover_source, configuration)
    signal_rows, signal_decision = build_sensitivity_signal_assessment(
        disease_rows, therapeutic_rows, configuration
    )
    datasets: dict[str, Any] = {
        "disease_ranking_panels": disease_rows,
        "fixed_target_rank_displacement": displacement_dataset,
        "fixed_ranking_metrics": fixed_metrics,
        "release_native_top_set_metrics": native_metrics,
        "top_set_members": top_member_dataset,
        "top_ranked_target_sets": fixed_top_rows + native_top_rows,
        "turnover_reason_decomposition": turnover_rows,
        "therapeutic_area_ranking_summary": therapeutic_rows,
        "sensitivity_signal_assessment": signal_rows,
    }
    flush_duckdb_datasets(connection, datasets.values())
    diagnostics = (
        fixed_diagnostics
        | native_diagnostics
        | area_source_diagnostics
        | therapeutic_diagnostics
    )
    failures = {name: count for ((name, count)) in (diagnostics.items()) if (count)}
    if (failures):
        raise RankingSensitivityBuildError(
            f"Input or calculation contract failed: {failures}"
        )
    return datasets
