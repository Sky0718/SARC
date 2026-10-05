from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from .comparators import apply_minimum_roster, compare_panels
from .contracts import PanelSnapshot
from .io import check, quoted, read_json, write_json
from .string_transfer import (
    compact_panel,
    grouped_panels,
    leading_sets,
    panel_schema,
    summarise,
)
from .string_transfer_inputs import copy_query, info_rows, single_row

def canonicalise_links(source, destination):
    destination = Path(destination)
    destination.parent.mkdir(parents = True, exist_ok = True)
    with duckdb.connect(config = {"threads": "2", "memory_limit": "3GB"}) as connection:
        connection.execute(
            f"CREATE VIEW raw AS SELECT * FROM read_csv({quoted(source)}, delim = ' ', header = true, compression = 'gzip', parallel = false)"
        )
        invalid = connection.execute(
            "SELECT count(*) FROM raw WHERE protein1 IS NULL OR protein2 IS NULL "
            "OR NOT starts_with(protein1, '9606.') OR NOT starts_with(protein2, '9606.') "
            "OR combined_score IS NULL OR combined_score < 0 OR combined_score > 1000 "
            "OR combined_score != floor(combined_score)"
        ).fetchone()[0]
        duplicate = connection.execute(
            "SELECT count(*) FROM (SELECT protein1, protein2 FROM raw GROUP BY ALL HAVING count(*) > 1)"
        ).fetchone()[0]
        check(
            (invalid == 0 and duplicate == 0),
            "Invalid score or duplicate directed edge",
        )
        connection.execute(
            "CREATE TEMP TABLE edges AS SELECT least(protein1, protein2) AS protein_a, greatest(protein1, protein2) AS protein_b, min(combined_score)::SMALLINT AS combined_score, max(combined_score) AS maximum_score, count(*) AS source_rows FROM raw WHERE protein1 != protein2 GROUP BY protein_a, protein_b"
        )
        invalid = connection.execute(
            "SELECT count(*) FROM edges WHERE combined_score != maximum_score OR source_rows > 2"
        ).fetchone()[0]
        check((invalid == 0), "Reciprocal scores or multiplicities disagree")
        copy_query(
            connection,
            "SELECT protein_a, protein_b, combined_score FROM edges ORDER BY protein_a, protein_b",
            destination,
        )

def prepare_string(config_path):
    config = read_json(config_path)
    base = Path(config_path).resolve().parent
    source = {
        key: (base / value).resolve() for ((key, value)) in (config["inputs"].items())
    }
    output = (base / config["output"]).resolve()
    output.mkdir(parents = True, exist_ok = True)
    with duckdb.connect(config = {"threads": "2", "memory_limit": "3GB"}) as connection:
        for (state) in (("baseline", "followup")):
            connection.execute(
                f"CREATE VIEW {state}_edges AS SELECT * FROM read_parquet({quoted(source[state + '_edges'])})"
            )
            connection.register(
                state + "_info",
                pa.Table.from_pylist(info_rows(source[state + "_info"])),
            )
            invalid = connection.execute(
                f"SELECT count(*) FROM {state}_edges e "
                f"LEFT JOIN {state}_info a ON e.protein_a = a.anchor_id "
                f"LEFT JOIN {state}_info b ON e.protein_b = b.anchor_id "
                "WHERE a.anchor_id IS NULL OR b.anchor_id IS NULL"
            ).fetchone()[0]
            check(
                (invalid == 0), "An edge endpoint is absent from its release metadata"
            )
        below_baseline_threshold = connection.execute(
            "SELECT count(*) FROM baseline_edges WHERE combined_score < 150"
        ).fetchone()[0]
        check(
            (below_baseline_threshold == 0),
            "The baseline archive must retain its 150 score floor",
        )
        connection.execute(
            "CREATE TEMP TABLE baseline_degrees AS SELECT anchor_id, count(*)::INTEGER AS baseline_degree_150 FROM (SELECT protein_a AS anchor_id FROM baseline_edges UNION ALL SELECT protein_b AS anchor_id FROM baseline_edges) GROUP BY anchor_id"
        )
        query = "SELECT i.anchor_id, i.label, coalesce(d.baseline_degree_150, 0)::INTEGER AS baseline_degree_150, coalesce(d.baseline_degree_150, 0) >= 20 AS cohort_20, coalesce(d.baseline_degree_150, 0) >= 30 AS cohort_30, coalesce(d.baseline_degree_150, 0) >= 50 AS cohort_50, CASE WHEN coalesce(d.baseline_degree_150, 0) < 20 THEN 'below_20' WHEN d.baseline_degree_150 < 30 THEN '20_29' WHEN d.baseline_degree_150 < 100 THEN '30_99' WHEN d.baseline_degree_150 < 500 THEN '100_499' ELSE '500_plus' END AS baseline_degree_stratum FROM baseline_info i LEFT JOIN baseline_degrees d USING(anchor_id) ORDER BY i.anchor_id"
        copy_query(connection, query, output / "baseline_anchor_registry.parquet")
        connection.execute(
            f"CREATE VIEW baseline_registry AS SELECT * FROM read_parquet({quoted(output / 'baseline_anchor_registry.parquet')})"
        )
        total_degree = connection.execute(
            "SELECT sum(baseline_degree_150) FROM baseline_registry"
        ).fetchone()[0]
        edges = connection.execute("SELECT count(*) FROM baseline_edges").fetchone()[0]
        check(
            (total_degree == (2 * edges)),
            "Edge identifiers are absent from baseline metadata",
        )
        copy_query(
            connection,
            "SELECT b.anchor_id, f.anchor_id IS NOT NULL AS present_in_followup_metadata, CASE WHEN f.anchor_id IS NOT NULL THEN 'EXACT_ID' ELSE 'UNRESOLVED' END AS anchor_identity_state FROM baseline_registry b LEFT JOIN followup_info f USING(anchor_id) ORDER BY b.anchor_id",
            output / "anchor_identity_coverage.parquet",
        )
        connection.execute(
            "CREATE TEMP TABLE edge_transitions AS SELECT coalesce(b.protein_a, f.protein_a) AS protein_a, coalesce(b.protein_b, f.protein_b) AS protein_b, b.combined_score AS baseline_score, f.combined_score AS followup_score, b.protein_a IS NOT NULL AS baseline_observed, f.protein_a IS NOT NULL AS followup_observed FROM baseline_edges b FULL OUTER JOIN followup_edges f USING(protein_a, protein_b)"
        )
        copy_query(
            connection,
            "SELECT * FROM edge_transitions",
            output / "edge_transitions.parquet",
        )
        write_json(
            output / "global_edge_summaries.json", global_edge_summaries(connection)
        )
        query = "WITH directed AS (SELECT protein_a AS anchor_id, protein_b AS neighbour_id, baseline_score, followup_score FROM edge_transitions UNION ALL SELECT protein_b AS anchor_id, protein_a AS neighbour_id, baseline_score, followup_score FROM edge_transitions) SELECT d.* FROM directed d INNER JOIN baseline_registry b USING(anchor_id) WHERE b.cohort_20 ORDER BY d.anchor_id, d.neighbour_id"
        copy_query(connection, query, output / "anchor_neighbour_transitions.parquet")

def global_edge_summaries(connection):
    summaries = []
    for (threshold) in ((150, 400, 700)):
        row = single_row(
            connection,
            f"SELECT count(*) FILTER(WHERE baseline_score >= {threshold}) AS baseline_eligible_edges, count(*) FILTER(WHERE followup_score >= {threshold}) AS followup_eligible_edges, count(*) FILTER(WHERE baseline_score >= {threshold} AND followup_score >= {threshold}) AS common_eligible_edges, count(*) FILTER(WHERE baseline_score >= {threshold} AND followup_score < {threshold}) AS baseline_eligible_followup_observed_below_threshold, count(*) FILTER(WHERE baseline_score >= {threshold} AND NOT followup_observed) AS baseline_eligible_followup_archive_absent, count(*) FILTER(WHERE followup_score >= {threshold} AND baseline_score < {threshold}) AS followup_eligible_baseline_observed_below_threshold, count(*) FILTER(WHERE followup_score >= {threshold} AND NOT baseline_observed) AS followup_eligible_baseline_archive_absent, quantile_cont(abs(followup_score-baseline_score), [0.5,0.9,0.95,0.99]) FILTER(WHERE baseline_score >= {threshold} AND followup_score >= {threshold}) AS absolute_change_quantiles, avg(followup_score-baseline_score) FILTER(WHERE baseline_score >= {threshold} AND followup_score >= {threshold}) AS mean_signed_change, max(abs(followup_score-baseline_score)) FILTER(WHERE baseline_score >= {threshold} AND followup_score >= {threshold}) AS maximum_absolute_change FROM edge_transitions",
        )
        row.update(
            {
                "resource": "STRING",
                "baseline_release": "11.5",
                "followup_release": "12.0",
                "threshold": threshold,
                "unit": "canonical_undirected_edge",
                "quantiles": [0.5, 0.9, 0.95, 0.99],
                "score_scale": "native_combined_score_integer_0_1000",
                "confidence_interval": None,
            }
        )
        check(
            row["baseline_eligible_edges"]
            == row["common_eligible_edges"]
            + row["baseline_eligible_followup_observed_below_threshold"]
            + row["baseline_eligible_followup_archive_absent"],
            "Threshold baseline decomposition differs",
        )
        check(
            row["followup_eligible_edges"]
            == row["common_eligible_edges"]
            + row["followup_eligible_baseline_observed_below_threshold"]
            + row["followup_eligible_baseline_archive_absent"],
            "Threshold followup decomposition differs",
        )
        summaries.append(row)
    return summaries

def strategy_records(report, threshold):
    for (strategy, details) in (report["strategies"].items()):
        for (task, answer) in (details["answers"].items()):
            reference = report["strategies"]["support_explicit"]["answers"][task]
            yield {
                "anchor_id": report["panel_id"],
                "threshold": threshold,
                "minimum_roster": report["minimum_roster"],
                "strategy": strategy,
                "task": task,
                "eligible": reference["state"] == "ANSWERED",
                "state": answer["state"],
                "boolean_answer": answer["answer"]
                if (isinstance(answer["answer"], bool))
                else None,
                "support_denominator": answer["support_denominator"],
                "reason": answer["reason"],
                "agrees_with_full_report": answer["answer"] == reference["answer"]
                if (answer["state"] == "ANSWERED")
                else None,
            }

def compare_string(config_path):
    config = read_json(config_path)
    base = Path(config_path).resolve().parent
    source = (base / config["source"]).resolve()
    output = (base / config["output"]).resolve()
    output.mkdir(parents = True, exist_ok = True)
    settings = {
        "thresholds": [150, 400, 700],
        "minimum_rosters": [20, 30, 50],
        "threshold_contrasts": [[150, 400], [150, 700], [400, 700]],
        "quantiles": [0.5, 0.9, 0.95, 0.99],
    }
    registry = {
        row["anchor_id"]: row
        for (row) in (
            pq.read_table(source / "baseline_anchor_registry.parquet").to_pylist()
        )
        if (row["cohort_20"])
    }
    identities = {
        row["anchor_id"]: row["present_in_followup_metadata"]
        for (row) in (
            pq.read_table(source / "anchor_identity_coverage.parquet").to_pylist()
        )
    }
    writers = {}
    pending = {"panels": [], "sets": [], "strategies": []}
    seen = set()
    for (anchor, records) in (grouped_panels(
        source / "anchor_neighbour_transitions.parquet"
    )):
        check(
            (anchor in registry and anchor not in seen),
            "Unexpected or duplicate anchor",
        )
        seen.add(anchor)
        for (threshold) in (settings["thresholds"]):
            first = tuple(
                (row["neighbour_id"], row["baseline_score"])
                for (row) in (records)
                if (
                    (row["baseline_score"] is not None)
                    and (row["baseline_score"] >= threshold)
                )
            )
            second = tuple(
                (row["neighbour_id"], row["followup_score"])
                for (row) in (records)
                if (
                    (row["followup_score"] is not None)
                    and (row["followup_score"] >= threshold)
                )
            )
            report = compare_panels(
                PanelSnapshot("STRING", "11.5", anchor, first),
                PanelSnapshot("STRING", "12.0", anchor, second, identities[anchor]),
                minimum_roster = 20,
            )
            pending["sets"].extend(leading_sets(report, threshold))
            for (minimum) in (settings["minimum_rosters"]):
                if (registry[anchor]["baseline_degree_150"] >= minimum):
                    view = (
                        report
                        if (minimum == 20)
                        else apply_minimum_roster(report, minimum)
                    )
                    pending["panels"].append(
                        compact_panel(view, threshold, registry[anchor])
                    )
                    pending["strategies"].extend(strategy_records(view, threshold))
        if ((len(seen) % 250) == 0):
            flush_records(output, pending, writers)
    flush_records(output, pending, writers)
    for (writer) in (writers.values()):
        writer.close()
    check((seen == set(registry)), "The full baseline cohort is required")
    summarise(output, settings)
    from .string_minimum_contrasts import contrast_records

    with duckdb.connect() as connection:
        connection.execute(
            f"CREATE VIEW panels AS SELECT * FROM read_parquet({quoted(output / 'panel_endpoints.parquet')})"
        )
        write_json(
            output / "minimum_roster_contrasts.json",
            contrast_records(
                connection, settings["minimum_rosters"], settings["thresholds"]
            ),
        )

def flush_records(output, pending, writers):
    names = {
        "panels": "panel_endpoints.parquet",
        "sets": "leading_sets.parquet",
        "strategies": "strategy_reports.parquet",
    }
    schemas = {
        "panels": panel_schema(),
        "sets": pa.schema(
            [
                ("anchor_id", pa.string()),
                ("threshold", pa.int32()),
                ("minimum_roster", pa.int32()),
                ("view", pa.string()),
                ("endpoint", pa.string()),
                ("baseline_set", pa.list_(pa.string())),
                ("followup_set", pa.list_(pa.string())),
                ("reason", pa.string()),
            ]
        ),
        "strategies": pa.schema(
            [
                ("anchor_id", pa.string()),
                ("threshold", pa.int32()),
                ("minimum_roster", pa.int32()),
                ("strategy", pa.string()),
                ("task", pa.string()),
                ("eligible", pa.bool_()),
                ("state", pa.string()),
                ("boolean_answer", pa.bool_()),
                ("support_denominator", pa.int32()),
                ("reason", pa.string()),
                ("agrees_with_full_report", pa.bool_()),
            ]
        ),
    }
    for (name, rows) in (pending.items()):
        if (rows):
            if (name not in writers):
                writers[name] = pq.ParquetWriter(
                    output / names[name], schemas[name], compression = "zstd"
                )
            writers[name].write_table(pa.Table.from_pylist(rows, schema = schemas[name]))
            rows.clear()
