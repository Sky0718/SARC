import itertools
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from .comparators import apply_minimum_roster, compare_panels
from .contracts import PanelSnapshot
from .io import check, read_json, write_json, write_parquet
from .temporal_analysis import (
    MINIMUMS,
    PRIMARY_MINIMUM,
    TOP_KS,
    flatten_report,
    minimum_comparisons,
    register_row,
    rows_by_panel,
    summarise,
    warning_summaries,
)
from .temporal_identity import (
    DISEASE_FIELDS,
    TARGET_FIELDS,
    disease_alignment,
    target_alignment,
)
from .temporal_prepare import literal, parquet_source, read_entities

def prepare_temporal(config_path):
    config = read_json(config_path)
    base = Path(config_path).resolve().parent
    groups = {
        release: {
            family: [(base / name).resolve() for (name) in (names)]
            for ((family, names)) in (families.items())
        }
        for ((release, families)) in (config["sources"].items())
    }
    output = (base / config["output"]).resolve()
    output.mkdir(parents = True, exist_ok = True)
    diseases = {
        release: read_entities(families["disease"], DISEASE_FIELDS)
        for ((release, families)) in (groups.items())
    }
    targets = {
        release: read_entities(families["target"], TARGET_FIELDS)
        for ((release, families)) in (groups.items())
    }
    mappings, evidence = disease_alignment(diseases["26.06"], diseases["26.09"])
    target_rows = target_alignment(targets["26.06"], targets["26.09"])
    write_parquet(output / "disease_identity.parquet", mappings)
    write_parquet(output / "target_identity.parquet", target_rows)
    write_json(output / "mapping_evidence.json", evidence)
    registry = pq.read_table((base / config["baseline_registry"]).resolve()).to_pylist()
    by_id = {row["baseline_disease_id"]: row for (row) in (mappings)}
    cohort = [
        {**row, **by_id[row["disease_id"]]}
        for (row) in (registry)
        if (row["eligible_disease"] and row["eligible_min_20"])
    ]
    check(
        (sum(row["eligible_min_30"] for (row) in (cohort)) == 11897),
        "The complete baseline cohort is required",
    )
    write_parquet(output / "cohort_followup_identity.parquet", cohort)
    with duckdb.connect(config = {"threads": "2", "memory_limit": "3GB"}) as connection:
        connection.register("cohort", pa.Table.from_pylist(cohort))
        for (release, label) in ((("26.06", "baseline"), ("26.09", "followup"))):
            connection.register(
                label + "_targets", pa.table({"target_id": sorted(targets[release])})
            )
            join = (
                "a.diseaseId = c.disease_id"
                if (label == "baseline")
                else "a.diseaseId = c.followup_disease_id"
            )
            source = parquet_source(groups[release]["association_overall_direct"])
            query = f"SELECT c.disease_id, a.targetId AS target_id, a.associationScore AS score, a.evidenceCount AS evidence_count, true AS record_present, t.target_id IS NOT NULL AS target_present, CASE WHEN a.associationScore IS NULL THEN 'PRESENT_NULL' WHEN a.associationScore = 0 THEN 'PRESENT_ZERO' ELSE 'PRESENT_POSITIVE' END AS record_state, coalesce(a.associationScore > 0 AND isfinite(a.associationScore) AND a.evidenceCount >= 1 AND t.target_id IS NOT NULL, false) AS eligible FROM {source} a JOIN cohort c ON {join} LEFT JOIN {label}_targets t ON a.targetId = t.target_id WHERE a.aggregationType = 'overall' AND a.aggregationValue = 'None'"
            connection.execute(f"CREATE VIEW {label}_rows AS {query}")
        query = "SELECT coalesce(b.disease_id, n.disease_id) AS disease_id, coalesce(b.target_id, n.target_id) AS target_id, coalesce(b.record_present, false) AS baseline_record_present, coalesce(n.record_present, false) AS followup_record_present, b.score AS baseline_score, n.score AS followup_score, b.evidence_count AS baseline_evidence_count, n.evidence_count AS followup_evidence_count, coalesce(b.eligible, false) AS baseline_eligible, coalesce(n.eligible, false) AS followup_eligible, coalesce(b.record_state, 'ABSENT_OR_CENSORED') AS baseline_record_state, CASE WHEN c.followup_disease_id IS NULL THEN 'MAPPING_UNRESOLVED' ELSE coalesce(n.record_state, 'ABSENT_OR_CENSORED') END AS followup_record_state, c.mapping_state AS disease_mapping_state FROM baseline_rows b FULL OUTER JOIN followup_rows n USING(disease_id, target_id) JOIN cohort c ON c.disease_id = coalesce(b.disease_id, n.disease_id)"
        connection.execute(
            f"COPY ({query} ORDER BY disease_id, target_id) TO {literal(output / 'association_support.parquet')} (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 65536)"
        )
        connection.execute(
            f"CREATE VIEW pairs AS SELECT * FROM read_parquet({literal(output / 'association_support.parquet')})"
        )
        duplicate = connection.execute(
            "SELECT count(*) FROM (SELECT disease_id, target_id FROM pairs GROUP BY ALL HAVING count(*) > 1)"
        ).fetchone()[0]
        check((duplicate == 0), "Duplicate association key")
        sizes = dict(
            connection.execute(
                "SELECT disease_id, count(*) FILTER (WHERE baseline_eligible) FROM pairs GROUP BY disease_id"
            ).fetchall()
        )
        check(
            (set(sizes) == {row["disease_id"] for (row) in (cohort)}),
            "Missing baseline disease",
        )
        check(
            all(
                sizes[row["disease_id"]] == row["native_roster_size"]
                for (row) in (cohort)
            ),
            "Baseline roster size differs",
        )

def compare_temporal(source, output):
    source, output = Path(source), Path(output)
    output.mkdir(parents = True, exist_ok = True)
    cohort = {
        row["disease_id"]: row
        for (row) in (
            pq.read_table(source / "cohort_followup_identity.parquet").to_pylist()
        )
    }
    rows, sets, rbo, strategies = [], [], [], []
    seen = set()
    for (panel_id, records) in (rows_by_panel(source / "association_support.parquet")):
        check(
            (panel_id in cohort and panel_id not in seen),
            "Unexpected or duplicate disease panel",
        )
        seen.add(panel_id)
        baseline = cohort[panel_id]
        first = tuple(
            (row["target_id"], row["baseline_score"])
            for (row) in (records)
            if (row["baseline_eligible"])
        )
        second = tuple(
            (row["target_id"], row["followup_score"])
            for (row) in (records)
            if (row["followup_eligible"])
        )
        left = PanelSnapshot("OpenTargets", "26.06", panel_id, first)
        right = PanelSnapshot(
            "OpenTargets",
            "26.09",
            panel_id,
            second,
            baseline["followup_disease_id"] is not None,
        )
        base_report = compare_panels(
            left,
            right,
            minimum_roster = 20,
            top_ks = TOP_KS,
            include_rbo = len(first) >= PRIMARY_MINIMUM,
        )
        for (minimum) in (MINIMUMS):
            if (len(first) < minimum):
                continue
            report = (
                base_report
                if (minimum == 20)
                else apply_minimum_roster(base_report, minimum)
            )
            row, selected = flatten_report(report, baseline, minimum)
            rows.append(row)
            sets.extend(selected)
            if (minimum == PRIMARY_MINIMUM):
                for (view, variant, persistence) in (itertools.product(
                    ("native", "fixed"), ("w", "b"), (0.9, 0.8, 0.95)
                )):
                    estimate = report[view]["rbo"]
                    values = (
                        estimate["value"][variant][str(persistence)]
                        if (estimate["estimable"])
                        else {key: None for (key) in (("min", "max", "res", "ext"))}
                    )
                    rbo.append(
                        {
                            "disease_id": panel_id,
                            "support": view,
                            "variant": variant,
                            "persistence": persistence,
                            "reason": estimate["reason"],
                            **values,
                        }
                    )
                for (strategy, details) in (report["strategies"].items()):
                    for (task, answer) in (details["answers"].items()):
                        strategies.append(
                            {
                                "disease_id": panel_id,
                                "strategy": strategy,
                                "task": task,
                                "state": answer["state"],
                                "boolean_answer": answer["answer"]
                                if (isinstance(answer["answer"], bool))
                                else None,
                                "support_denominator": answer["support_denominator"],
                                "reason": answer["reason"],
                            }
                        )
    check((seen == set(cohort)), "Incomplete baseline cohort")
    check(
        (sum(row["minimum_roster"] == 30 for (row) in (rows)) == 11897),
        "Incomplete primary endpoint population",
    )
    for (name, records) in ((
        ("panel_endpoints", rows),
        ("leading_sets", sets),
        ("rank_biased_overlap", rbo),
        ("strategy_answers", strategies),
    )):
        write_parquet(output / (name + ".parquet"), records)
    registry = summarise(rows)
    for (view, variant, persistence, quantity) in (itertools.product(
        ("native", "fixed"), ("w", "b"), (0.9, 0.8, 0.95), ("min", "max", "res", "ext")
    )):
        selected = [
            row
            for (row) in (rbo)
            if (
                (row["support"] == view)
                and (row["variant"] == variant)
                and (row["persistence"] == persistence)
            )
        ]
        registry.append(
            register_row(
                f"{view}_rbo_{variant}_p{persistence}_{quantity}_mean",
                [row[quantity] for (row) in (selected)],
                11897,
                kind = "scalar",
                reasons = [row["reason"] for (row) in (selected)],
            )
        )
    warnings, confusion = warning_summaries(rows)
    registry.extend(warnings)
    write_json(output / "endpoint_register.json", registry)
    write_json(output / "warning_rule_assessment.json", confusion)
    write_json(output / "minimum_roster_comparison.json", minimum_comparisons(rows))
