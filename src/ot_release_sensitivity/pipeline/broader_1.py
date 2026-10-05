from __future__ import annotations
from types import SimpleNamespace
import math
from pathlib import Path
import duckdb
import numpy as np
import pandas as pd
from scipy.stats import kendalltau
from .broader_types import (
    EPSILONS,
    EXCLUDED_CONTEXT_ROOTS,
    PANEL_THRESHOLDS,
    POPULATIONS,
    RELEASES,
    RELEASE_PAIRS,
    TAXONOMY,
)

def sql_path(path):
    return path.resolve().as_posix().replace("'", "''")

def score_column(release):
    return "score" if (release == "25.12") else "associationScore"

def association_path(project_root, release):
    return (
        project_root
        / "data"
        / "raw"
        / release
        / "association_overall_direct"
        / "*.parquet"
    )

def write_frame(frame, path, order):
    result = frame.copy()
    if (order and (not result.empty)):
        result = result.sort_values(order, kind = "mergesort").reset_index(drop = True)
    result.to_csv(
        path,
        index = False,
        encoding = "utf-8",
        lineterminator = "\n",
        na_rep = "",
        float_format = "%.17g",
    )
    return len(result)

def configure_connection(project_root, threads, memory_limit):
    temporary = project_root / "data" / "interim" / "broader_efo_duckdb_tmp"
    temporary.mkdir(parents = True, exist_ok = True)
    connection = duckdb.connect()
    connection.execute(f"SET threads = {int(threads)}")
    connection.execute(f"SET memory_limit = '{str(memory_limit).upper()}'")
    connection.execute(f"SET temp_directory = '{sql_path(temporary)}'")
    return connection

def split_identifiers(value):
    return {item for (item) in (str(value or "").split(";")) if (item)}

def load_broader_entities(project_root):
    path = project_root / "data" / "derived" / "disease_crosswalk.csv"
    frame = pd.read_csv(path, dtype = str, keep_default_na = False)
    mask = (
        (frame["source_namespace"] == "EFO")
        & (frame["exclusion_reason"] == "no_explicit_disease_authority_cross_reference")
        & (frame["source_entity_type"] == "EFO_WITHOUT_DISEASE_AUTHORITY")
    )
    result = frame.loc[
        mask,
        [
            "release",
            "source_disease_id",
            "source_label",
            "source_therapeutic_area_ids",
            "canonical_disease_id",
            "canonical_therapeutic_area_ids",
            "identifier_mapping_status",
            "ambiguity_status",
            "exclusion_reason",
        ],
    ].copy()
    result["excluded_ta_context_ids"] = result["source_therapeutic_area_ids"].map(
        lambda value: ";".join(
            sorted(split_identifiers(value) & EXCLUDED_CONTEXT_ROOTS)
        )
    )
    result["excluded_ta_context"] = result["excluded_ta_context_ids"].ne("")
    observed = result.groupby("release").size().to_dict()
    expected = {"25.12": 451, "26.03": 464, "26.06": 462}
    if (observed != expected):
        raise RuntimeError(
            f"Broader EFO entity counts differ from the frozen contract: {observed}"
        )
    duplicates = result.duplicated(["release", "source_disease_id"]).sum()
    if (duplicates):
        raise RuntimeError(
            "Broader EFO registry contains duplicate release-entity keys"
        )
    return result

def add_positive_support_counts(connection, project_root, entities):
    frames = []
    for (release) in (RELEASES):
        relation = entities.loc[
            entities["release"] == release, ["source_disease_id"]
        ].copy()
        connection.register("release_entities", relation)
        path = sql_path(association_path(project_root, release))
        score = score_column(release)
        query = f"\n            SELECT\n                a.diseaseId AS source_disease_id,\n                COUNT(*) AS positive_target_count\n            FROM read_parquet('{path}') a\n            JOIN release_entities e ON a.diseaseId = e.source_disease_id\n            WHERE isfinite(a.{score})\n              AND a.{score} > 0\n              AND a.evidenceCount IS NOT NULL\n              AND a.evidenceCount >= 1\n            GROUP BY a.diseaseId\n        "
        counts = connection.execute(query).df()
        counts["release"] = release
        frames.append(counts)
        connection.unregister("release_entities")
    support = pd.concat(frames, ignore_index = True)
    result = entities.merge(
        support, on = ["release", "source_disease_id"], how = "left", validate = "one_to_one"
    )
    result["positive_target_count"] = (
        result["positive_target_count"].fillna(0).astype(int)
    )
    return result

def build_pair_entities(entities):
    rows = []
    for (old_release, new_release) in (RELEASE_PAIRS):
        old = entities.loc[entities["release"] == old_release].set_index(
            "source_disease_id", drop = False
        )
        new = entities.loc[entities["release"] == new_release].set_index(
            "source_disease_id", drop = False
        )
        for (identifier) in (sorted(set(old.index) & set(new.index))):
            old_row = old.loc[identifier]
            new_row = new.loc[identifier]
            if (
                old_row["identifier_mapping_status"] != "EXACT_ID"
                or new_row["identifier_mapping_status"] != "EXACT_ID"
            ):
                raise RuntimeError(
                    f"Non-exact identifier entered the broader EFO pair population: {identifier}"
                )
            rows.append(
                {
                    "release_pair": f"{old_release}_to_{new_release}",
                    "old_release": old_release,
                    "new_release": new_release,
                    "source_disease_id": identifier,
                    "source_label": old_row["source_label"],
                    "old_source_therapeutic_area_ids": old_row[
                        "source_therapeutic_area_ids"
                    ],
                    "new_source_therapeutic_area_ids": new_row[
                        "source_therapeutic_area_ids"
                    ],
                    "old_canonical_therapeutic_area_ids": old_row[
                        "canonical_therapeutic_area_ids"
                    ],
                    "new_canonical_therapeutic_area_ids": new_row[
                        "canonical_therapeutic_area_ids"
                    ],
                    "label_stable": old_row["source_label"] == new_row["source_label"],
                    "canonical_disease_ta_stable": old_row[
                        "canonical_therapeutic_area_ids"
                    ]
                    == new_row["canonical_therapeutic_area_ids"],
                    "excluded_ta_context": bool(
                        old_row["excluded_ta_context"] or new_row["excluded_ta_context"]
                    ),
                }
            )
    return pd.DataFrame(rows)

def build_fixed_members(connection, project_root, pair_entities):
    frames = []
    for (old_release, new_release) in (RELEASE_PAIRS):
        release_pair = f"{old_release}_to_{new_release}"
        candidates = pair_entities.loc[
            pair_entities["release_pair"] == release_pair
        ].copy()
        connection.register("pair_entities", candidates)
        old_path = sql_path(association_path(project_root, old_release))
        new_path = sql_path(association_path(project_root, new_release))
        old_score = score_column(old_release)
        new_score = score_column(new_release)
        query = f"\n            WITH old_rows AS (\n                SELECT\n                    a.diseaseId AS source_disease_id,\n                    a.targetId AS canonical_target_id,\n                    CAST(a.{old_score} AS DOUBLE) AS old_score\n                FROM read_parquet('{old_path}') a\n                JOIN pair_entities p ON a.diseaseId = p.source_disease_id\n                WHERE isfinite(a.{old_score})\n                  AND a.{old_score} > 0\n                  AND a.evidenceCount IS NOT NULL\n                  AND a.evidenceCount >= 1\n            ), new_rows AS (\n                SELECT\n                    a.diseaseId AS source_disease_id,\n                    a.targetId AS canonical_target_id,\n                    CAST(a.{new_score} AS DOUBLE) AS new_score\n                FROM read_parquet('{new_path}') a\n                JOIN pair_entities p ON a.diseaseId = p.source_disease_id\n                WHERE isfinite(a.{new_score})\n                  AND a.{new_score} > 0\n                  AND a.evidenceCount IS NOT NULL\n                  AND a.evidenceCount >= 1\n            )\n            SELECT\n                p.release_pair,\n                p.old_release,\n                p.new_release,\n                p.source_disease_id,\n                p.source_label,\n                p.old_source_therapeutic_area_ids,\n                p.new_source_therapeutic_area_ids,\n                p.old_canonical_therapeutic_area_ids,\n                p.new_canonical_therapeutic_area_ids,\n                p.label_stable,\n                p.canonical_disease_ta_stable,\n                p.excluded_ta_context,\n                o.canonical_target_id,\n                o.old_score,\n                n.new_score\n            FROM old_rows o\n            JOIN new_rows n USING (source_disease_id, canonical_target_id)\n            JOIN pair_entities p USING (source_disease_id)\n            ORDER BY p.release_pair, p.source_disease_id, o.canonical_target_id\n        "
        frames.append(connection.execute(query).df())
        connection.unregister("pair_entities")
    result = pd.concat(frames, ignore_index = True)
    if (result.duplicated(
        ["release_pair", "source_disease_id", "canonical_target_id"]
    ).any()):
        raise RuntimeError("Broader EFO fixed support contains duplicate panel members")
    return result

def top_set(targets, scores, cutoff):
    if (len(scores) < cutoff):
        return set()
    threshold = sorted(scores, reverse = True)[cutoff - 1]
    return {
        target for ((target, score)) in (zip(targets, scores)) if (score >= threshold)
    }

def leader_set(targets, scores, epsilon = 0.0):
    maximum = max(scores)
    return {
        target
        for ((target, score)) in (zip(targets, scores))
        if (score >= maximum - epsilon)
    }

def leader_taxonomy(old_leaders, new_leaders):
    shared = old_leaders & new_leaders
    if (old_leaders == new_leaders):
        return "IDENTICAL"
    if (shared == old_leaders and len(new_leaders) > len(old_leaders)):
        return "TIE_EXPANSION"
    if (shared == new_leaders and len(old_leaders) > len(new_leaders)):
        return "TIE_CONTRACTION"
    if (shared):
        return "OVERLAPPING_LEADER_SET"
    if (len(old_leaders) == 1 and len(new_leaders) == 1):
        return "UNIQUE_LEADER_REPLACEMENT"
    return "DISJOINT_REPLACEMENT"

def jaccard(left, right):
    union = left | right
    return len(left & right) / len(union) if (union) else math.nan

def eligible_panel_groups(fixed_members, threshold):
    panel = (
        fixed_members.groupby(["release_pair", "source_disease_id"], sort = True)
        .agg(
            fixed_target_count = ("canonical_target_id", "size"),
            old_distinct_score_count = ("old_score", "nunique"),
            new_distinct_score_count = ("new_score", "nunique"),
        )
        .reset_index()
    )
    panel = panel.loc[
        (panel["fixed_target_count"] >= threshold)
        & (panel["old_distinct_score_count"] >= 2)
        & (panel["new_distinct_score_count"] >= 2)
    ]
    keys = pd.MultiIndex.from_frame(panel[["release_pair", "source_disease_id"]])
    member_keys = pd.MultiIndex.from_frame(
        fixed_members[["release_pair", "source_disease_id"]]
    )
    return fixed_members.loc[member_keys.isin(keys)].copy()

def build_panel_metrics(fixed_members):
    eligible = eligible_panel_groups(fixed_members, 30)
    rows = []
    for ((release_pair, disease), group) in (eligible.groupby(
        ["release_pair", "source_disease_id"], sort = True
    )):
        group = group.sort_values("canonical_target_id", kind = "mergesort")
        targets = group["canonical_target_id"].astype(str).tolist()
        old_scores = group["old_score"].astype(float).tolist()
        new_scores = group["new_score"].astype(float).tolist()
        old_leaders = leader_set(targets, old_scores)
        new_leaders = leader_set(targets, new_scores)
        old_top_10 = top_set(targets, old_scores, 10)
        new_top_10 = top_set(targets, new_scores, 10)
        tau = kendalltau(
            old_scores, new_scores, variant = "b", nan_policy = "raise"
        ).statistic
        absolute_change = np.abs(np.asarray(new_scores) - np.asarray(old_scores))
        first = group.iloc[0]
        rows.append(
            {
                "release_pair": release_pair,
                "old_release": first["old_release"],
                "new_release": first["new_release"],
                "source_disease_id": disease,
                "source_label": first["source_label"],
                "old_source_therapeutic_area_ids": first[
                    "old_source_therapeutic_area_ids"
                ],
                "new_source_therapeutic_area_ids": first[
                    "new_source_therapeutic_area_ids"
                ],
                "label_stable": bool(first["label_stable"]),
                "canonical_disease_ta_stable": bool(
                    first["canonical_disease_ta_stable"]
                ),
                "excluded_ta_context": bool(first["excluded_ta_context"]),
                "fixed_target_count": len(group),
                "old_distinct_score_count": len(set(old_scores)),
                "new_distinct_score_count": len(set(new_scores)),
                "median_absolute_score_change": float(
                    np.quantile(absolute_change, 0.5, method = "linear")
                ),
                "q90_absolute_score_change": float(
                    np.quantile(absolute_change, 0.9, method = "linear")
                ),
                "kendall_tau_b": float(tau),
                "old_top_ranked_target_ids": ";".join(sorted(old_leaders)),
                "new_top_ranked_target_ids": ";".join(sorted(new_leaders)),
                "top_ranked_set_changed": old_leaders != new_leaders,
                "leader_change_taxonomy": leader_taxonomy(old_leaders, new_leaders),
                "old_top_10_effective_size": len(old_top_10),
                "new_top_10_effective_size": len(new_top_10),
                "fixed_top_10_jaccard": jaccard(old_top_10, new_top_10),
                "fixed_top_10_changed": old_top_10 != new_top_10,
            }
        )
    return (pd.DataFrame(rows), eligible)

def population_frame(panel_metrics, population):
    if (population == "BROADER_EFO"):
        return panel_metrics
    return panel_metrics.loc[~panel_metrics["excluded_ta_context"]]

def build_summary(panel_metrics):
    rows = []
    for (release_pair) in ((f"{old}_to_{new}" for ((old, new)) in (RELEASE_PAIRS))):
        pair = panel_metrics.loc[panel_metrics["release_pair"] == release_pair]
        for (population) in (POPULATIONS):
            frame = population_frame(pair, population)
            changed = int(frame["top_ranked_set_changed"].sum())
            top_10_changed = int(frame["fixed_top_10_changed"].sum())
            rows.append(
                {
                    "release_pair": release_pair,
                    "population": population,
                    "minimum_fixed_target_count": 30,
                    "eligible_panel_count": len(frame),
                    "top_ranked_set_changed_count": changed,
                    "top_ranked_set_changed_fraction": changed / len(frame),
                    "fixed_top_10_changed_count": top_10_changed,
                    "fixed_top_10_changed_fraction": top_10_changed / len(frame),
                    "median_kendall_tau_b": float(frame["kendall_tau_b"].median()),
                    "median_panel_median_absolute_score_change": float(
                        frame["median_absolute_score_change"].median()
                    ),
                    "median_fixed_top_10_jaccard": float(
                        frame["fixed_top_10_jaccard"].median()
                    ),
                }
            )
    return pd.DataFrame(rows)

def build_epsilon_results(panel_metrics, fixed_members):
    eligible_keys = pd.MultiIndex.from_frame(
        panel_metrics[["release_pair", "source_disease_id"]]
    )
    member_keys = pd.MultiIndex.from_frame(
        fixed_members[["release_pair", "source_disease_id"]]
    )
    eligible = fixed_members.loc[member_keys.isin(eligible_keys)].copy()
    panel_rows = []
    for ((release_pair, disease), group) in (eligible.groupby(
        ["release_pair", "source_disease_id"], sort = True
    )):
        group = group.sort_values("canonical_target_id", kind = "mergesort")
        targets = group["canonical_target_id"].astype(str).tolist()
        old_scores = group["old_score"].astype(float).tolist()
        new_scores = group["new_score"].astype(float).tolist()
        excluded = bool(group.iloc[0]["excluded_ta_context"])
        for (epsilon) in (EPSILONS):
            old_leaders = leader_set(targets, old_scores, epsilon)
            new_leaders = leader_set(targets, new_scores, epsilon)
            panel_rows.append(
                {
                    "release_pair": release_pair,
                    "source_disease_id": disease,
                    "excluded_ta_context": excluded,
                    "epsilon": epsilon,
                    "leader_set_changed": old_leaders != new_leaders,
                    "leader_change_taxonomy": leader_taxonomy(old_leaders, new_leaders),
                }
            )
    panel_frame = pd.DataFrame(panel_rows)
    summary_rows = []
    taxonomy_rows = []
    for (release_pair) in ((f"{old}_to_{new}" for ((old, new)) in (RELEASE_PAIRS))):
        pair = panel_frame.loc[panel_frame["release_pair"] == release_pair]
        for (population) in (POPULATIONS):
            population_rows = (
                pair
                if (population == "BROADER_EFO")
                else pair.loc[~pair["excluded_ta_context"]]
            )
            for (epsilon) in (EPSILONS):
                frame = population_rows.loc[population_rows["epsilon"] == epsilon]
                changed = int(frame["leader_set_changed"].sum())
                summary_rows.append(
                    {
                        "release_pair": release_pair,
                        "population": population,
                        "minimum_fixed_target_count": 30,
                        "epsilon": epsilon,
                        "eligible_panel_count": len(frame),
                        "leader_set_changed_count": changed,
                        "leader_set_changed_fraction": changed / len(frame),
                    }
                )
                if (epsilon == 0.0):
                    counts = frame["leader_change_taxonomy"].value_counts().to_dict()
                    for (category) in (TAXONOMY):
                        count = int(counts.get(category, 0))
                        taxonomy_rows.append(
                            {
                                "release_pair": release_pair,
                                "population": population,
                                "minimum_fixed_target_count": 30,
                                "leader_change_taxonomy": category,
                                "panel_count": count,
                                "panel_fraction": count / len(frame),
                            }
                        )
    return (pd.DataFrame(summary_rows), pd.DataFrame(taxonomy_rows))

def build_support_flow(entities, pair_entities, fixed_members):
    rows = []
    for (release) in (RELEASES):
        frame = entities.loc[entities["release"] == release]
        rows.extend(
            [
                {
                    "scope": "RELEASE",
                    "release": release,
                    "release_pair": "",
                    "population": "BROADER_EFO",
                    "stage": "ELIGIBLE_ENTITY",
                    "count": len(frame),
                },
                {
                    "scope": "RELEASE",
                    "release": release,
                    "release_pair": "",
                    "population": "BROADER_EFO",
                    "stage": "POSITIVE_SUPPORT_ENTITY",
                    "count": int((frame["positive_target_count"] > 0).sum()),
                },
                {
                    "scope": "RELEASE",
                    "release": release,
                    "release_pair": "",
                    "population": "BROADER_EFO",
                    "stage": "EXCLUDED_TA_CONTEXT_ENTITY",
                    "count": int(frame["excluded_ta_context"].sum()),
                },
                {
                    "scope": "RELEASE",
                    "release": release,
                    "release_pair": "",
                    "population": "BROADER_EFO_NO_EXCLUDED_TA_CONTEXT",
                    "stage": "ELIGIBLE_ENTITY",
                    "count": int((~frame["excluded_ta_context"]).sum()),
                },
            ]
        )
    for (old_release, new_release) in (RELEASE_PAIRS):
        release_pair = f"{old_release}_to_{new_release}"
        candidates = pair_entities.loc[pair_entities["release_pair"] == release_pair]
        rows.append(
            {
                "scope": "PAIR",
                "release": "",
                "release_pair": release_pair,
                "population": "BROADER_EFO",
                "stage": "EXACT_ID_ENTITY",
                "count": len(candidates),
            }
        )
        rows.append(
            {
                "scope": "PAIR",
                "release": "",
                "release_pair": release_pair,
                "population": "BROADER_EFO_NO_EXCLUDED_TA_CONTEXT",
                "stage": "EXACT_ID_ENTITY",
                "count": int((~candidates["excluded_ta_context"]).sum()),
            }
        )
        pair_members = fixed_members.loc[fixed_members["release_pair"] == release_pair]
        for (population) in (POPULATIONS):
            members = (
                pair_members
                if (population == "BROADER_EFO")
                else pair_members.loc[~pair_members["excluded_ta_context"]]
            )
            for (threshold) in (PANEL_THRESHOLDS):
                eligible = eligible_panel_groups(members, threshold)
                count = (
                    eligible[["release_pair", "source_disease_id"]]
                    .drop_duplicates()
                    .shape[0]
                )
                rows.append(
                    {
                        "scope": "PAIR",
                        "release": "",
                        "release_pair": release_pair,
                        "population": population,
                        "stage": f"FIXED_PANEL_N_GE_{threshold}",
                        "count": count,
                    }
                )
    exact_all = set.intersection(
        *(
            set(entities.loc[entities["release"] == release, "source_disease_id"])
            for (release) in (RELEASES)
        )
    )
    rows.append(
        {
            "scope": "THREE_RELEASE",
            "release": "",
            "release_pair": "",
            "population": "BROADER_EFO",
            "stage": "EXACT_ID_ENTITY",
            "count": len(exact_all),
        }
    )
    return pd.DataFrame(rows)

def run(project_root, output_dir = None, threads = 8, memory_limit = "8GB"):
    args = SimpleNamespace(
        project_root = Path(project_root),
        output_dir = output_dir,
        threads = threads,
        memory_limit = memory_limit,
    )
    project_root = args.project_root.resolve()
    output_dir = (
        args.output_dir or project_root / "data" / "derived" / "robustness"
    ).resolve()
    output_dir.mkdir(parents = True, exist_ok = True)
    connection = configure_connection(project_root, args.threads, args.memory_limit)
    try:
        entities = load_broader_entities(project_root)
        entities = add_positive_support_counts(connection, project_root, entities)
        pair_entities = build_pair_entities(entities)
        fixed_members = build_fixed_members(connection, project_root, pair_entities)
        panel_metrics, eligible_members = build_panel_metrics(fixed_members)
        summary = build_summary(panel_metrics)
        epsilon_summary, taxonomy_summary = build_epsilon_results(
            panel_metrics, fixed_members
        )
        support_flow = build_support_flow(entities, pair_entities, fixed_members)
    finally:
        connection.close()
    entity_columns = [
        "release",
        "source_disease_id",
        "source_label",
        "source_therapeutic_area_ids",
        "canonical_disease_id",
        "canonical_therapeutic_area_ids",
        "identifier_mapping_status",
        "ambiguity_status",
        "exclusion_reason",
        "excluded_ta_context_ids",
        "excluded_ta_context",
        "positive_target_count",
    ]
    member_columns = [
        "release_pair",
        "old_release",
        "new_release",
        "source_disease_id",
        "source_label",
        "old_source_therapeutic_area_ids",
        "new_source_therapeutic_area_ids",
        "old_canonical_therapeutic_area_ids",
        "new_canonical_therapeutic_area_ids",
        "label_stable",
        "canonical_disease_ta_stable",
        "excluded_ta_context",
        "canonical_target_id",
        "old_score",
        "new_score",
    ]
    outputs = {
        "broader_efo_entity_registry.csv": (
            entities[entity_columns],
            ["release", "source_disease_id"],
        ),
        "broader_efo_support_flow.csv": (
            support_flow,
            ["scope", "release", "release_pair", "population", "stage"],
        ),
        "broader_efo_fixed_panel_members.csv": (
            eligible_members[member_columns],
            ["release_pair", "source_disease_id", "canonical_target_id"],
        ),
        "broader_efo_panel_metrics.csv": (
            panel_metrics,
            ["release_pair", "source_disease_id"],
        ),
        "broader_efo_summary.csv": (summary, ["release_pair", "population"]),
        "broader_efo_epsilon_leader_sensitivity.csv": (
            epsilon_summary,
            ["release_pair", "population", "epsilon"],
        ),
        "broader_efo_exact_leader_change_taxonomy.csv": (
            taxonomy_summary,
            ["release_pair", "population", "leader_change_taxonomy"],
        ),
    }
    record_counts = {}
    for (name, (frame, order)) in (outputs.items()):
        record_counts[name] = write_frame(frame, output_dir / name, order)
