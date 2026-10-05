from __future__ import annotations
import math
import shutil
import duckdb
import numpy as np
import pandas as pd
from scipy.stats import kendalltau
from scipy.stats import rankdata
from scipy.stats import spearmanr
from .robustness_types import EPSILONS, MARGIN_LABELS, TAXONOMY

def sql_path(path):
    return path.resolve().as_posix().replace("'", "''")

def write_frame(frame, path, order, written):
    result = frame.copy()
    if (order and (not result.empty)):
        result = result.sort_values(order, kind = "mergesort").reset_index(drop = True)
    result.to_csv(
        path, index = False, lineterminator = "\n", na_rep = "", float_format = "%.17g"
    )
    written.append((path, len(result)))

def configure_connection(project_root, threads = 4, memory_limit = "8GB"):
    temporary = project_root / "data" / "interim" / "robustness_duckdb_tmp"
    temporary.mkdir(parents = True, exist_ok = True)
    connection = duckdb.connect()
    connection.execute(f"SET threads = {int(threads)}")
    connection.execute(
        "SET memory_limit = '" + str(memory_limit).replace("'", "''") + "'"
    )
    connection.execute(f"SET temp_directory = '{sql_path(temporary)}'")
    return (connection, temporary)

def create_views(connection, derived):
    paths = {
        "fixed_members": derived / "fixed_ranking_panel_members.parquet",
        "primary_metrics": derived / "ranking_sensitivity_fixed_ranking_metrics.csv",
        "native_metrics": derived
        / "ranking_sensitivity_release_native_top_set_metrics.csv",
        "top_members": derived / "ranking_sensitivity_top_set_members.parquet",
        "stable_members": derived
        / "cross_layer_sensitivity_stable_entity_pair_revisions.parquet",
        "stable_metrics": derived
        / "cross_layer_sensitivity_stable_entity_ranking_metrics.csv",
        "genetic_members": derived
        / "cross_layer_sensitivity_evidence_domain_pair_revisions.parquet",
        "genetic_metrics": derived
        / "cross_layer_sensitivity_evidence_domain_ranking_metrics.csv",
        "disease_crosswalk": derived / "disease_crosswalk.csv",
        "target_crosswalk": derived / "target_crosswalk.csv",
    }
    for (path) in (paths.values()):
        if (not path.is_file()):
            raise FileNotFoundError(path)
    connection.execute(
        f"CREATE VIEW primary_metric_rows AS SELECT * FROM read_csv_auto('{sql_path(paths['primary_metrics'])}', header = true)"
    )
    connection.execute(
        f"CREATE VIEW native_metric_rows AS SELECT * FROM read_csv_auto('{sql_path(paths['native_metrics'])}', header = true)"
    )
    connection.execute(
        f"CREATE VIEW stable_metric_rows AS SELECT * FROM read_csv_auto('{sql_path(paths['stable_metrics'])}', header = true)"
    )
    connection.execute(
        f"CREATE VIEW genetic_metric_rows AS SELECT * FROM read_csv_auto('{sql_path(paths['genetic_metrics'])}', header = true)"
    )
    connection.execute(
        f"CREATE VIEW fixed_member_source AS SELECT * FROM read_parquet('{sql_path(paths['fixed_members'])}')"
    )
    connection.execute(
        f"CREATE VIEW top_member_rows AS SELECT * FROM read_parquet('{sql_path(paths['top_members'])}')"
    )
    connection.execute(
        f"CREATE VIEW stable_member_rows AS SELECT release_pair, canonical_disease_id, canonical_target_id, CAST(old_score AS DOUBLE) AS old_score, CAST(new_score AS DOUBLE) AS new_score FROM read_parquet('{sql_path(paths['stable_members'])}')"
    )
    connection.execute(
        f"CREATE VIEW genetic_member_rows AS SELECT release_pair, canonical_disease_id, canonical_target_id, CAST(old_score AS DOUBLE) AS old_score, CAST(new_score AS DOUBLE) AS new_score FROM read_parquet('{sql_path(paths['genetic_members'])}')"
    )
    connection.execute(
        f"CREATE VIEW disease_crosswalk_rows AS SELECT * FROM read_csv_auto('{sql_path(paths['disease_crosswalk'])}', header = true, all_varchar = true)"
    )
    connection.execute(
        f"CREATE VIEW target_crosswalk_rows AS SELECT * FROM read_csv_auto('{sql_path(paths['target_crosswalk'])}', header = true, all_varchar = true)"
    )
    connection.execute(
        "\n        CREATE VIEW primary_member_rows AS\n        SELECT\n            m.release_pair,\n            m.canonical_disease_id,\n            m.canonical_target_id,\n            CAST(m.old_score AS DOUBLE) AS old_score,\n            CAST(m.new_score AS DOUBLE) AS new_score,\n            CAST(m.stable_entity_member AS BOOLEAN) AS stable_entity_member\n        FROM fixed_member_source m\n        JOIN primary_metric_rows p USING (release_pair, canonical_disease_id)\n        WHERE m.support_name = 'PAIR_SPECIFIC_PERSISTENT'\n          AND p.panel_population = 'PRIMARY_FIXED'\n          AND CAST(p.minimum_target_threshold AS BIGINT) = 30\n        "
    )
    return paths

def layer_specifications():
    return {
        "PRIMARY_FIXED": {
            "members": "primary_member_rows",
            "metrics": "primary_metric_rows",
            "metric_filter": "panel_population = 'PRIMARY_FIXED' AND CAST(minimum_target_threshold AS BIGINT) = 30",
            "old_ids": "old_fixed_top_ranked_target_ids",
            "new_ids": "new_fixed_top_ranked_target_ids",
            "changed": "fixed_top_ranked_set_changed",
        },
        "STABLE_ENTITY_FIXED": {
            "members": "stable_member_rows",
            "metrics": "stable_metric_rows",
            "metric_filter": "TRUE",
            "old_ids": "old_top_ranked_target_ids",
            "new_ids": "new_top_ranked_target_ids",
            "changed": "top_ranked_set_changed",
        },
        "GENETIC_ASSOCIATION_FIXED": {
            "members": "genetic_member_rows",
            "metrics": "genetic_metric_rows",
            "metric_filter": "TRUE",
            "old_ids": "old_top_ranked_target_ids",
            "new_ids": "new_top_ranked_target_ids",
            "changed": "top_ranked_set_changed",
        },
    }

def epsilon_panel_frame(connection, relation, layer):
    epsilon_values = ",".join((f"({value:.17g})" for (value) in (EPSILONS)))
    query = f"\n        WITH epsilon_values(epsilon) AS (\n            VALUES {epsilon_values}\n        ), maxima AS (\n            SELECT\n                release_pair,\n                canonical_disease_id,\n                MAX(old_score) AS old_maximum,\n                MAX(new_score) AS new_maximum\n            FROM {relation}\n            GROUP BY release_pair, canonical_disease_id\n        ), leaders AS (\n            SELECT\n                r.release_pair,\n                r.canonical_disease_id,\n                r.canonical_target_id,\n                e.epsilon,\n                r.old_score >= m.old_maximum - e.epsilon AS old_leader,\n                r.new_score >= m.new_maximum - e.epsilon AS new_leader\n            FROM {relation} r\n            JOIN maxima m USING (release_pair, canonical_disease_id)\n            CROSS JOIN epsilon_values e\n            WHERE r.old_score >= m.old_maximum - e.epsilon\n               OR r.new_score >= m.new_maximum - e.epsilon\n        ), aggregated AS (\n            SELECT\n                release_pair,\n                canonical_disease_id,\n                epsilon,\n                COUNT(*) FILTER (WHERE old_leader) AS old_leader_count,\n                COUNT(*) FILTER (WHERE new_leader) AS new_leader_count,\n                COUNT(*) FILTER (WHERE old_leader AND new_leader) AS shared_leader_count,\n                STRING_AGG(canonical_target_id, ';' ORDER BY canonical_target_id) FILTER (WHERE old_leader) AS old_leader_target_ids,\n                STRING_AGG(canonical_target_id, ';' ORDER BY canonical_target_id) FILTER (WHERE new_leader) AS new_leader_target_ids\n            FROM leaders\n            GROUP BY release_pair, canonical_disease_id, epsilon\n        )\n        SELECT\n            '{layer}' AS analysis_layer,\n            release_pair,\n            canonical_disease_id,\n            epsilon,\n            old_leader_target_ids,\n            new_leader_target_ids,\n            old_leader_count,\n            new_leader_count,\n            shared_leader_count,\n            old_leader_count <> new_leader_count OR shared_leader_count <> old_leader_count AS leader_set_changed,\n            shared_leader_count::DOUBLE / (old_leader_count + new_leader_count - shared_leader_count) AS leader_set_jaccard,\n            CASE\n                WHEN old_leader_count = new_leader_count AND shared_leader_count = old_leader_count THEN 'IDENTICAL'\n                WHEN shared_leader_count = old_leader_count AND new_leader_count > old_leader_count THEN 'TIE_EXPANSION'\n                WHEN shared_leader_count = new_leader_count AND old_leader_count > new_leader_count THEN 'TIE_CONTRACTION'\n                WHEN old_leader_count = 1 AND new_leader_count = 1 AND shared_leader_count = 0 THEN 'UNIQUE_LEADER_REPLACEMENT'\n                WHEN shared_leader_count = 0 THEN 'DISJOINT_REPLACEMENT'\n                ELSE 'OVERLAPPING_LEADER_SET'\n            END AS leader_change_taxonomy\n        FROM aggregated\n        ORDER BY release_pair, canonical_disease_id, epsilon\n    "
    return connection.sql(query).df()

def build_epsilon_outputs(connection, specifications, output_dir, written):
    frames = []
    validation_rows = []
    exact_frames = {}
    for (layer, specification) in (specifications.items()):
        frame = epsilon_panel_frame(connection, specification["members"], layer)
        frames.append(frame)
        exact = frame.loc[frame["epsilon"] == 0.0].copy()
        exact_frames[layer] = exact
        expected = connection.sql(
            f"SELECT release_pair, canonical_disease_id, {specification['old_ids']} AS expected_old_ids, {specification['new_ids']} AS expected_new_ids, CAST({specification['changed']} AS BOOLEAN) AS expected_changed FROM {specification['metrics']} WHERE {specification['metric_filter']}"
        ).df()
        merged = exact.merge(
            expected,
            on = ["release_pair", "canonical_disease_id"],
            how = "outer",
            validate = "one_to_one",
        )
        old_match = (
            merged["old_leader_target_ids"]
            .fillna("")
            .eq(merged["expected_old_ids"].fillna(""))
        )
        new_match = (
            merged["new_leader_target_ids"]
            .fillna("")
            .eq(merged["expected_new_ids"].fillna(""))
        )
        change_match = merged["leader_set_changed"].eq(merged["expected_changed"])
        for (release_pair, group) in (merged.assign(
            match = old_match & new_match & change_match
        ).groupby("release_pair", dropna = False)):
            validation_rows.append(
                {
                    "analysis_layer": layer,
                    "release_pair": release_pair,
                    "panel_count": len(group),
                    "mismatch_count": int((~group["match"]).sum()),
                    "passed": bool(group["match"].all()),
                }
            )
    epsilon = pd.concat(frames, ignore_index = True)
    exact = epsilon.loc[epsilon["epsilon"] == 0.0].drop(columns = ["epsilon"]).copy()
    write_frame(
        exact,
        output_dir / "exact_leader_change_taxonomy.csv",
        ["analysis_layer", "release_pair", "canonical_disease_id"],
        written,
    )
    summary_rows = []
    for ((layer, release_pair, epsilon_value), group) in (epsilon.groupby(
        ["analysis_layer", "release_pair", "epsilon"], sort = True
    )):
        row = {
            "analysis_layer": layer,
            "release_pair": release_pair,
            "epsilon": epsilon_value,
            "eligible_panel_count": len(group),
            "changed_leader_set_count": int(group["leader_set_changed"].sum()),
            "changed_leader_set_fraction": float(group["leader_set_changed"].mean()),
            "median_leader_set_jaccard": float(group["leader_set_jaccard"].median()),
            "median_old_leader_set_size": float(group["old_leader_count"].median()),
            "median_new_leader_set_size": float(group["new_leader_count"].median()),
        }
        for (category) in (TAXONOMY):
            row[f"{category.casefold()}_count"] = int(
                (group["leader_change_taxonomy"] == category).sum()
            )
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    write_frame(
        summary,
        output_dir / "epsilon_leader_sensitivity.csv",
        ["analysis_layer", "release_pair", "epsilon"],
        written,
    )
    validation = pd.DataFrame(validation_rows)
    write_frame(
        validation,
        output_dir / "exact_leader_reproduction_validation.csv",
        ["analysis_layer", "release_pair"],
        written,
    )
    if (not validation["passed"].all()):
        raise RuntimeError("Exact leader sets do not reproduce the frozen metrics")
    return exact_frames

def assign_margin_bins(values):
    bins = [-np.inf, 0.0, 1e-12, 1e-08, 1e-06, 0.0001, 0.001, 0.01, 0.05, np.inf]
    return pd.cut(
        values, bins = bins, labels = MARGIN_LABELS, right = True, include_lowest = True
    )

def build_margin_outputs(connection, specifications, exact_frames, output_dir, written):
    panel_frames = []
    for (layer, specification) in (specifications.items()):
        relation = specification["members"]
        margins = connection.sql(
            f"\n            WITH maxima AS (\n                SELECT\n                    release_pair,\n                    canonical_disease_id,\n                    MAX(old_score) AS old_top_score,\n                    MAX(new_score) AS new_top_score\n                FROM {relation}\n                GROUP BY release_pair, canonical_disease_id\n            ), details AS (\n                SELECT\n                    r.release_pair,\n                    r.canonical_disease_id,\n                    MAX(m.old_top_score) AS old_top_score,\n                    MAX(m.new_top_score) AS new_top_score,\n                    COUNT(*) FILTER (WHERE r.old_score = m.old_top_score) AS old_exact_leader_count,\n                    COUNT(*) FILTER (WHERE r.new_score = m.new_top_score) AS new_exact_leader_count,\n                    MAX(r.old_score) FILTER (WHERE r.old_score < m.old_top_score) AS old_second_distinct_score,\n                    MAX(r.new_score) FILTER (WHERE r.new_score < m.new_top_score) AS new_second_distinct_score\n                FROM {relation} r\n                JOIN maxima m USING (release_pair, canonical_disease_id)\n                GROUP BY r.release_pair, r.canonical_disease_id\n            )\n            SELECT\n                '{layer}' AS analysis_layer,\n                *,\n                CASE WHEN old_exact_leader_count > 1 THEN 0.0 ELSE old_top_score - old_second_distinct_score END AS old_top_1_to_2_margin,\n                CASE WHEN new_exact_leader_count > 1 THEN 0.0 ELSE new_top_score - new_second_distinct_score END AS new_top_1_to_2_margin\n            FROM details\n            ORDER BY release_pair, canonical_disease_id\n            "
        ).df()
        exact = exact_frames[layer][
            [
                "release_pair",
                "canonical_disease_id",
                "leader_set_changed",
                "leader_change_taxonomy",
            ]
        ]
        margins = margins.merge(
            exact, on = ["release_pair", "canonical_disease_id"], validate = "one_to_one"
        )
        margins["old_margin_bin"] = assign_margin_bins(margins["old_top_1_to_2_margin"])
        panel_frames.append(margins)
    panels = pd.concat(panel_frames, ignore_index = True)
    write_frame(
        panels,
        output_dir / "leader_margin_panels.csv",
        ["analysis_layer", "release_pair", "canonical_disease_id"],
        written,
    )
    summary_rows = []
    for ((layer, release_pair, margin_bin), group) in (panels.groupby(
        ["analysis_layer", "release_pair", "old_margin_bin"], observed = True, sort = True
    )):
        summary_rows.append(
            {
                "analysis_layer": layer,
                "release_pair": release_pair,
                "old_margin_bin": str(margin_bin),
                "eligible_panel_count": len(group),
                "changed_leader_set_count": int(group["leader_set_changed"].sum()),
                "changed_leader_set_fraction": float(
                    group["leader_set_changed"].mean()
                ),
                "old_tied_leader_count": int(
                    (group["old_exact_leader_count"] > 1).sum()
                ),
                "old_tied_leader_fraction": float(
                    (group["old_exact_leader_count"] > 1).mean()
                ),
                "median_old_top_1_to_2_margin": float(
                    group["old_top_1_to_2_margin"].median()
                ),
                "p90_old_top_1_to_2_margin": float(
                    group["old_top_1_to_2_margin"].quantile(0.9)
                ),
            }
        )
    write_frame(
        pd.DataFrame(summary_rows),
        output_dir / "leader_margin_summary.csv",
        ["analysis_layer", "release_pair", "old_margin_bin"],
        written,
    )
    return panels

def build_endpoint_matrix(paths, output_dir, written):
    metrics = pd.read_csv(paths["primary_metrics"])
    rows = []

    def add(
        group, threshold, release_pair, endpoint, numerator, denominator, value, unit
    ):
        rows.append(
            {
                "analysis_layer": "PRIMARY_FIXED",
                "support_definition": "pair-specific persistent disease-target pairs",
                "minimum_panel_size": threshold,
                "analysis_status": "primary" if (threshold == 30) else "sensitivity",
                "release_pair": release_pair,
                "endpoint": endpoint,
                "numerator": numerator,
                "denominator": denominator,
                "value": value,
                "unit": unit,
                "source_file": "data/derived/ranking_sensitivity_fixed_ranking_metrics.csv",
            }
        )

    for ((threshold, release_pair), group) in (metrics.groupby(
        ["minimum_target_threshold", "release_pair"], sort = True
    )):
        threshold = int(threshold)
        count = len(group)
        add(
            group,
            threshold,
            release_pair,
            "eligible_panel_count",
            count,
            count,
            count,
            "disease panels",
        )
        changed = group["fixed_top_ranked_set_changed"].astype(bool)
        mapping_only = (
            group["fixed_top_ranked_change_reason"].eq("MAPPING_AFFECTED") & changed
        )
        non_mapping = changed & ~mapping_only
        add(
            group,
            threshold,
            release_pair,
            "top_ranked_set_changed_raw",
            int(changed.sum()),
            count,
            float(changed.mean()),
            "fraction",
        )
        add(
            group,
            threshold,
            release_pair,
            "top_ranked_set_changed_mapping_only",
            int(mapping_only.sum()),
            count,
            float(mapping_only.mean()),
            "fraction",
        )
        add(
            group,
            threshold,
            release_pair,
            "top_ranked_set_changed_non_mapping_only",
            int(non_mapping.sum()),
            count,
            float(non_mapping.mean()),
            "fraction",
        )
        for (cutoff) in ([5, 10, 20]):
            estimable = group[f"fixed_top_{cutoff}_changed"].notna()
            values = group.loc[estimable, f"fixed_top_{cutoff}_changed"].astype(bool)
            denominator = int(estimable.sum())
            add(
                group,
                threshold,
                release_pair,
                f"top_{cutoff}_set_changed",
                int(values.sum()),
                denominator,
                float(values.mean()),
                "fraction",
            )
            jaccard = group.loc[estimable, f"fixed_top_{cutoff}_jaccard"]
            add(
                group,
                threshold,
                release_pair,
                f"median_top_{cutoff}_jaccard",
                pd.NA,
                denominator,
                float(jaccard.median()),
                "Jaccard index",
            )
        add(
            group,
            threshold,
            release_pair,
            "median_panel_kendall_tau_b",
            pd.NA,
            count,
            float(group["kendall_tau_b"].median()),
            "Kendall tau-b",
        )
        add(
            group,
            threshold,
            release_pair,
            "median_panel_p90_normalised_rank_displacement",
            pd.NA,
            count,
            float(group["p90_normalised_rank_displacement"].median()),
            "normalised rank displacement",
        )
        high = group["p90_normalised_rank_displacement"] >= 0.1
        add(
            group,
            threshold,
            release_pair,
            "panels_with_p90_normalised_rank_displacement_at_least_0.10",
            int(high.sum()),
            count,
            float(high.mean()),
            "fraction",
        )
    frame = pd.DataFrame(rows)
    write_frame(
        frame,
        output_dir / "panel_size_endpoint_matrix.csv",
        ["minimum_panel_size", "release_pair", "endpoint"],
        written,
    )

def build_native_distribution(connection, output_dir, written):
    panels = connection.sql(
        "\n        WITH members AS (\n            SELECT\n                release_pair,\n                canonical_disease_id,\n                COUNT(*) FILTER (WHERE CAST(old_member AS BOOLEAN)) AS old_top_10_size,\n                COUNT(*) FILTER (WHERE CAST(new_member AS BOOLEAN)) AS new_top_10_size,\n                COUNT(*) FILTER (WHERE CAST(old_member AS BOOLEAN) AND CAST(new_member AS BOOLEAN)) AS shared_member_count,\n                COUNT(*) FILTER (WHERE CAST(old_member AS BOOLEAN) AND NOT CAST(new_member AS BOOLEAN)) AS old_lost_count,\n                COUNT(*) FILTER (WHERE CAST(new_member AS BOOLEAN) AND NOT CAST(old_member AS BOOLEAN)) AS new_added_count\n            FROM top_member_rows\n            WHERE analysis_type = 'RELEASE_NATIVE'\n              AND set_definition = 'TOP_10'\n            GROUP BY release_pair, canonical_disease_id\n        ), joined AS (\n            SELECT\n                m.*,\n                n.primary_fixed_anchored_eligible,\n                n.all_native_supplementary_eligible,\n                CAST(n.native_top_10_jaccard AS DOUBLE) AS reported_jaccard,\n                CAST(n.old_native_top_10_effective_size AS BIGINT) AS reported_old_size,\n                CAST(n.new_native_top_10_effective_size AS BIGINT) AS reported_new_size,\n                CAST(n.native_top_10_changed AS BOOLEAN) AS reported_changed\n            FROM members m\n            JOIN native_metric_rows n USING (release_pair, canonical_disease_id)\n        ), populations AS (\n            SELECT 'PRIMARY_FIXED_ANCHORED' AS population, * FROM joined WHERE CAST(primary_fixed_anchored_eligible AS BOOLEAN)\n            UNION ALL\n            SELECT 'ALL_NATIVE_ELIGIBLE' AS population, * FROM joined WHERE CAST(all_native_supplementary_eligible AS BOOLEAN)\n        )\n        SELECT\n            population,\n            release_pair,\n            canonical_disease_id,\n            old_top_10_size,\n            new_top_10_size,\n            shared_member_count,\n            old_lost_count,\n            new_added_count,\n            GREATEST(old_lost_count, new_added_count) AS replacement_count,\n            old_lost_count + new_added_count AS total_membership_change_count,\n            shared_member_count::DOUBLE / NULLIF(old_top_10_size, 0) AS retained_old_fraction,\n            shared_member_count::DOUBLE / NULLIF(old_top_10_size + new_top_10_size - shared_member_count, 0) AS calculated_jaccard,\n            old_lost_count = 0 AND new_added_count = 0 AS exact_unchanged,\n            reported_jaccard,\n            reported_old_size,\n            reported_new_size,\n            reported_changed\n        FROM populations\n        ORDER BY population, release_pair, canonical_disease_id\n        "
    ).df()
    panels["validation_passed"] = (
        panels["old_top_10_size"].eq(panels["reported_old_size"])
        & panels["new_top_10_size"].eq(panels["reported_new_size"])
        & np.isclose(
            panels["calculated_jaccard"],
            panels["reported_jaccard"],
            atol = 1e-12,
            rtol = 0.0,
        )
        & panels["exact_unchanged"].eq(~panels["reported_changed"])
    )
    write_frame(
        panels,
        output_dir / "release_native_top10_panels.csv",
        ["population", "release_pair", "canonical_disease_id"],
        written,
    )
    distribution_rows = []
    distribution_fields = [
        "shared_member_count",
        "old_lost_count",
        "new_added_count",
        "replacement_count",
        "total_membership_change_count",
    ]
    for ((population, release_pair), group) in (panels.groupby(
        ["population", "release_pair"], sort = True
    )):
        for (field) in (distribution_fields):
            counts = group[field].value_counts().sort_index()
            for (value, count) in (counts.items()):
                distribution_rows.append(
                    {
                        "population": population,
                        "release_pair": release_pair,
                        "metric": field,
                        "value": int(value),
                        "panel_count": int(count),
                        "panel_fraction": float(count / len(group)),
                    }
                )
    write_frame(
        pd.DataFrame(distribution_rows),
        output_dir / "release_native_top10_shared_replacement_distribution.csv",
        ["population", "release_pair", "metric", "value"],
        written,
    )
    summary_rows = []
    for ((population, release_pair), group) in (panels.groupby(
        ["population", "release_pair"], sort = True
    )):
        summary_rows.append(
            {
                "population": population,
                "release_pair": release_pair,
                "eligible_panel_count": len(group),
                "exact_unchanged_count": int(group["exact_unchanged"].sum()),
                "exact_unchanged_fraction": float(group["exact_unchanged"].mean()),
                "median_shared_member_count": float(
                    group["shared_member_count"].median()
                ),
                "p10_shared_member_count": float(
                    group["shared_member_count"].quantile(0.1)
                ),
                "p90_shared_member_count": float(
                    group["shared_member_count"].quantile(0.9)
                ),
                "median_replacement_count": float(group["replacement_count"].median()),
                "p90_replacement_count": float(
                    group["replacement_count"].quantile(0.9)
                ),
                "mean_jaccard": float(group["calculated_jaccard"].mean()),
                "median_jaccard": float(group["calculated_jaccard"].median()),
                "validation_failure_count": int((~group["validation_passed"]).sum()),
            }
        )
    summary = pd.DataFrame(summary_rows)
    write_frame(
        summary,
        output_dir / "release_native_top10_summary.csv",
        ["population", "release_pair"],
        written,
    )
    if (not panels["validation_passed"].all()):
        raise RuntimeError(
            "Release-native top-10 membership does not reproduce frozen metrics"
        )

def build_genetic_diagnostics(connection, margin_panels, output_dir, written):
    genetic = connection.sql(
        "SELECT * FROM genetic_metric_rows ORDER BY release_pair, canonical_disease_id"
    ).df()
    primary = connection.sql(
        "\n        SELECT\n            release_pair,\n            canonical_disease_id,\n            CAST(fixed_target_count AS BIGINT) AS primary_fixed_target_count,\n            CAST(fixed_top_ranked_set_changed AS BOOLEAN) AS primary_top_ranked_set_changed,\n            CAST(fixed_top_10_changed AS BOOLEAN) AS primary_top_10_changed,\n            CAST(kendall_tau_b AS DOUBLE) AS primary_kendall_tau_b\n        FROM primary_metric_rows\n        WHERE panel_population = 'PRIMARY_FIXED'\n          AND CAST(minimum_target_threshold AS BIGINT) = 30\n        "
    ).df()
    cutoffs = connection.sql(
        "\n        WITH positioned AS (\n            SELECT\n                release_pair,\n                canonical_disease_id,\n                old_score,\n                new_score,\n                ROW_NUMBER() OVER (PARTITION BY release_pair, canonical_disease_id ORDER BY old_score DESC, canonical_target_id) AS old_position,\n                ROW_NUMBER() OVER (PARTITION BY release_pair, canonical_disease_id ORDER BY new_score DESC, canonical_target_id) AS new_position\n            FROM genetic_member_rows\n        )\n        SELECT\n            release_pair,\n            canonical_disease_id,\n            MAX(old_score) FILTER (WHERE old_position = 10) AS old_top_10_cutoff_score,\n            MAX(old_score) FILTER (WHERE old_position = 11) AS old_top_10_next_score,\n            MAX(new_score) FILTER (WHERE new_position = 10) AS new_top_10_cutoff_score,\n            MAX(new_score) FILTER (WHERE new_position = 11) AS new_top_10_next_score\n        FROM positioned\n        GROUP BY release_pair, canonical_disease_id\n        "
    ).df()
    margins = margin_panels.loc[
        margin_panels["analysis_layer"] == "GENETIC_ASSOCIATION_FIXED",
        [
            "release_pair",
            "canonical_disease_id",
            "old_top_1_to_2_margin",
            "new_top_1_to_2_margin",
        ],
    ]
    frame = genetic.merge(
        primary,
        on = ["release_pair", "canonical_disease_id"],
        how = "left",
        validate = "one_to_one",
    )
    frame = frame.merge(
        cutoffs, on = ["release_pair", "canonical_disease_id"], validate = "one_to_one"
    )
    frame = frame.merge(
        margins, on = ["release_pair", "canonical_disease_id"], validate = "one_to_one"
    )
    frame["genetic_to_primary_target_ratio"] = (
        frame["fixed_target_count"] / frame["primary_fixed_target_count"]
    )
    frame["old_positive_score_fraction"] = (
        frame["old_positive_score_count"] / frame["fixed_target_count"]
    )
    frame["new_positive_score_fraction"] = (
        frame["new_positive_score_count"] / frame["fixed_target_count"]
    )
    frame["old_score_collision_fraction"] = (
        1.0 - frame["old_distinct_score_count"] / frame["fixed_target_count"]
    )
    frame["new_score_collision_fraction"] = (
        1.0 - frame["new_distinct_score_count"] / frame["fixed_target_count"]
    )
    frame["old_leader_tie"] = frame["old_top_ranked_set_size"] > 1
    frame["new_leader_tie"] = frame["new_top_ranked_set_size"] > 1
    frame["old_top_10_cutoff_margin"] = (
        frame["old_top_10_cutoff_score"] - frame["old_top_10_next_score"]
    )
    frame["new_top_10_cutoff_margin"] = (
        frame["new_top_10_cutoff_score"] - frame["new_top_10_next_score"]
    )
    frame["old_top_10_boundary_tie"] = frame["old_top_10_cutoff_margin"] == 0.0
    frame["new_top_10_boundary_tie"] = frame["new_top_10_cutoff_margin"] == 0.0
    selected = frame[
        [
            "release_pair",
            "canonical_disease_id",
            "fixed_target_count",
            "primary_fixed_target_count",
            "genetic_to_primary_target_ratio",
            "old_positive_score_fraction",
            "new_positive_score_fraction",
            "old_score_collision_fraction",
            "new_score_collision_fraction",
            "old_leader_tie",
            "new_leader_tie",
            "old_top_1_to_2_margin",
            "new_top_1_to_2_margin",
            "old_top_10_effective_size",
            "new_top_10_effective_size",
            "old_top_10_cutoff_margin",
            "new_top_10_cutoff_margin",
            "old_top_10_boundary_tie",
            "new_top_10_boundary_tie",
            "kendall_tau_b",
            "top_ranked_set_changed",
            "top_10_changed",
            "primary_kendall_tau_b",
            "primary_top_ranked_set_changed",
            "primary_top_10_changed",
        ]
    ].copy()
    write_frame(
        selected,
        output_dir / "genetic_coverage_tie_diagnostics.csv",
        ["release_pair", "canonical_disease_id"],
        written,
    )
    summary_rows = []
    for (release_pair, group) in (selected.groupby("release_pair", sort = True)):
        summary_rows.append(
            {
                "release_pair": release_pair,
                "eligible_genetic_panel_count": len(group),
                "matched_primary_panel_count": int(
                    group["primary_fixed_target_count"].notna().sum()
                ),
                "median_genetic_target_count": float(
                    group["fixed_target_count"].median()
                ),
                "median_primary_target_count": float(
                    group["primary_fixed_target_count"].median()
                ),
                "median_genetic_to_primary_target_ratio": float(
                    group["genetic_to_primary_target_ratio"].median()
                ),
                "median_old_positive_score_fraction": float(
                    group["old_positive_score_fraction"].median()
                ),
                "median_new_positive_score_fraction": float(
                    group["new_positive_score_fraction"].median()
                ),
                "old_leader_tie_fraction": float(group["old_leader_tie"].mean()),
                "new_leader_tie_fraction": float(group["new_leader_tie"].mean()),
                "old_top_10_boundary_tie_fraction": float(
                    group["old_top_10_boundary_tie"].mean()
                ),
                "new_top_10_boundary_tie_fraction": float(
                    group["new_top_10_boundary_tie"].mean()
                ),
                "median_old_top_1_to_2_margin": float(
                    group["old_top_1_to_2_margin"].median()
                ),
                "median_new_top_1_to_2_margin": float(
                    group["new_top_1_to_2_margin"].median()
                ),
                "median_old_top_10_cutoff_margin": float(
                    group["old_top_10_cutoff_margin"].median()
                ),
                "median_new_top_10_cutoff_margin": float(
                    group["new_top_10_cutoff_margin"].median()
                ),
                "median_genetic_kendall_tau_b": float(group["kendall_tau_b"].median()),
                "genetic_top_ranked_change_fraction": float(
                    group["top_ranked_set_changed"].astype(bool).mean()
                ),
                "genetic_top_10_change_fraction": float(
                    group["top_10_changed"].astype(bool).mean()
                ),
            }
        )
    write_frame(
        pd.DataFrame(summary_rows),
        output_dir / "genetic_coverage_tie_summary.csv",
        ["release_pair"],
        written,
    )
    concordance = (
        selected.groupby(
            [
                "release_pair",
                "primary_top_ranked_set_changed",
                "top_ranked_set_changed",
            ],
            dropna = False,
            sort = True,
        )
        .size()
        .reset_index(name = "panel_count")
        .rename(columns = {"top_ranked_set_changed": "genetic_top_ranked_set_changed"})
    )
    concordance["panel_fraction_within_release_pair"] = concordance[
        "panel_count"
    ] / concordance.groupby("release_pair")["panel_count"].transform("sum")
    write_frame(
        concordance,
        output_dir / "genetic_primary_leader_change_concordance.csv",
        [
            "release_pair",
            "primary_top_ranked_set_changed",
            "genetic_top_ranked_set_changed",
        ],
        written,
    )

def build_ablation_outputs(connection, output_dir, written):
    connection.execute(
        "\n        CREATE VIEW stable_disease_ids AS\n        SELECT DISTINCT canonical_disease_id\n        FROM disease_crosswalk_rows\n        WHERE TRY_CAST(stable_entity_eligible AS BOOLEAN)\n        "
    )
    connection.execute(
        "\n        CREATE VIEW stable_target_ids AS\n        SELECT DISTINCT canonical_target_id\n        FROM target_crosswalk_rows\n        WHERE TRY_CAST(stable_entity_eligible AS BOOLEAN)\n        "
    )
    discrepancy = connection.sql(
        "\n        SELECT COUNT(*) AS discrepancy_count\n        FROM primary_member_rows m\n        LEFT JOIN stable_disease_ids d USING (canonical_disease_id)\n        LEFT JOIN stable_target_ids t USING (canonical_target_id)\n        WHERE CAST(m.stable_entity_member AS BOOLEAN) <> (d.canonical_disease_id IS NOT NULL AND t.canonical_target_id IS NOT NULL)\n        "
    ).fetchone()[0]
    if (int(discrepancy) != 0):
        raise RuntimeError(
            "Staged stable-entity membership does not reproduce the frozen composite flag"
        )
    connection.execute(
        "\n        CREATE VIEW ablation_members AS\n        SELECT 'PRIMARY_FIXED' AS restriction_stage, m.* FROM primary_member_rows m\n        UNION ALL\n        SELECT 'STABLE_DISEASE_ONLY' AS restriction_stage, m.*\n        FROM primary_member_rows m\n        JOIN stable_disease_ids d USING (canonical_disease_id)\n        UNION ALL\n        SELECT 'STABLE_TARGET_ONLY' AS restriction_stage, m.*\n        FROM primary_member_rows m\n        JOIN stable_target_ids t USING (canonical_target_id)\n        UNION ALL\n        SELECT 'STABLE_DISEASE_AND_TARGET' AS restriction_stage, m.*\n        FROM primary_member_rows m\n        JOIN stable_disease_ids d USING (canonical_disease_id)\n        JOIN stable_target_ids t USING (canonical_target_id)\n        "
    )
    panels = connection.sql(
        "\n        WITH eligibility AS (\n            SELECT\n                restriction_stage,\n                release_pair,\n                canonical_disease_id,\n                COUNT(*) AS fixed_target_count,\n                COUNT(DISTINCT old_score) AS old_distinct_score_count,\n                COUNT(DISTINCT new_score) AS new_distinct_score_count\n            FROM ablation_members\n            GROUP BY restriction_stage, release_pair, canonical_disease_id\n            HAVING COUNT(*) >= 30\n               AND COUNT(DISTINCT old_score) >= 2\n               AND COUNT(DISTINCT new_score) >= 2\n        ), positioned AS (\n            SELECT\n                m.*,\n                e.fixed_target_count,\n                e.old_distinct_score_count,\n                e.new_distinct_score_count,\n                MAX(m.old_score) OVER (PARTITION BY m.restriction_stage, m.release_pair, m.canonical_disease_id) AS old_maximum,\n                MAX(m.new_score) OVER (PARTITION BY m.restriction_stage, m.release_pair, m.canonical_disease_id) AS new_maximum,\n                ROW_NUMBER() OVER (PARTITION BY m.restriction_stage, m.release_pair, m.canonical_disease_id ORDER BY m.old_score DESC, m.canonical_target_id) AS old_position,\n                ROW_NUMBER() OVER (PARTITION BY m.restriction_stage, m.release_pair, m.canonical_disease_id ORDER BY m.new_score DESC, m.canonical_target_id) AS new_position\n            FROM ablation_members m\n            JOIN eligibility e USING (restriction_stage, release_pair, canonical_disease_id)\n        ), thresholds AS (\n            SELECT\n                restriction_stage,\n                release_pair,\n                canonical_disease_id,\n                MAX(fixed_target_count) AS fixed_target_count,\n                MAX(old_distinct_score_count) AS old_distinct_score_count,\n                MAX(new_distinct_score_count) AS new_distinct_score_count,\n                MAX(old_maximum) AS old_maximum,\n                MAX(new_maximum) AS new_maximum,\n                MAX(old_score) FILTER (WHERE old_position = 10) AS old_top_10_threshold,\n                MAX(new_score) FILTER (WHERE new_position = 10) AS new_top_10_threshold\n            FROM positioned\n            GROUP BY restriction_stage, release_pair, canonical_disease_id\n        ), aggregated AS (\n            SELECT\n                p.restriction_stage,\n                p.release_pair,\n                p.canonical_disease_id,\n                MAX(t.fixed_target_count) AS fixed_target_count,\n                MAX(t.old_distinct_score_count) AS old_distinct_score_count,\n                MAX(t.new_distinct_score_count) AS new_distinct_score_count,\n                COUNT(*) FILTER (WHERE p.old_score = t.old_maximum) AS old_leader_count,\n                COUNT(*) FILTER (WHERE p.new_score = t.new_maximum) AS new_leader_count,\n                COUNT(*) FILTER (WHERE p.old_score = t.old_maximum AND p.new_score = t.new_maximum) AS shared_leader_count,\n                COUNT(*) FILTER (WHERE p.old_score >= t.old_top_10_threshold) AS old_top_10_size,\n                COUNT(*) FILTER (WHERE p.new_score >= t.new_top_10_threshold) AS new_top_10_size,\n                COUNT(*) FILTER (WHERE p.old_score >= t.old_top_10_threshold AND p.new_score >= t.new_top_10_threshold) AS shared_top_10_count\n            FROM positioned p\n            JOIN thresholds t USING (restriction_stage, release_pair, canonical_disease_id)\n            GROUP BY p.restriction_stage, p.release_pair, p.canonical_disease_id\n        )\n        SELECT\n            *,\n            old_leader_count <> new_leader_count OR shared_leader_count <> old_leader_count AS top_ranked_set_changed,\n            shared_leader_count::DOUBLE / (old_leader_count + new_leader_count - shared_leader_count) AS top_ranked_jaccard,\n            old_top_10_size <> new_top_10_size OR shared_top_10_count <> old_top_10_size AS top_10_changed,\n            shared_top_10_count::DOUBLE / (old_top_10_size + new_top_10_size - shared_top_10_count) AS top_10_jaccard\n        FROM aggregated\n        ORDER BY restriction_stage, release_pair, canonical_disease_id\n        "
    ).df()
    write_frame(
        panels,
        output_dir / "stable_entity_ablation_panels.csv",
        ["restriction_stage", "release_pair", "canonical_disease_id"],
        written,
    )
    primary_counts = (
        panels.loc[panels["restriction_stage"] == "PRIMARY_FIXED"]
        .groupby("release_pair")
        .size()
    )
    summary_rows = []
    for ((stage, release_pair), group) in (panels.groupby(
        ["restriction_stage", "release_pair"], sort = True
    )):
        primary_count = int(primary_counts.loc[release_pair])
        summary_rows.append(
            {
                "restriction_stage": stage,
                "release_pair": release_pair,
                "eligible_panel_count": len(group),
                "eligible_target_pair_count": int(group["fixed_target_count"].sum()),
                "panels_removed_from_primary_count": primary_count - len(group),
                "panels_removed_from_primary_fraction": float(
                    (primary_count - len(group)) / primary_count
                ),
                "median_fixed_target_count": float(
                    group["fixed_target_count"].median()
                ),
                "top_ranked_set_changed_count": int(
                    group["top_ranked_set_changed"].sum()
                ),
                "top_ranked_set_changed_fraction": float(
                    group["top_ranked_set_changed"].mean()
                ),
                "top_10_changed_count": int(group["top_10_changed"].sum()),
                "top_10_changed_fraction": float(group["top_10_changed"].mean()),
                "median_top_10_jaccard": float(group["top_10_jaccard"].median()),
            }
        )
    summary = pd.DataFrame(summary_rows)
    write_frame(
        summary,
        output_dir / "stable_entity_ablation_summary.csv",
        ["restriction_stage", "release_pair"],
        written,
    )
    reported = connection.sql(
        "\n        SELECT\n            release_pair,\n            COUNT(*) AS reported_panel_count,\n            SUM(CASE WHEN CAST(top_ranked_set_changed AS BOOLEAN) THEN 1 ELSE 0 END) AS reported_changed_count,\n            SUM(CASE WHEN CAST(top_10_changed AS BOOLEAN) THEN 1 ELSE 0 END) AS reported_top_10_changed_count\n        FROM stable_metric_rows\n        GROUP BY release_pair\n        "
    ).df()
    reproduced = summary.loc[
        summary["restriction_stage"] == "STABLE_DISEASE_AND_TARGET",
        [
            "release_pair",
            "eligible_panel_count",
            "top_ranked_set_changed_count",
            "top_10_changed_count",
        ],
    ]
    validation = reproduced.merge(reported, on = "release_pair", validate = "one_to_one")
    validation["membership_discrepancy_count"] = int(discrepancy)
    validation["passed"] = (
        validation["eligible_panel_count"].eq(validation["reported_panel_count"])
        & validation["top_ranked_set_changed_count"].eq(
            validation["reported_changed_count"]
        )
        & validation["top_10_changed_count"].eq(
            validation["reported_top_10_changed_count"]
        )
    )
    write_frame(
        validation,
        output_dir / "stable_entity_ablation_validation.csv",
        ["release_pair"],
        written,
    )
    if (not validation["passed"].all()):
        raise RuntimeError(
            "Staged stable-entity restriction does not reproduce frozen stable metrics"
        )

def top_set(targets, scores, cutoff = None):
    targets = np.asarray(targets, dtype = object)
    scores = np.asarray(scores, dtype = float)
    eligible = scores > 0.0
    targets = targets[eligible]
    scores = scores[eligible]
    if (len(scores) == 0):
        return None
    if (cutoff is not None and len(scores) < cutoff):
        return None
    threshold = (
        float(scores.max()) if (cutoff is None) else float(np.sort(scores)[-cutoff])
    )
    return set(targets[scores >= threshold].tolist())

def jaccard(left, right):
    if (left is None or right is None):
        return math.nan
    return len(left & right) / len(left | right)

def calculate_independent_metrics(group, old_field = "old_score", new_field = "new_score"):
    ordered = group.sort_values("canonical_target_id", kind = "mergesort")
    targets = ordered["canonical_target_id"].astype(str).to_numpy()
    old_scores = ordered[old_field].astype(float).to_numpy()
    new_scores = ordered[new_field].astype(float).to_numpy()
    tau = float(
        kendalltau(old_scores, new_scores, variant = "b", method = "auto").statistic
    )
    rho = float(spearmanr(old_scores, new_scores).statistic)
    old_ranks = rankdata(-old_scores, method = "average")
    new_ranks = rankdata(-new_scores, method = "average")
    displacement = np.abs(new_ranks - old_ranks) / (len(targets) - 1)
    result = {
        "fixed_target_count": len(targets),
        "kendall_tau_b": tau,
        "spearman_correlation": rho,
        "median_normalised_rank_displacement": float(
            np.quantile(displacement, 0.5, method = "linear")
        ),
        "p90_normalised_rank_displacement": float(
            np.quantile(displacement, 0.9, method = "linear")
        ),
        "maximum_normalised_rank_displacement": float(displacement.max()),
    }
    for (cutoff) in ([5, 10, 20]):
        old = top_set(targets, old_scores, cutoff)
        new = top_set(targets, new_scores, cutoff)
        result[f"top_{cutoff}_jaccard"] = jaccard(old, new)
        result[f"top_{cutoff}_changed"] = old != new
        result[f"old_top_{cutoff}_effective_size"] = len(old)
        result[f"new_top_{cutoff}_effective_size"] = len(new)
    old = top_set(targets, old_scores)
    new = top_set(targets, new_scores)
    result["old_top_ranked_target_ids"] = ";".join(sorted(old))
    result["new_top_ranked_target_ids"] = ";".join(sorted(new))
    result["top_ranked_set_changed"] = old != new
    result["top_ranked_jaccard"] = jaccard(old, new)
    return result

def build_matched_support(connection, output_dir, written):
    matched = connection.sql(
        "\n        SELECT\n            g.release_pair,\n            g.canonical_disease_id,\n            g.canonical_target_id,\n            p.old_score AS overall_old_score,\n            p.new_score AS overall_new_score,\n            g.old_score AS genetic_old_score,\n            g.new_score AS genetic_new_score\n        FROM genetic_member_rows g\n        JOIN primary_member_rows p USING (release_pair, canonical_disease_id, canonical_target_id)\n        ORDER BY g.release_pair, g.canonical_disease_id, g.canonical_target_id\n        "
    ).df()
    panel_rows = []
    for ((release_pair, disease), group) in (matched.groupby(
        ["release_pair", "canonical_disease_id"], sort = True
    )):
        if (len(group) < 30):
            continue
        endpoint_fields = {
            "OVERALL_SCORE_ON_MATCHED_GENETIC_SUPPORT": (
                "overall_old_score",
                "overall_new_score",
            ),
            "GENETIC_SCORE_ON_MATCHED_SUPPORT": (
                "genetic_old_score",
                "genetic_new_score",
            ),
        }
        if (any(
            (
                group[old].nunique() < 2 or group[new].nunique() < 2
                for ((old, new)) in (endpoint_fields.values())
            )
        )):
            continue
        for (endpoint, (old, new)) in (endpoint_fields.items()):
            metrics = calculate_independent_metrics(group, old, new)
            panel_rows.append(
                {
                    "analysis_layer": endpoint,
                    "release_pair": release_pair,
                    "canonical_disease_id": disease,
                    "support_denominator": "identical matched disease-target support",
                    **metrics,
                }
            )
    panels = pd.DataFrame(panel_rows)
    write_frame(
        panels,
        output_dir / "matched_overall_genetic_panel_metrics.csv",
        ["analysis_layer", "release_pair", "canonical_disease_id"],
        written,
    )
    summary_rows = []
    for ((layer, release_pair), group) in (panels.groupby(
        ["analysis_layer", "release_pair"], sort = True
    )):
        summary_rows.append(
            {
                "analysis_layer": layer,
                "release_pair": release_pair,
                "eligible_panel_count": len(group),
                "eligible_target_pair_count": int(group["fixed_target_count"].sum()),
                "median_fixed_target_count": float(
                    group["fixed_target_count"].median()
                ),
                "median_kendall_tau_b": float(group["kendall_tau_b"].median()),
                "top_ranked_set_changed_count": int(
                    group["top_ranked_set_changed"].sum()
                ),
                "top_ranked_set_changed_fraction": float(
                    group["top_ranked_set_changed"].mean()
                ),
                "top_10_changed_count": int(group["top_10_changed"].sum()),
                "top_10_changed_fraction": float(group["top_10_changed"].mean()),
                "median_top_10_jaccard": float(group["top_10_jaccard"].median()),
                "support_denominator": "identical matched disease-target support",
            }
        )
    summary = pd.DataFrame(summary_rows)
    write_frame(
        summary,
        output_dir / "matched_overall_genetic_summary.csv",
        ["analysis_layer", "release_pair"],
        written,
    )
    pivot = panels.pivot(
        index = ["release_pair", "canonical_disease_id"],
        columns = "analysis_layer",
        values = ["top_ranked_set_changed", "top_10_changed"],
    ).reset_index()
    pivot.columns = [
        "_".join((part for (part) in (column) if (part)))
        if (isinstance(column, tuple))
        else column
        for (column) in (pivot.columns)
    ]
    write_frame(
        pivot,
        output_dir / "matched_overall_genetic_concordance.csv",
        ["release_pair", "canonical_disease_id"],
        written,
    )

def run(project_root, output_dir = None, threads = 4, memory_limit = "8GB"):
    project_root = project_root.resolve()
    derived = project_root / "data" / "derived"
    output_dir = (output_dir or derived / "robustness").resolve()
    output_dir.mkdir(parents = True, exist_ok = True)
    connection, temporary = configure_connection(project_root, threads, memory_limit)
    written = []
    try:
        paths = create_views(connection, derived)
        specifications = layer_specifications()
        exact_frames = build_epsilon_outputs(
            connection, specifications, output_dir, written
        )
        margin_panels = build_margin_outputs(
            connection, specifications, exact_frames, output_dir, written
        )
        build_endpoint_matrix(paths, output_dir, written)
        build_native_distribution(connection, output_dir, written)
        build_genetic_diagnostics(connection, margin_panels, output_dir, written)
        build_ablation_outputs(connection, output_dir, written)
        build_matched_support(connection, output_dir, written)
    finally:
        connection.close()
        shutil.rmtree(temporary, ignore_errors = True)
    return output_dir
