from .io import quoted, check, write_json
import itertools
import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
from .contracts import Measure, ResultContext

SET_ENDPOINTS = ("leader", "top5", "top10", "top20")
SCORE_FIELDS = (
    "mean_signed_change",
    "mean_absolute_change",
    "median_absolute_change",
    "q90_absolute_change",
    "maximum_absolute_change",
)
TASKS = (
    "fixed_leading_set_change",
    "native_leading_set_change",
    "persistent_versus_roster_membership_description",
    "mapping_ambiguity_requires_abstention",
)

def grouped_panels(path):
    current = None
    rows = []
    for (batch) in (pq.ParquetFile(path).iter_batches(batch_size = 65536)):
        for (row) in (batch.to_pylist()):
            anchor = row["anchor_id"]
            if (current is not None and anchor != current):
                check(
                    anchor > current,
                    "anchor input must be strictly grouped in ascending order",
                )
                yield (current, rows)
                rows = []
            current = anchor
            rows.append(row)
    if (current is not None):
        yield (current, rows)

def set_endpoint(report, view, name):
    return (
        report[view]["leader"]
        if (name == "leader")
        else report[view]["top_k"][name[3:]]
    )

def panel_schema():
    fields = [
        ("anchor_id", pa.string()),
        ("minimum_roster", pa.int32()),
        ("threshold", pa.int32()),
        ("baseline_degree_150", pa.int32()),
        ("degree_stratum", pa.string()),
        ("identity_state", pa.string()),
    ]
    fields.extend(
        (
            (name + "_count", pa.int32())
            for (name) in (("baseline", "followup", "fixed", "entering", "exiting"))
        )
    )
    for (view, endpoint) in (itertools.product(("native", "fixed"), SET_ENDPOINTS)):
        prefix = view + "_" + endpoint
        fields.extend(
            (
                (prefix + "_changed", pa.bool_()),
                (prefix + "_reason", pa.string()),
                (prefix + "_jaccard", pa.float64()),
                (prefix + "_intersection", pa.int32()),
                (prefix + "_union", pa.int32()),
            )
        )
    fields.extend(
        (
            ("fixed_rank_reason", pa.string()),
            ("kendall_tau_b", pa.float64()),
            ("spearman", pa.float64()),
            ("fixed_score_reason", pa.string()),
            ("score_pair_count", pa.int32()),
            ("changed_score_pair_count", pa.int32()),
        )
    )
    fields.extend(((name, pa.float64()) for (name) in (SCORE_FIELDS)))
    for (view) in (("native", "fixed")):
        fields.append((view + "_rbo_reason", pa.string()))
        fields.extend(
            (
                (
                    f"{view}_rbo_{variant}_p{str(persistence).replace('.', '_')}_{bound}",
                    pa.float64(),
                )
                for ((variant, persistence, bound)) in (
                    itertools.product(
                        ("w", "b"), (0.9, 0.8, 0.95), ("min", "max", "res", "ext")
                    )
                )
            )
        )
    return pa.schema(fields)

def compact_panel(report, threshold, baseline):
    row = {
        "anchor_id": report["panel_id"],
        "minimum_roster": report["minimum_roster"],
        "threshold": threshold,
        "baseline_degree_150": baseline["baseline_degree_150"],
        "degree_stratum": baseline["baseline_degree_stratum"],
        "identity_state": "EXACT_ID"
        if (report["identity"]["decidable"])
        else "UNRESOLVED",
    }
    row.update(
        {name + "_count": value for ((name, value)) in (report["support"].items())}
    )
    for (view, name) in (itertools.product(("native", "fixed"), SET_ENDPOINTS)):
        endpoint = set_endpoint(report, view, name)
        prefix = view + "_" + name
        jaccard = endpoint.get("jaccard") or {}
        row.update(
            {
                prefix + "_changed": endpoint["value"],
                prefix + "_reason": endpoint["reason"],
                prefix + "_jaccard": jaccard.get("value"),
                prefix + "_intersection": jaccard.get("numerator"),
                prefix + "_union": jaccard.get("denominator"),
            }
        )
    (rank, score) = (report["fixed"]["rank"], report["fixed"]["score"])
    row.update(
        {
            "fixed_rank_reason": rank["reason"],
            "kendall_tau_b": rank.get("kendall_tau_b"),
            "spearman": rank.get("spearman"),
            "fixed_score_reason": score["reason"],
            "score_pair_count": score.get("pair_support_denominator"),
            "changed_score_pair_count": score.get("changed_pair_count"),
        }
    )
    row.update({name: score.get(name) for (name) in (SCORE_FIELDS)})
    for (view) in (("native", "fixed")):
        rbo = report[view]["rbo"]
        row[view + "_rbo_reason"] = rbo["reason"]
        for (variant, persistence, bound) in (itertools.product(
            ("w", "b"), (0.9, 0.8, 0.95), ("min", "max", "res", "ext")
        )):
            key = f"{view}_rbo_{variant}_p{str(persistence).replace('.', '_')}_{bound}"
            row[key] = (
                rbo["value"][variant][str(persistence)][bound]
                if (rbo["estimable"])
                else None
            )
    return row

def leading_sets(report, threshold):
    for (view, name) in (itertools.product(("native", "fixed"), SET_ENDPOINTS)):
        endpoint = set_endpoint(report, view, name)
        yield {
            "anchor_id": report["panel_id"],
            "threshold": threshold,
            "minimum_roster": report["minimum_roster"],
            "view": view,
            "endpoint": name,
            "baseline_set": endpoint.get("baseline_set"),
            "followup_set": endpoint.get("followup_set"),
            "reason": endpoint["reason"],
        }

def fraction(metric, numerator, denominator, context, **fields):
    (numerator, denominator) = (int(numerator), int(denominator))
    value = numerator / denominator if (denominator) else None
    reason = None if (denominator) else "NO_ELIGIBLE_SUPPORT"
    checked = Measure(
        ResultContext(
            "STRING",
            ("11.5", "12.0"),
            context,
            context,
            metric,
            "anchor_panel",
            "anchor_neighbour_transitions.parquet",
            "STRING 11.5 to 12.0",
        ),
        "fraction",
        value,
        denominator,
        numerator if (denominator) else None,
        denominator if (denominator) else None,
        reason,
    )
    return {
        "metric": metric,
        "kind": checked.kind,
        "value": checked.value,
        "numerator": checked.numerator,
        "denominator": checked.denominator,
        "support_denominator": checked.support_denominator,
        "reason": checked.reason,
        **fields,
    }

def summarise(output, settings):
    connection = duckdb.connect(config = {"threads": "2", "memory_limit": "1GB"})
    connection.execute(
        f"CREATE VIEW panels AS SELECT * FROM read_parquet({quoted((output / 'panel_endpoints.parquet').as_posix())})"
    )
    connection.execute(
        f"CREATE VIEW strategies AS SELECT * FROM read_parquet({quoted((output / 'strategy_reports.parquet').as_posix())})"
    )
    columns = [field.name for (field) in (panel_schema())]
    groupings = connection.execute(
        "SELECT DISTINCT minimum_roster,threshold,degree_stratum FROM panels ORDER BY ALL"
    ).fetchall()
    groupings += [
        (minimum, threshold, "ALL")
        for ((minimum, threshold)) in (
            itertools.product(settings["minimum_rosters"], settings["thresholds"])
        )
    ]
    summaries = []
    scalar_columns = [
        name
        for (name) in (columns)
        if (
            name in SCORE_FIELDS
            or name in ("kendall_tau_b", "spearman")
            or name.endswith("_jaccard")
            or ("_rbo_" in name and (not name.endswith("_reason")))
        )
    ]
    for (minimum, threshold, stratum) in (groupings):
        condition = f"minimum_roster = {minimum} AND threshold = {threshold}" + (
            "" if (stratum == "ALL") else f" AND degree_stratum = {quoted(stratum)}"
        )
        cohort = connection.execute(
            f"SELECT count(*) FROM panels WHERE {condition}"
        ).fetchone()[0]
        context = {
            "minimum_roster": minimum,
            "threshold": threshold,
            "degree_stratum": stratum,
            "cohort_count": cohort,
        }
        for (view, endpoint) in (itertools.product(("native", "fixed"), SET_ENDPOINTS)):
            column = view + "_" + endpoint + "_changed"
            (eligible, changed) = connection.execute(
                f"SELECT count({column}), count(*) FILTER(WHERE {column}) FROM panels WHERE {condition}"
            ).fetchone()
            summaries.append(
                fraction(
                    view + "_" + endpoint + "_change_fraction",
                    changed,
                    eligible,
                    "endpoint_evaluable_baseline_anchors",
                    **context,
                )
            )
            summaries.append(
                fraction(
                    view + "_" + endpoint + "_coverage",
                    eligible,
                    cohort,
                    "all_baseline_anchors",
                    **context,
                )
            )
        for (endpoint) in (SET_ENDPOINTS):
            (native, fixed) = (
                "native_" + endpoint + "_changed",
                "fixed_" + endpoint + "_changed",
            )
            (eligible, different) = connection.execute(
                f"SELECT count(*) FILTER(WHERE {native} IS NOT NULL AND {fixed} IS NOT NULL), count(*) FILTER(WHERE {native} != {fixed}) FROM panels WHERE {condition}"
            ).fetchone()
            summaries.append(
                fraction(
                    "native_fixed_" + endpoint + "_disagreement",
                    different,
                    eligible,
                    "common_evaluable_anchors",
                    **context,
                )
            )
            summaries.append(
                fraction(
                    "native_fixed_" + endpoint + "_common_coverage",
                    eligible,
                    cohort,
                    "all_baseline_anchors",
                    **context,
                )
            )
        for (column) in (scalar_columns):
            (count, mean, quantiles) = connection.execute(
                f"SELECT count({column}), avg({column}), quantile_cont({column}, [0.5,0.9,0.95,0.99]) FROM panels WHERE {condition}"
            ).fetchone()
            summaries.append(
                {
                    "metric": column,
                    "kind": "scalar_summary",
                    "mean": mean,
                    "quantile_probabilities": settings["quantiles"],
                    "quantiles": quantiles,
                    "support_denominator": count,
                    "numerator": None,
                    "denominator": None,
                    "reason": None if (count) else "NO_ELIGIBLE_SUPPORT",
                    **context,
                }
            )
    strategy_rows = connection.execute(
        "SELECT minimum_roster,threshold,strategy,task,count(*) AS cohort,count(*) FILTER(WHERE eligible) AS eligible_count,count(*) FILTER(WHERE eligible AND state = 'ANSWERED') AS answered_eligible,count(*) FILTER(WHERE NOT eligible) AS ineligible,count(*) FILTER(WHERE NOT eligible AND state = 'ABSTAINED') AS abstained_ineligible,count(*) FILTER(WHERE NOT eligible AND state = 'ANSWERED') AS attempted_ineligible,count(*) FILTER(WHERE state = 'ANSWERED' AND agrees_with_full_report) AS full_agreement,count(*) FILTER(WHERE state = 'ANSWERED') AS answered FROM strategies GROUP BY ALL ORDER BY ALL"
    ).fetchall()
    strategy_summary = []
    for ((
        minimum,
        threshold,
        strategy,
        task,
        cohort,
        eligible,
        answered_eligible,
        ineligible,
        abstained,
        attempted,
        agreement,
        answered,
    )) in (strategy_rows):
        fields = {
            "minimum_roster": minimum,
            "threshold": threshold,
            "strategy": strategy,
            "task": task,
            "cohort_count": cohort,
        }
        for (name, numerator, denominator) in ((
            ("evaluable_answer_coverage", answered_eligible, eligible),
            ("ineligible_answer_fraction", attempted, ineligible),
            ("explicit_abstention_fraction", abstained, ineligible),
            ("agreement_with_full_report_among_answers", agreement, answered),
        )):
            strategy_summary.append(
                fraction(
                    name, numerator, denominator, "declared_task_support", **fields
                )
            )
    contrasts = []
    for (minimum, (first, second), view, endpoint) in (itertools.product(
        settings["minimum_rosters"],
        settings["threshold_contrasts"],
        ("native", "fixed"),
        SET_ENDPOINTS,
    )):
        column = view + "_" + endpoint + "_changed"
        (cohort, common, first_changed, second_changed, discordant) = (
            connection.execute(
                f"SELECT count(*),count(*) FILTER(WHERE a.{column} IS NOT NULL AND b.{column} IS NOT NULL),count(*) FILTER(WHERE a.{column} AND b.{column} IS NOT NULL),count(*) FILTER(WHERE b.{column} AND a.{column} IS NOT NULL),count(*) FILTER(WHERE a.{column} != b.{column}) FROM panels a JOIN panels b USING(anchor_id,minimum_roster) WHERE a.minimum_roster = {minimum} AND a.threshold = {first} AND b.threshold = {second}"
            ).fetchone()
        )
        fields = {
            "minimum_roster": minimum,
            "thresholds": [first, second],
            "view": view,
            "endpoint": endpoint,
            "cohort_count": cohort,
        }
        contrasts.append(
            {
                **fields,
                "common_evaluable_count": common,
                "first_changed_count": first_changed,
                "second_changed_count": second_changed,
                "difference_second_minus_first": (second_changed - first_changed)
                / common
                if (common)
                else None,
                "common_coverage": fraction(
                    "common_evaluable_coverage",
                    common,
                    cohort,
                    "identical_baseline_anchors",
                ),
                "discordance": fraction(
                    "threshold_change_indicator_discordance",
                    discordant,
                    common,
                    "identical_common_evaluable_anchors",
                ),
            }
        )
    missingness = []
    for (view, endpoint) in (itertools.product(("native", "fixed"), SET_ENDPOINTS)):
        column = view + "_" + endpoint + "_reason"
        for (minimum, threshold, reason, count) in (connection.execute(
            f"SELECT minimum_roster,threshold,{column},count(*) FROM panels GROUP BY ALL ORDER BY ALL"
        ).fetchall()):
            missingness.append(
                {
                    "minimum_roster": minimum,
                    "threshold": threshold,
                    "view": view,
                    "endpoint": endpoint,
                    "state": "ESTIMABLE" if (reason is None) else reason,
                    "count": count,
                }
            )
    write_json(output / "endpoint_summaries.json", summaries)
    write_json(output / "strategy_summaries.json", strategy_summary)
    write_json(output / "threshold_contrasts.json", contrasts)
    write_json(output / "endpoint_coverage_flow.json", missingness)
    connection.close()
