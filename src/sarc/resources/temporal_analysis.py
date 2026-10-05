from collections import Counter
import numpy as np
import pyarrow.parquet as pq
from .contracts import Measure, ResultContext

PRIMARY_MINIMUM = 30
MINIMUMS = (20, 30, 50)
TOP_KS = (5, 10, 20)

def rows_by_panel(path):
    current = None
    rows = []
    for (batch) in (pq.ParquetFile(path).iter_batches(batch_size = 32768, use_threads = False)):
        for (row) in (batch.to_pylist()):
            if (current is not None and row["disease_id"] != current):
                yield (current, rows)
                rows = []
            current = row["disease_id"]
            rows.append(row)
    if (current is not None):
        yield (current, rows)

def flatten_report(report, baseline, minimum):
    result = {
        "disease_id": baseline["disease_id"],
        "resource": "OpenTargets",
        "baseline_release": "26.06",
        "followup_release": "26.09",
        "minimum_roster": minimum,
        "baseline_native_roster_size": baseline["native_roster_size"],
        "mapping_state": baseline["mapping_state"],
        "baseline_native_leader_margin": baseline["native_leader_margin"],
        "baseline_exact_leader_tie_count": baseline["exact_leader_tie_count"],
        "baseline_mapping_status": baseline["baseline_mapping_status"],
        "identity_decidable": report["identity"]["decidable"],
        "comparator_version": report["comparator_version"],
    }
    result.update({f"{key}_n": value for ((key, value)) in (report["support"].items())})
    sets = []
    for (view) in (("native", "fixed")):
        for (endpoint, value) in ([("leader", report[view]["leader"])] + [
            (f"top_{k}", report[view]["top_k"][str(k)]) for (k) in (TOP_KS)
        ]):
            key = f"{view}_{endpoint}"
            result[f"{key}_changed"] = value["value"]
            result[f"{key}_reason"] = value["reason"]
            jaccard = value.get("jaccard")
            result[f"{key}_jaccard"] = (
                jaccard["value"] if (jaccard is not None) else None
            )
            result[f"{key}_jaccard_numerator"] = (
                jaccard["numerator"] if (jaccard is not None) else None
            )
            result[f"{key}_jaccard_denominator"] = (
                jaccard["denominator"] if (jaccard is not None) else None
            )
            sets.append(
                {
                    "disease_id": baseline["disease_id"],
                    "minimum_roster": minimum,
                    "support": view,
                    "endpoint": endpoint,
                    "estimable": value["estimable"],
                    "reason": value["reason"],
                    "baseline_set": value.get("baseline_set"),
                    "followup_set": value.get("followup_set"),
                    "jaccard_numerator": result[f"{key}_jaccard_numerator"],
                    "jaccard_denominator": result[f"{key}_jaccard_denominator"],
                }
            )
    rank = report["fixed"]["rank"]
    result["fixed_kendall_tau_b"] = rank.get("kendall_tau_b")
    result["fixed_spearman"] = rank.get("spearman")
    result["fixed_rank_reason"] = rank["reason"]
    score = report["fixed"]["score"]
    result["score_pair_support"] = score.get("pair_support_denominator", 0)
    result["fixed_score_reason"] = score["reason"]
    for (field) in ((
        "mean_signed_change",
        "mean_absolute_change",
        "median_absolute_change",
        "q90_absolute_change",
        "maximum_absolute_change",
        "changed_pair_count",
    )):
        result[f"fixed_{field}"] = score.get(field)
    result["native_fixed_leader_disagreement"] = (
        None
        if (
            result["native_leader_changed"] is None
            or result["fixed_leader_changed"] is None
        )
        else result["native_leader_changed"] != result["fixed_leader_changed"]
    )
    result["native_fixed_top_10_disagreement"] = (
        None
        if (
            result["native_top_10_changed"] is None
            or result["fixed_top_10_changed"] is None
        )
        else result["native_top_10_changed"] != result["fixed_top_10_changed"]
    )
    result["native_fixed_leader_disagreement_reason"] = (
        result["native_leader_reason"] or result["fixed_leader_reason"]
    )
    result["native_fixed_top_10_disagreement_reason"] = (
        result["native_top_10_reason"] or result["fixed_top_10_reason"]
    )
    return (result, sets)

def resampling_interval(values, *, seed = 960625, replicates = 2000):
    array = np.asarray(
        [np.nan if (value is None) else float(value) for (value) in (values)],
        dtype = float,
    )
    finite = np.isfinite(array)
    if (not finite.any()):
        return (None, None)
    generator = np.random.default_rng(seed)
    estimates = []
    for (offset) in (range(0, replicates, 40)):
        sample = array[
            generator.integers(
                0, len(array), size = (min(40, replicates - offset), len(array))
            )
        ]
        denominator = np.count_nonzero(np.isfinite(sample), axis = 1)
        numerator = np.nansum(sample, axis = 1)
        estimates.extend(
            (numerator[denominator > 0] / denominator[denominator > 0]).tolist()
        )
    if (not estimates):
        return (None, None)
    (lower, upper) = np.quantile(estimates, [0.025, 0.975], method = "linear")
    return (float(lower), float(upper))

def register_row(
    metric,
    values,
    cohort,
    *,
    minimum = 30,
    stratum = "all",
    kind = "fraction",
    interval = False,
    reasons = None,
):
    observed = [value for (value) in (values) if (value is not None)]
    support = len(observed)
    numerator = int(sum(observed)) if (kind == "fraction" and support) else None
    denominator = support if (kind == "fraction" and support) else None
    value = (
        numerator / denominator
        if (kind == "fraction" and support)
        else float(np.mean(observed))
        if (support)
        else None
    )
    reason = None if (support) else "NO_EVALUABLE_UNITS"
    context = ResultContext(
        "OpenTargets",
        ("26.06", "26.09"),
        f"baseline_native_min_{minimum}",
        stratum,
        metric,
        "disease_panel",
        "association_support.parquet;cohort_followup_identity.parquet",
        "OpenTargets 26.06 to 26.09",
    )
    Measure(context, kind, value, support, numerator, denominator, reason)
    (low, high) = (
        resampling_interval(values) if (interval and support) else (None, None)
    )
    return {
        "resource": "OpenTargets",
        "release_pair": "26.06_to_26.09",
        "cohort": f"baseline_native_min_{minimum}",
        "minimum_roster": minimum,
        "stratum": stratum,
        "estimand": metric,
        "measure_type": kind,
        "unit": "disease_panel",
        "value": value,
        "numerator": numerator,
        "denominator": denominator,
        "support_denominator": support,
        "baseline_cohort_size": cohort,
        "stratum_cohort_size": len(values),
        "outside_stratum_count": cohort - len(values),
        "unevaluable_count": len(values) - support,
        "reason": reason,
        "unevaluable_reasons": dict(
            Counter((item for (item) in (reasons or []) if (item is not None)))
        ),
        "descriptive_bootstrap_lower": low,
        "descriptive_bootstrap_upper": high,
        "bootstrap_replicates": 2000 if (interval and support) else 0,
        "bootstrap_seed": 960625 if (interval and support) else None,
        "bootstrap_unit": "baseline_disease_panel_within_declared_stratum"
        if (interval and support)
        else None,
        "interval_interpretation": "descriptive_resampling_not_independent_population_inference"
        if (interval and support)
        else None,
        "source": "panel_endpoints.parquet",
    }

def summarise(rows):
    registry = []
    set_columns = [
        f"{view}_{item}_changed"
        for (view) in (("native", "fixed"))
        for (item) in (("leader", "top_5", "top_10", "top_20"))
    ]
    for (minimum) in (MINIMUMS):
        selected = [row for (row) in (rows) if (row["minimum_roster"] == minimum)]
        for (column) in (set_columns + [
            "native_fixed_leader_disagreement",
            "native_fixed_top_10_disagreement",
        ]):
            reasons = [
                row.get(column.removesuffix("_changed") + "_reason")
                if (row[column] is None)
                else None
                for (row) in (selected)
            ]
            registry.append(
                register_row(
                    column + "_fraction",
                    [row[column] for (row) in (selected)],
                    len(selected),
                    minimum = minimum,
                    interval = True,
                    reasons = reasons,
                )
            )
            registry.append(
                register_row(
                    column + "_coverage",
                    [row[column] is not None for (row) in (selected)],
                    len(selected),
                    minimum = minimum,
                    interval = False,
                )
            )
        for (column, reason) in ((
            ("fixed_kendall_tau_b", "fixed_rank_reason"),
            ("fixed_spearman", "fixed_rank_reason"),
            ("fixed_mean_absolute_change", "fixed_score_reason"),
            ("fixed_median_absolute_change", "fixed_score_reason"),
        )):
            registry.append(
                register_row(
                    column + "_panel_mean",
                    [row[column] for (row) in (selected)],
                    len(selected),
                    minimum = minimum,
                    kind = "scalar",
                    interval = True,
                    reasons = [row[reason] for (row) in (selected)],
                )
            )
        for (view) in (("native", "fixed")):
            for (item) in (("leader", "top_5", "top_10", "top_20")):
                registry.append(
                    register_row(
                        f"{view}_{item}_jaccard_mean",
                        [row[f"{view}_{item}_jaccard"] for (row) in (selected)],
                        len(selected),
                        minimum = minimum,
                        kind = "scalar",
                        interval = False,
                        reasons = [row[f"{view}_{item}_reason"] for (row) in (selected)],
                    )
                )
    return registry

def minimum_comparisons(rows):
    results = []
    for (first_minimum, second_minimum) in (((20, 30), (30, 50), (20, 50))):
        first = {
            row["disease_id"]: row
            for (row) in (rows)
            if (row["minimum_roster"] == first_minimum)
        }
        second = {
            row["disease_id"]: row
            for (row) in (rows)
            if (row["minimum_roster"] == second_minimum)
        }
        for (endpoint) in ((
            "native_leader_changed",
            "native_top_10_changed",
            "fixed_leader_changed",
            "fixed_top_10_changed",
            "fixed_kendall_tau_b",
        )):
            common = sorted(
                (
                    identifier
                    for (identifier) in (first.keys() & second.keys())
                    if (
                        first[identifier][endpoint] is not None
                        and second[identifier][endpoint] is not None
                    )
                )
            )
            differences = [
                float(second[identifier][endpoint]) - float(first[identifier][endpoint])
                for (identifier) in (common)
            ]
            results.append(
                {
                    "first_minimum": first_minimum,
                    "second_minimum": second_minimum,
                    "endpoint": endpoint,
                    "first_baseline_cohort_size": len(first),
                    "second_baseline_cohort_size": len(second),
                    "first_evaluable_count": sum(
                        (row[endpoint] is not None for (row) in (first.values()))
                    ),
                    "second_evaluable_count": sum(
                        (row[endpoint] is not None for (row) in (second.values()))
                    ),
                    "identical_evaluable_panel_count": len(common),
                    "paired_mean_difference": float(np.mean(differences))
                    if (differences)
                    else None,
                    "maximum_absolute_difference": max(
                        (abs(value) for (value) in (differences))
                    )
                    if (differences)
                    else None,
                    "reason": None if (differences) else "NO_COMMON_EVALUABLE_PANELS",
                    "source": "panel_endpoints.parquet",
                }
            )
    return results

def margin_bin(value):
    bounds = (0, 1e-06, 0.0001, 0.001, 0.01, 0.05, 1)
    for (index, upper) in (enumerate(bounds)):
        if (value <= upper):
            return "0" if (index == 0) else f"({bounds[index - 1]:g},{upper:g}]"
    raise ValueError("Baseline leader margin outside frozen bins")

def size_bin(value):
    return "30-99" if (value < 100) else "100-499" if (value < 500) else "500+"

def warning_summaries(rows):
    selected = [row for (row) in (rows) if (row["minimum_roster"] == 30)]
    strata = {
        "margin:" + label: []
        for (label) in (
            (
                "0",
                "(0,1e-06]",
                "(1e-06,0.0001]",
                "(0.0001,0.001]",
                "(0.001,0.01]",
                "(0.01,0.05]",
                "(0.05,1]",
            )
        )
    }
    strata.update({"size:" + label: [] for (label) in (("30-99", "100-499", "500+"))})
    strata.update({"ties:" + label: [] for (label) in (("unique", "tied"))})
    strata.update(
        {"warning:" + label: [] for (label) in (("margin_le_0.01", "margin_gt_0.01"))}
    )
    for (row) in (selected):
        margin = margin_bin(row["baseline_native_leader_margin"])
        size = size_bin(row["baseline_native_roster_size"])
        tie = "tied" if (row["baseline_exact_leader_tie_count"] > 1) else "unique"
        warning = (
            "margin_le_0.01"
            if (row["baseline_native_leader_margin"] <= 0.01)
            else "margin_gt_0.01"
        )
        for (key) in ((
            "margin:" + margin,
            "size:" + size,
            "ties:" + tie,
            "warning:" + warning,
        )):
            strata[key].append(row)
    registry = []
    for (stratum, panels) in (strata.items()):
        registry.append(
            register_row(
                "native_leader_changed_fraction",
                [row["native_leader_changed"] for (row) in (panels)],
                len(selected),
                stratum = stratum,
                interval = True,
                reasons = [row["native_leader_reason"] for (row) in (panels)],
            )
        )
        registry.append(
            register_row(
                "native_leader_coverage",
                [row["native_leader_changed"] is not None for (row) in (panels)],
                len(selected),
                stratum = stratum,
            )
        )
    eligible = [
        row for (row) in (selected) if (row["native_leader_changed"] is not None)
    ]
    confusion = Counter()
    for (row) in (eligible):
        warning = row["baseline_native_leader_margin"] <= 0.01
        event = row["native_leader_changed"]
        confusion[
            ("warned" if (warning) else "unwarned")
            + ("_changed" if (event) else "_unchanged")
        ] += 1
    for (key) in ((
        "warned_changed",
        "warned_unchanged",
        "unwarned_changed",
        "unwarned_unchanged",
    )):
        confusion.setdefault(key, 0)
    event_rows = [row for (row) in (selected) if (row["native_leader_changed"] is True)]
    stable_rows = [
        row for (row) in (selected) if (row["native_leader_changed"] is False)
    ]
    warned = [
        row for (row) in (selected) if (row["baseline_native_leader_margin"] <= 0.01)
    ]
    registry.append(
        register_row(
            "warning_fraction_baseline_cohort",
            [row["baseline_native_leader_margin"] <= 0.01 for (row) in (selected)],
            len(selected),
        )
    )
    registry.append(
        register_row(
            "warning_event_sensitivity",
            [row["baseline_native_leader_margin"] <= 0.01 for (row) in (event_rows)],
            len(selected),
            stratum = "evaluable_changed_panels",
            interval = True,
        )
    )
    registry.append(
        register_row(
            "warning_false_alert_fraction",
            [row["baseline_native_leader_margin"] <= 0.01 for (row) in (stable_rows)],
            len(selected),
            stratum = "evaluable_unchanged_panels",
            interval = True,
        )
    )
    registry.append(
        register_row(
            "warning_missed_event_fraction",
            [row["baseline_native_leader_margin"] > 0.01 for (row) in (event_rows)],
            len(selected),
            stratum = "evaluable_changed_panels",
            interval = True,
        )
    )
    registry.append(
        register_row(
            "warning_positive_predictive_fraction",
            [row["native_leader_changed"] for (row) in (warned)],
            len(selected),
            stratum = "baseline_warned_panels",
            interval = True,
            reasons = [row["native_leader_reason"] for (row) in (warned)],
        )
    )
    return (
        registry,
        {
            "rule": "baseline_native_leader_margin <= 0.01",
            "baseline_feature_release": "26.06",
            "baseline_total": len(selected),
            "evaluable_total": len(eligible),
            "unevaluable_total": len(selected) - len(eligible),
            "counts": dict(confusion),
            "interpretation": "descriptive warning rule, not a fitted or calibrated predictive model",
        },
    )
