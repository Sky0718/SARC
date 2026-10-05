from .io import check
import itertools
from .string_transfer import SET_ENDPOINTS, fraction

def contrast_records(connection, minimums, thresholds):
    records = []
    endpoints = [
        view + "_" + name + "_changed"
        for ((view, name)) in (itertools.product(("native", "fixed"), SET_ENDPOINTS))
    ] + ["kendall_tau_b", "spearman", "mean_absolute_change"]
    for ((first, second), threshold, endpoint) in (itertools.product(
        itertools.combinations(minimums, 2), thresholds, endpoints
    )):
        (first_cohort, first_eligible) = connection.execute(
            f"SELECT count(*),count({endpoint}) FROM panels WHERE minimum_roster = {first} AND threshold = {threshold}"
        ).fetchone()
        (second_cohort, second_eligible) = connection.execute(
            f"SELECT count(*),count({endpoint}) FROM panels WHERE minimum_roster = {second} AND threshold = {threshold}"
        ).fetchone()
        (shared, common, difference, maximum) = connection.execute(
            f"SELECT count(*),count(*) FILTER(WHERE a.{endpoint} IS NOT NULL AND b.{endpoint} IS NOT NULL),avg(b.{endpoint}::DOUBLE - a.{endpoint}::DOUBLE) FILTER(WHERE a.{endpoint} IS NOT NULL AND b.{endpoint} IS NOT NULL),max(abs(b.{endpoint}::DOUBLE - a.{endpoint}::DOUBLE)) FILTER(WHERE a.{endpoint} IS NOT NULL AND b.{endpoint} IS NOT NULL) FROM panels a JOIN panels b USING(anchor_id,threshold) WHERE a.threshold = {threshold} AND a.minimum_roster = {first} AND b.minimum_roster = {second}"
        ).fetchone()
        check(
            shared == second_cohort,
            "stricter baseline minimum is not nested in the lower cohort",
        )
        records.append(
            {
                "resource": "STRING",
                "releases": ["11.5", "12.0"],
                "threshold": threshold,
                "first_minimum": first,
                "second_minimum": second,
                "endpoint": endpoint,
                "first_baseline_cohort_count": first_cohort,
                "second_baseline_cohort_count": second_cohort,
                "shared_baseline_cohort_count": shared,
                "first_evaluable_count": first_eligible,
                "second_evaluable_count": second_eligible,
                "common_evaluable_count": common,
                "paired_mean_difference_second_minus_first": difference,
                "maximum_absolute_paired_difference": maximum,
                "reason": None if (common) else "NO_COMMON_EVALUABLE_PANELS",
                "common_evaluable_coverage": fraction(
                    "common_evaluable_coverage",
                    common,
                    shared,
                    "identical_baseline_anchors",
                ),
                "interpretation": "Minimum changes eligibility, not scores or ranking on the same jointly evaluable anchor; cohort-level marginal rates are not used as paired effects",
            }
        )
    return records
