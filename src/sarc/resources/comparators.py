from collections import defaultdict
from copy import deepcopy
from math import isfinite, log1p
from time import perf_counter
import numpy as np
from scipy.stats import kendalltau, spearmanr
from .contracts import ContractError, PanelSnapshot, align_identities, require_count

COMPARATOR_VERSION = "1.0.0"
STRATEGIES = (
    "native_only",
    "fixed_only",
    "fixed_and_native",
    "support_explicit",
    "head_weighted_overlap",
)

def unavailable(reason, **fields):
    return {
        "estimable": False,
        "value": None,
        "support_denominator": 0,
        "reason": reason,
        **fields,
    }

def available(value, **fields):
    return {
        "estimable": True,
        "value": value,
        "support_denominator": 1,
        "reason": None,
        **fields,
    }

def ranked_groups(scores):
    groups = defaultdict(list)
    for (member, score) in (scores.items()):
        groups[score].append(member)
    return [tuple(sorted(groups[value])) for (value) in (sorted(groups, reverse = True))]

def rank_intervals(groups):
    result = {}
    rank = 1
    for (group) in (groups):
        bottom = rank + len(group) - 1
        for (member) in (group):
            result[member] = (rank, bottom)
        rank = bottom + 1
    return result

def add_polynomial(differences, start, stop, coefficients, limit):
    start = max(1, int(start))
    stop = min(limit, int(stop))
    if (start > stop):
        return
    for (order, coefficient) in (enumerate(coefficients)):
        differences[order, start] += coefficient
        differences[order, stop + 1] -= coefficient

def polynomial_values(differences, limit):
    depths = np.arange(1, limit + 1, dtype = float)
    values = np.cumsum(differences, axis = 1)[:, 1 : limit + 1]
    return values[0] + values[1] * depths + values[2] * depths * depths

def contribution_sums(intervals, limit, variant, squared = False, starts = None):
    differences = np.zeros((3, limit + 2), dtype = float)
    for (member, (top, bottom)) in (intervals.items()):
        start = top if (starts is None) else max(top, starts[member])
        if (variant == "w"):
            add_polynomial(differences, start, limit, (1, 0, 0), limit)
        else:
            width = bottom - top + 1
            slope = 1 / width
            intercept = (1 - top) / width
            coefficients = (
                (intercept * intercept, 2 * intercept * slope, slope * slope)
                if (squared)
                else (intercept, slope, 0)
            )
            add_polynomial(differences, start, bottom - 1, coefficients, limit)
            add_polynomial(differences, max(start, bottom), limit, (1, 0, 0), limit)
    return polynomial_values(differences, limit)

def overlap_sums(short, long, limit, variant):
    differences = np.zeros((3, limit + 2), dtype = float)
    for (member) in (short.keys() & long.keys()):
        (first_top, first_bottom) = short[member]
        (second_top, second_bottom) = long[member]
        start = max(first_top, second_top)
        if (variant == "w"):
            add_polynomial(differences, start, limit, (1, 0, 0), limit)
            continue
        boundaries = sorted(
            {start, max(start, first_bottom), max(start, second_bottom), limit + 1}
        )
        for (lower, upper) in (zip(boundaries[:-1], boundaries[1:])):
            first_slope = (
                0 if (lower >= first_bottom) else 1 / (first_bottom - first_top + 1)
            )
            second_slope = (
                0 if (lower >= second_bottom) else 1 / (second_bottom - second_top + 1)
            )
            first_intercept = 1 if (first_slope == 0) else (1 - first_top) * first_slope
            second_intercept = (
                1 if (second_slope == 0) else (1 - second_top) * second_slope
            )
            coefficients = (
                first_intercept * second_intercept,
                first_slope * second_intercept + second_slope * first_intercept,
                first_slope * second_slope,
            )
            add_polynomial(differences, lower, upper - 1, coefficients, limit)
    return polynomial_values(differences, limit)

def rbo_agreements(first_groups, second_groups, variant):
    (first, second) = (rank_intervals(first_groups), rank_intervals(second_groups))
    (short, long) = (first, second) if (len(first) <= len(second)) else (second, first)
    (shorter, longer) = (len(short), len(long))
    if (shorter == 0):
        raise ContractError("RBO requires two nonempty observed prefixes")
    depths = np.arange(1, longer + 1, dtype = float)
    overlap = overlap_sums(short, long, longer, variant)
    short_norm = contribution_sums(short, longer, variant, squared = variant == "b")
    long_norm = contribution_sums(long, longer, variant, squared = variant == "b")
    if (variant == "w"):
        common_agreement = (
            2 * overlap[:shorter] / (short_norm[:shorter] + long_norm[:shorter])
        )
        tail_denominator = (depths[shorter:] + long_norm[shorter:]) / 2
    else:
        common_agreement = overlap[:shorter] / np.sqrt(
            short_norm[:shorter] * long_norm[:shorter]
        )
        tail_denominator = np.sqrt(depths[shorter:] * long_norm[shorter:])
    minimum = np.empty(longer, dtype = float)
    maximum = np.empty(longer, dtype = float)
    extrapolated = np.empty(longer, dtype = float)
    minimum[:shorter] = common_agreement
    maximum[:shorter] = common_agreement
    extrapolated[:shorter] = common_agreement
    unique_long = {
        member: positions
        for ((member, positions)) in (long.items())
        if (member not in short)
    }
    if (longer > shorter):
        starts = {
            member: shorter + offset + 1
            for ((offset, member)) in (enumerate(unique_long))
        }
        best_unseen = contribution_sums(unique_long, longer, variant, starts = starts)
        active_unique = contribution_sums(unique_long, longer, "w")
        unique_contribution = contribution_sums(unique_long, longer, variant)
        unseen_mean = unique_contribution[shorter:] / active_unique[shorter:]
        extrapolated_unseen = (
            (depths[shorter:] - shorter) * common_agreement[-1] * unseen_mean
        )
        minimum[shorter:] = overlap[shorter:] / tail_denominator
        maximum[shorter:] = (
            overlap[shorter:] + best_unseen[shorter:]
        ) / tail_denominator
        extrapolated[shorter:] = (
            overlap[shorter:] + extrapolated_unseen
        ) / tail_denominator
    return (
        shorter,
        longer,
        len(short.keys() & long.keys()),
        common_agreement[-1],
        minimum,
        maximum,
        extrapolated,
    )

def rbo_from_agreements(agreement, persistence):
    (shorter, longer, shared, last_agreement, minimum, maximum, extrapolated) = (
        agreement
    )
    if (not isfinite(persistence) or not 0 < persistence < 1):
        raise ContractError("RBO persistence must lie strictly between zero and one")
    depths = np.arange(1, longer + 1, dtype = float)
    weights = np.power(persistence, depths)
    prefix_log = -log1p(-persistence) - float(np.sum(weights / depths))
    if (prefix_log < 1e-10):
        depth = longer + 1
        power = persistence**depth
        prefix_log = 0.0
        while (power / depth > 1e-18):
            prefix_log += power / depth
            depth += 1
            power *= persistence
    tail_minimum = shared * prefix_log
    join_depth = longer + shorter - shared
    future_depths = np.arange(longer + 1, join_depth + 1, dtype = float)
    tail_maximum = float(
        np.sum(
            (2 * future_depths - longer - shorter + shared)
            / future_depths
            * np.power(persistence, future_depths)
        )
    ) + persistence ** (join_depth + 1) / (1 - persistence)
    tail_extrapolated = (
        (shared + (longer - shorter) * last_agreement)
        / longer
        * persistence ** (longer + 1)
        / (1 - persistence)
    )
    factor = (1 - persistence) / persistence
    result = {
        "min": factor * (float(np.dot(minimum, weights)) + tail_minimum),
        "max": factor * (float(np.dot(maximum, weights)) + tail_maximum),
        "ext": factor * (float(np.dot(extrapolated, weights)) + tail_extrapolated),
    }
    for (key, value) in (result.items()):
        if (not isfinite(value) or value < -1e-10 or value > 1 + 1e-10):
            raise ContractError("RBO escaped its numerical domain")
        result[key] = min(1.0, max(0.0, value))
    result["res"] = result["max"] - result["min"]
    if (result["res"] < -1e-10):
        raise ContractError("RBO bounds are reversed")
    return result

def tie_aware_rbo(first, second, *, persistence = (0.9, 0.8, 0.95), variants = ("w", "b")):
    if (any((variant not in {"w", "b"} for (variant) in (variants)))):
        raise ContractError("only frozen RBO variants w and b are permitted")
    if (not first or not second):
        return unavailable("EMPTY_SUPPORT")
    groups = (ranked_groups(first), ranked_groups(second))
    results = {}
    for (variant) in (variants):
        agreement = rbo_agreements(*groups, variant)
        results[variant] = {
            str(value): rbo_from_agreements(agreement, value)
            for (value) in (persistence)
        }
    return available(
        results,
        prefix_lengths = [len(first), len(second)],
        interpretation = "observed_prefix_bounds_and_extrapolation_not_sampling_intervals",
    )

def set_comparison(left, right, *, top_k = 1, minimum_roster = 1):
    count = min(len(left), len(right))
    if (count < max(top_k, minimum_roster)):
        return unavailable(
            "EMPTY_SUPPORT" if (count == 0) else "INSUFFICIENT_MEMBERS",
            baseline_count = len(left),
            followup_count = len(right),
            baseline_set = None,
            followup_set = None,
            jaccard = None,
        )
    left_cutoff = sorted(left.values(), reverse = True)[top_k - 1]
    right_cutoff = sorted(right.values(), reverse = True)[top_k - 1]
    first = frozenset(
        (member for ((member, value)) in (left.items()) if (value >= left_cutoff))
    )
    second = frozenset(
        (member for ((member, value)) in (right.items()) if (value >= right_cutoff))
    )
    (numerator, denominator) = (len(first & second), len(first | second))
    jaccard = {
        "value": numerator / denominator,
        "numerator": numerator,
        "denominator": denominator,
        "support_denominator": denominator,
    }
    return available(
        first != second,
        baseline_count = len(left),
        followup_count = len(right),
        baseline_set = sorted(first),
        followup_set = sorted(second),
        jaccard = jaccard,
    )

def rank_comparison(left, right, minimum_roster):
    count = len(left)
    if (count < minimum_roster):
        return unavailable(
            "EMPTY_SUPPORT" if (count == 0) else "INSUFFICIENT_MEMBERS",
            member_count = count,
            kendall_tau_b = None,
            spearman = None,
        )
    if (len(set(left.values())) < 2 or len(set(right.values())) < 2):
        return unavailable(
            "CONSTANT_SCORES", member_count = count, kendall_tau_b = None, spearman = None
        )
    members = sorted(left)
    (first, second) = (
        [left[member] for (member) in (members)],
        [right[member] for (member) in (members)],
    )
    tau = float(kendalltau(first, second, variant = "b").statistic)
    spearman = float(spearmanr(first, second).statistic)
    if (not isfinite(tau) or not isfinite(spearman)):
        raise ContractError("eligible rank comparison produced a nonfinite value")
    return available(tau, member_count = count, kendall_tau_b = tau, spearman = spearman)

def score_comparison(left, right):
    if (not left):
        return unavailable("EMPTY_SUPPORT", member_count = 0)
    difference = np.array([right[member] - left[member] for (member) in (sorted(left))])
    absolute = np.abs(difference)
    return available(
        float(np.median(absolute)),
        member_count = len(left),
        pair_support_denominator = len(left),
        mean_signed_change = float(np.mean(difference)),
        mean_absolute_change = float(np.mean(absolute)),
        median_absolute_change = float(np.quantile(absolute, 0.5, method = "linear")),
        q90_absolute_change = float(np.quantile(absolute, 0.9, method = "linear")),
        maximum_absolute_change = float(np.max(absolute)),
        changed_pair_count = int(np.count_nonzero(difference)),
    )

def task_answer(endpoint, *, reported = True):
    if (not reported):
        return {
            "state": "NOT_REPORTED",
            "answer": None,
            "support_denominator": None,
            "reason": "OUTSIDE_DECLARED_REPORT",
        }
    return {
        "state": "ANSWERED" if (endpoint["estimable"]) else "ABSTAINED",
        "answer": endpoint["value"],
        "support_denominator": endpoint["support_denominator"],
        "reason": endpoint["reason"],
    }

def strategy_reports(report):
    result = {}
    membership = report["membership"]
    ambiguity = available(not report["identity"]["decidable"])
    for (strategy) in (STRATEGIES):
        native = strategy in {"native_only", "fixed_and_native", "support_explicit"}
        fixed = strategy in {"fixed_only", "fixed_and_native", "support_explicit"}
        result[strategy] = {
            "answers": {
                "fixed_leading_set_change": task_answer(
                    report["fixed"]["leader"], reported = fixed
                ),
                "native_leading_set_change": task_answer(
                    report["native"]["leader"], reported = native
                ),
                "persistent_versus_roster_membership_description": task_answer(
                    membership, reported = strategy != "head_weighted_overlap"
                ),
                "mapping_ambiguity_requires_abstention": task_answer(ambiguity),
            },
            "information_policy": "same_records_and_mapping; membership_derivation_allowed; continuous_RBO_not_binary_warning",
        }
    return result

def compare_panels(
    baseline: PanelSnapshot,
    followup: PanelSnapshot,
    *,
    mapping_proposals = (),
    minimum_roster = 30,
    top_ks = (5, 10, 20),
    rbo_p = (0.9, 0.8, 0.95),
    rbo_variants = ("w", "b"),
    include_rbo = True,
):
    started = perf_counter()
    if (not isinstance(baseline, PanelSnapshot) or not isinstance(
        followup, PanelSnapshot
    )):
        raise ContractError("both inputs must be validated PanelSnapshot records")
    if (baseline.resource != followup.resource or baseline.panel_id != followup.panel_id):
        raise ContractError(
            "paired snapshots require the same resource and canonical panel identifier"
        )
    if (require_count(minimum_roster, "minimum_roster") == 0):
        raise ContractError("minimum_roster must be positive")
    for (top_k) in (top_ks):
        if (require_count(top_k, "top_k") == 0):
            raise ContractError("top_k must be positive")
    (left, raw_right) = (dict(baseline.members), dict(followup.members))
    alignment = align_identities(left, raw_right, mapping_proposals)
    decidable = (
        baseline.identity_decidable
        and followup.identity_decidable
        and (not alignment.ambiguous_baseline)
        and (not alignment.ambiguous_followup)
    )
    right_to_left = {later: earlier for ((earlier, later, _)) in (alignment.matches)}
    right = {
        right_to_left.get(member, member): value
        for ((member, value)) in (raw_right.items())
    }
    (persistent, entering, exiting) = (
        set(left) & set(right),
        set(right) - set(left),
        set(left) - set(right),
    )
    report = {
        "comparator_version": COMPARATOR_VERSION,
        "resource": baseline.resource,
        "releases": [baseline.release, followup.release],
        "panel_id": baseline.panel_id,
        "minimum_roster": minimum_roster,
        "identity": {
            "decidable": bool(decidable),
            "matches": list(alignment.matches),
            "ambiguous_baseline": sorted(alignment.ambiguous_baseline),
            "ambiguous_followup": sorted(alignment.ambiguous_followup),
        },
        "support": {
            "baseline": len(left),
            "followup": len(right),
            "fixed": len(persistent) if (decidable) else None,
            "entering": len(entering) if (decidable) else None,
            "exiting": len(exiting) if (decidable) else None,
        },
    }
    if (not decidable):
        endpoint = unavailable("IDENTITY_UNRESOLVED")
        report["native"] = {
            "leader": dict(endpoint),
            "top_k": {str(value): dict(endpoint) for (value) in (top_ks)},
            "rbo": dict(endpoint),
        }
        report["fixed"] = {
            "leader": dict(endpoint),
            "top_k": {str(value): dict(endpoint) for (value) in (top_ks)},
            "rank": dict(endpoint),
            "score": dict(endpoint),
            "rbo": dict(endpoint),
        }
        report["membership"] = dict(endpoint)
    else:
        (fixed_left, fixed_right) = (
            {member: left[member] for (member) in (persistent)},
            {member: right[member] for (member) in (persistent)},
        )
        report["native"] = {
            "leader": set_comparison(left, right),
            "top_k": {
                str(value): set_comparison(left, right, top_k = value)
                for (value) in (top_ks)
            },
        }
        report["fixed"] = {
            "leader": set_comparison(
                fixed_left, fixed_right, minimum_roster = minimum_roster
            ),
            "top_k": {
                str(value): set_comparison(
                    fixed_left, fixed_right, top_k = value, minimum_roster = minimum_roster
                )
                for (value) in (top_ks)
            },
            "rank": rank_comparison(fixed_left, fixed_right, minimum_roster),
            "score": score_comparison(fixed_left, fixed_right),
        }
        report["membership"] = available(
            {
                "persistent": sorted(persistent),
                "entering": sorted(entering),
                "exiting": sorted(exiting),
            },
            persistent_count = len(persistent),
            entering_count = len(entering),
            exiting_count = len(exiting),
        )
        report["native"]["rbo"] = (
            tie_aware_rbo(left, right, persistence = rbo_p, variants = rbo_variants)
            if (include_rbo)
            else unavailable("NOT_REQUESTED")
        )
        report["fixed"]["rbo"] = (
            tie_aware_rbo(
                fixed_left, fixed_right, persistence = rbo_p, variants = rbo_variants
            )
            if (include_rbo and len(persistent) >= minimum_roster)
            else unavailable(
                "NOT_REQUESTED" if (not include_rbo) else "INSUFFICIENT_MEMBERS"
            )
        )
    report["strategies"] = strategy_reports(report)
    report["wall_seconds"] = perf_counter() - started
    return report

def apply_minimum_roster(report, minimum_roster):
    if (require_count(minimum_roster, "minimum_roster") < report["minimum_roster"]):
        raise ContractError(
            "a stricter saved report cannot restore previously unavailable estimates"
        )
    result = deepcopy(report)
    result["minimum_roster"] = minimum_roster
    count = result["support"]["fixed"]
    if (count is not None and count < minimum_roster):
        fixed = result["fixed"]
        for (key) in (("leader", "rank", "rbo")):
            fixed[key] = unavailable(
                "EMPTY_SUPPORT" if (count == 0) else "INSUFFICIENT_MEMBERS",
                member_count = count,
            )
        for (top_k) in (fixed["top_k"]):
            fixed["top_k"][top_k] = unavailable(
                "EMPTY_SUPPORT" if (count == 0) else "INSUFFICIENT_MEMBERS",
                member_count = count,
            )
    result["strategies"] = strategy_reports(result)
    return result
