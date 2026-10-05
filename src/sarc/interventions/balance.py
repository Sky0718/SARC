import math
import statistics
from collections import Counter

DEGREE_LIMITS = (
    0,
    1,
    2,
    3,
    4,
    5,
    6,
    8,
    10,
    13,
    16,
    20,
    25,
    32,
    40,
    50,
    64,
    80,
    100,
    125,
    160,
    200,
    250,
    320,
    400,
    500,
    640,
    800,
    1000,
    1250,
    1600,
    2000,
    2500,
    3200,
    4000,
    5000,
    6400,
    8000,
    10000,
    12800,
    16000,
    20000,
    25600,
    32000,
    40000,
    51200,
    64000,
    80000,
    100000,
)
BALANCE_LIMIT = 0.1

def degree_bin(value):
    if (value < 0 or value > DEGREE_LIMITS[-1]):
        raise ValueError("Degree outside fixed matching strata")
    return next(
        (index for ((index, limit)) in (enumerate(DEGREE_LIMITS)) if (value <= limit))
    )

def standardised_difference(real, matched):
    if (len(real) != len(matched)):
        raise ValueError("Balance arrays have unequal multiplicity")
    if (not real):
        return {
            "n": 0,
            "absolute_smd": None,
            "pass": True,
            "state": "EMPTY_MATCHED_SET",
        }
    (first_mean, second_mean) = (statistics.fmean(real), statistics.fmean(matched))
    (first_var, second_var) = (
        statistics.pvariance(real),
        statistics.pvariance(matched),
    )
    scale = math.sqrt((first_var + second_var) / 2)
    if (scale == 0):
        score = 0.0 if (first_mean == second_mean) else None
        passed = first_mean == second_mean
    else:
        score = abs(first_mean - second_mean) / scale
        passed = score <= BALANCE_LIMIT
    return {
        "n": len(real),
        "real_mean_log_degree": first_mean,
        "control_mean_log_degree": second_mean,
        "absolute_smd": score,
        "pass": passed,
        "state": "DEFINED" if (score is not None) else "ZERO_VARIANCE_DIFFERENT_MEANS",
    }

def endpoint_arrays(edges, degrees):
    pairs = [
        sorted((math.log1p(degrees[a]), math.log1p(degrees[b]))) for ((a, b)) in (edges)
    ]
    return {
        "pooled": [value for (pair) in (pairs) for (value) in (pair)],
        "lower": [pair[0] for (pair) in (pairs)],
        "higher": [pair[1] for (pair) in (pairs)],
    }

def verify_matches(
    baseline,
    observed,
    targets_remove,
    targets_add,
    controls_remove,
    controls_add,
    degrees,
):
    from . import matching as controls

    result = controls.verify_matches(
        baseline,
        observed,
        targets_remove,
        targets_add,
        controls_remove,
        controls_add,
        degrees,
    )
    available = Counter(
        (
            controls.removal_stratum(edge, score, degrees)
            for ((edge, score)) in (baseline.items())
        )
    )
    for (label) in ((*controls.CLASSES, "PB_PW")):
        labels = ("PB", "PW") if (label == "PB_PW") else (label,)
        selected_real_remove = set().union(
            *(targets_remove[item] for (item) in (labels))
        )
        selected_real_add = set().union(*(targets_add[item] for (item) in (labels)))
        selected_control_remove = set().union(
            *(controls_remove[item] for (item) in (labels))
        )
        selected_control_add = set().union(
            *(set(controls_add[item]) for (item) in (labels))
        )
        balance = {}
        for (direction, real, matched) in ((
            ("removed", selected_real_remove, selected_control_remove),
            ("added", selected_real_add, selected_control_add),
        )):
            real_arrays = endpoint_arrays(real, degrees)
            matched_arrays = endpoint_arrays(matched, degrees)
            balance[direction] = {
                role: standardised_difference(real_arrays[role], matched_arrays[role])
                for (role) in (real_arrays)
            }
        counts = Counter(
            (
                controls.removal_stratum(edge, baseline[edge], degrees)
                for (edge) in (selected_real_remove)
            )
        )
        forced = sum(
            (
                max(0, 2 * count - available[stratum])
                for ((stratum, count)) in (counts.items())
            )
        )
        entry = result.setdefault(
            label,
            {"removed": len(selected_real_remove), "added": len(selected_real_add)},
        )
        entry.update(
            {
                "continuous_log_degree_balance": balance,
                "conservative_forced_original_removal_lower_bound": forced,
                "forced_fraction": forced / len(selected_real_remove)
                if (selected_real_remove)
                else None,
                "actual_original_removal_overlap_fraction": len(
                    selected_real_remove & selected_control_remove
                )
                / len(selected_real_remove)
                if (selected_real_remove)
                else None,
                "actual_original_addition_overlap_fraction": len(
                    selected_real_add & selected_control_add
                )
                / len(selected_real_add)
                if (selected_real_add)
                else None,
            }
        )
    failures = [
        {"class": label, "direction": direction, "role": role, **entry}
        for ((label, item)) in (result.items())
        for ((direction, values)) in (item["continuous_log_degree_balance"].items())
        for ((role, entry)) in (values.items())
        if (not entry["pass"])
    ]
    if (failures):
        raise ValueError(
            {"message": "Continuous-degree balance failed", "failures": failures}
        )
        raise ValueError("Fixed continuous-degree balance criterion failed")
    return result
