from . import panels as PANELS
import math
from fractions import Fraction
import numpy as np

def bank(samples, measurements, indices):
    expected = {
        (method, version)
        for (method) in (PANELS.METHODS)
        for (version) in ((*PANELS.POLICIES.VERSIONS, "RC"))
    }
    if (set(measurements) != expected):
        raise ValueError("Complete two-method eight-candidate menu required")
    result = {}
    for (candidate, rows) in (measurements.items()):
        if (set(rows) != set(samples)):
            raise ValueError("Candidate population mismatch")
        values = [rows[sample]["hits"] for (sample) in (samples)]
        if (any(
            (
                value is not None and (not isinstance(value, (Fraction, int)))
                for (value) in (values)
            )
        )):
            raise ValueError("Exact accepted hit values required")
        if (any((value is None for (value) in (values)))):
            result[candidate] = (None, None)
        else:
            vector = np.asarray([float(value) for (value) in (values)])
            result[candidate] = (
                sum(values) / len(samples),
                vector[indices].mean(axis = 1),
            )
    return result

def interval(point, draws, total_draws):
    if (draws is None):
        return {
            "value": point,
            "lower_95": None,
            "upper_95": None,
            "draws": total_draws,
            "available_draws": 0,
            "status": "UNAVAILABLE_FULL_POPULATION",
        }
    draws = np.asarray(draws, dtype = float)
    if (draws.shape != (total_draws,)):
        raise ValueError("No resample may be removed")
    available = np.isfinite(draws)
    if (not available.all()):
        return {
            "value": point,
            "lower_95": None,
            "upper_95": None,
            "draws": total_draws,
            "available_draws": int(available.sum()),
            "status": "UNDEFINED_RESAMPLES_RETAINED_NO_CONDITIONAL_INTERVAL",
        }
    (lower, upper) = np.quantile(draws, [0.025, 0.975], method = "linear")
    return {
        "value": point,
        "lower_95": float(lower),
        "upper_95": float(upper),
        "draws": total_draws,
        "available_draws": total_draws,
        "status": "CONDITIONAL_FIXED_CHOICES_PAIRED_SAMPLE_PERCENTILE",
    }

def subtract(first, second):
    if (first[0] is None or second[0] is None or first[1] is None or (second[1] is None)):
        return (None, None)
    return (first[0] - second[0], first[1] - second[1])

def difference_of_gaps(source, target, first, second):
    if (first is None or second is None):
        return (None, None)
    return subtract(
        subtract(source[first], source[second]), subtract(target[first], target[second])
    )

def loss(values, selected, menu):
    if (
        selected is None
        or selected not in menu
        or any(
            (
                values[candidate][0] is None or values[candidate][1] is None
                for (candidate) in (menu)
            )
        )
    ):
        return (None, None)
    point = max((values[candidate][0] for (candidate) in (menu))) - values[selected][0]
    resampled = (
        np.maximum.reduce([values[candidate][1] for (candidate) in (menu)])
        - values[selected][1]
    )
    return (point, resampled)

def choice_menu(choices, policy):
    mapping = choices["policies"][policy]
    if (set(mapping) != set(PANELS.METHODS)):
        raise ValueError("Frozen two-method policy mapping required")
    menu = [
        tuple(mapping[method]) if (mapping[method] is not None) else None
        for (method) in (PANELS.METHODS)
    ]
    ordering = [
        tuple(candidate) for (candidate) in (choices["development_order"][policy])
    ]
    if (len(ordering) != len(set(ordering)) or any(
        (candidate not in menu for (candidate) in (ordering))
    )):
        raise ValueError("Frozen winner is not from this policy menu")
    return (menu, ordering)

def selection_intervals(samples, measurements, choices, bootstrap):
    indices = PANELS.checked_indices(samples, bootstrap)
    values = bank(samples, measurements, indices)
    native = [candidate for (candidate) in (values) if (candidate[1] != "RC")]
    result = {}
    for (policy) in (("V0", "V1", "V2", "MVS", "RC")):
        (menu, order) = choice_menu(choices, policy)
        selected = order[0] if (order) else None
        shared = loss(values, selected, menu) if (None not in menu) else (None, None)
        expanded = loss(values, selected, list(values) if (policy == "RC") else native)
        result[policy] = {
            "fixed_winner": selected,
            "same_two_candidate_policy_loss": interval(*shared, len(indices)),
            "enlarged_eight_loss" if (policy == "RC") else "six_native_loss": interval(
                *expanded, len(indices)
            ),
            "method_gap_PersonaDrive_minus_DawnRank": interval(
                *subtract(values[menu[1]], values[menu[0]]), len(indices)
            )
            if (None not in menu)
            else interval(None, None, len(indices)),
            "choice_reselected_in_draw": False,
        }
    return result

def transfer_intervals(
    samples, string_measurements, e_measurements, choices, bootstrap, support_contract
):
    if (
        support_contract.get("same_ordered_mutation_eligible_sets") is not True
        or support_contract.get(
            "STRING_values_reprojected_to_fixed_resource_intersection"
        )
        is not True
    ):
        raise ValueError(
            "Identical support re-projection is required before numerical transfer comparisons"
        )
    indices = PANELS.checked_indices(samples, bootstrap)
    source = bank(samples, string_measurements, indices)
    target = bank(samples, e_measurements, indices)
    result = {
        "per_candidate_E_minus_STRING": {
            method + "::" + version: interval(
                *subtract(target[method, version], source[method, version]),
                len(indices),
            )
            for ((method, version)) in (source)
        },
        "fixed_policy_gap_transfer": {},
    }
    for (policy) in (("V0", "V1", "V2", "MVS", "RC")):
        (menu, order) = choice_menu(choices, policy)
        result["fixed_policy_gap_transfer"][policy] = interval(
            *difference_of_gaps(source, target, menu[1], menu[0]), len(indices)
        )
    (menu, order) = choice_menu(choices, "MVS")
    (winner, runner) = order[:2] if (len(order) == 2) else (None, None)
    selected_transfer = difference_of_gaps(source, target, winner, runner)
    result["selected_advantage_STRING_minus_BioGRID"] = interval(
        *selected_transfer, len(indices)
    )
    result["excess_selected_advantage_loss_versus_fixed"] = {}
    for (version) in (PANELS.POLICIES.VERSIONS):
        fixed = difference_of_gaps(
            source,
            target,
            (winner[0], version) if (winner) else None,
            (runner[0], version) if (runner) else None,
        )
        result["excess_selected_advantage_loss_versus_fixed"][version] = interval(
            *subtract(selected_transfer, fixed), len(indices)
        )
    return result | {
        "shared_sample_indices_across_resources": True,
        "choices_fixed_from_STRING_development": True,
        "STRING_side_same_cohort_resubstitution_not_new_development_fit": True,
        "independent_biological_validation": False,
        "support_contract": support_contract,
    }

def capture_interval(samples, outcomes, probabilities, bootstrap):
    indices = PANELS.checked_indices(samples, bootstrap)
    if (set(outcomes) != set(samples) or set(probabilities) != set(samples)):
        raise ValueError("Complete ordered capture population required")
    if (any((value not in (None, 0, 1) for (value) in (outcomes.values())))):
        raise ValueError("Binary one-hit target required")
    if (any(
        (
            value is not None
            and (not math.isfinite(float(value)) or not 0 <= value <= 1)
            for (value) in (probabilities.values())
        )
    )):
        raise ValueError("Invalid frozen probability")
    if (any(
        (value is None for (value) in ((*outcomes.values(), *probabilities.values())))
    )):
        return interval(None, None, len(indices)) | {
            "allocation": None,
            "allocation_recomputed_in_draw": False,
        }
    order = sorted(
        samples, key = lambda sample: (-probabilities[sample], sample.encode("utf-8"))
    )
    allocated = set(order[: math.ceil(0.25 * len(samples))])
    vector = np.asarray([outcomes[sample] for (sample) in (samples)])
    captured = np.asarray(
        [outcomes[sample] * (sample in allocated) for (sample) in (samples)]
    )
    denominator = vector[indices].sum(axis = 1)
    numerator = captured[indices].sum(axis = 1)
    draws = np.full(len(indices), np.nan)
    np.divide(numerator, denominator, out = draws, where = denominator != 0)
    point = None if (vector.sum() == 0) else float(captured.sum() / vector.sum())
    return interval(point, draws, len(indices)) | {
        "allocation": [sample for (sample) in (order) if (sample in allocated)],
        "allocation_recomputed_in_draw": False,
        "zero_event_draws": int((denominator == 0).sum()),
        "budget_fraction": 0.25,
    }
