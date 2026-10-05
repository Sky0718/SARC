from . import scoring as CORE
from . import policies as POLICIES

RELEASES = ("4_1_190", "4_4_200", "4_4_223")
CONTROL_SEEDS = (271828, 314159, 161803, 141421, 173205)
import itertools
import math
from collections import Counter
from fractions import Fraction
import numpy as np

METHODS = ("DawnRank", "PersonaDrive")
REFERENCES = (
    "NCG6_primary_all",
    "NCG6_known_subgroup",
    "NCG6_candidate_subgroup",
    "CancerMine2019_secondary",
)
SUPPORTS = ("XE_common", "native", "fixed_resource_intersection")
KS = (1, 5, 10, 20)
TRANSITIONS = ("4_1_190_to_4_4_200", "4_4_200_to_4_4_223")
DELTAS = (Fraction(1, 2), Fraction(1), Fraction(2))
FAMILIES = (
    "strict_input_only",
    "edit_volume_only",
    "degree_incidence_only",
    "additive_volume_degree",
    "old_output_margin",
    "input_plus_margin",
)

def delta(left, right):
    return None if (left is None or right is None) else left - right

def scalar(vector, indices):
    n = len(vector)
    available = [value for (value) in (vector) if (value is not None)]
    result = {
        "population_count": n,
        "complete_count": len(available),
        "value": None,
        "complete_subset_value_secondary": sum(available) / len(available)
        if (available)
        else None,
        "lower_95": None,
        "upper_95": None,
    }
    if (len(available) != n):
        return result | {"status": "UNAVAILABLE_FULL_POPULATION_NO_RESAMPLE_DROPPED"}
    if (n == 0):
        raise ValueError("Empty planned population")
    values = np.asarray([float(value) for (value) in (vector)])
    draws = values[indices].mean(axis = 1)
    (lower, upper) = np.quantile(draws, [0.025, 0.975], method = "linear")
    return result | {
        "value": sum(vector) / n,
        "lower_95": float(lower),
        "upper_95": float(upper),
        "status": "CONDITIONAL_PAIRED_SAMPLE_DESCRIPTIVE_PERCENTILE",
        "draws": len(indices),
    }

def checked_indices(samples, record):
    if (
        record["sample_order"] != samples
        or record["seed"] != 20260926
        or record["draws"] != 2000
    ):
        raise ValueError("Inherited ordered paired draw contract changed")
    indices = np.asarray(record["indices"])
    if (
        indices.shape != (2000, len(samples))
        or indices.dtype.kind not in "iu"
        or np.any(indices < 0)
        or np.any(indices >= len(samples))
    ):
        raise ValueError("Invalid complete paired draw geometry")
    return indices.astype(np.int64, copy = False)

def score_population(samples, ranks, eligibility, positives, k):
    if (set(ranks) != set(samples) or set(eligibility) != set(samples)):
        raise ValueError("Full paired model population required")
    values = {}
    for (sample) in (samples):
        rank = ranks[sample]
        values[sample] = CORE.score_list(
            rank["genes"], rank["scores"], eligibility[sample], positives, k
        )
    return values

def summary(samples, values, indices):
    return CORE.mean_records([values[sample] for (sample) in (samples)]) | {
        "interval": scalar([values[sample]["hits"] for (sample) in (samples)], indices)
    }

def contrast(samples, changed, baseline, indices):
    result = CORE.paired_contrast(samples, [(1, changed), (-1, baseline)])
    differences = [row["hits"] for (row) in (result["sample_contrasts"])]
    return result | {
        "interval": scalar(differences, indices),
        "absolute_change": scalar(
            [None if (value is None) else abs(value) for (value) in (differences)],
            indices,
        ),
        "one_hit_event_fraction": scalar(
            [
                None if (value is None) else int(abs(value) >= 1)
                for (value) in (differences)
            ],
            indices,
        ),
    }

def membership_population(samples, left, right, eligibility, k, indices):
    rows = {}
    for (sample) in (samples):
        memberships = [
            CORE.leading_membership(
                rank[sample]["genes"], rank[sample]["scores"], eligibility[sample], k
            )
            for (rank) in ((left, right))
        ]
        values = CORE.leading_stability(*memberships)
        rows[sample] = values | {
            "old_candidate_count": None
            if (left[sample]["genes"] is None)
            else len(set(left[sample]["genes"]) & set(eligibility[sample])),
            "new_candidate_count": None
            if (right[sample]["genes"] is None)
            else len(set(right[sample]["genes"]) & set(eligibility[sample])),
            "eligible_count": len(eligibility[sample]),
        }
    return {
        "samples": rows,
        "summaries": {
            key: scalar([rows[sample][key] for (sample) in (samples)], indices)
            for (key) in (
                (
                    "jaccard",
                    "fractional_membership_jaccard",
                    "old_candidate_count",
                    "new_candidate_count",
                    "eligible_count",
                )
            )
        },
        "successful_both_empty_count": sum(
            (row["status"] == "BOTH_SUCCESSFUL_EMPTY" for (row) in (rows.values()))
        ),
    }

def candidate_ranks(samples, method, getter, eligibility):
    native = {
        version: getter(method, "XE__native_" + release)
        for ((version, release)) in (zip(POLICIES.VERSIONS, RELEASES))
    }
    rc = {}
    for (sample) in (samples):
        consensus = POLICIES.consensus(
            [
                (native[version][sample]["genes"], native[version][sample]["scores"])
                for (version) in (POLICIES.VERSIONS)
            ],
            eligibility[sample],
        )
        rc[sample] = {
            "genes": consensus["genes"],
            "scores": consensus["scores"],
            "status": consensus["status"],
        }
    return native | {"RC": rc}

def selection_panel(samples, candidate_scores, choices, indices):
    menu = {
        (method, version)
        for (method) in (METHODS)
        for (version) in ((*POLICIES.VERSIONS, "RC"))
    }
    if (set(candidate_scores) != menu or set(choices["policies"]) != {
        "V0",
        "V1",
        "V2",
        "MVS",
        "RC",
    }):
        raise ValueError("Complete frozen two-method candidate menu required")
    means = {
        candidate: CORE.mean_records([rows[sample] for (sample) in (samples)])["value"]
        for ((candidate, rows)) in (candidate_scores.items())
    }
    native = {
        candidate: value
        for ((candidate, value)) in (means.items())
        if (candidate[1] != "RC")
    }
    result = {}
    for (policy) in (("V0", "V1", "V2", "MVS", "RC")):
        mapping = choices["policies"][policy]
        if (set(mapping) != set(METHODS)):
            raise ValueError("A three-method choice cannot silently become two-method")
        selected = {
            method: tuple(candidate) if (candidate is not None) else None
            for ((method, candidate)) in (mapping.items())
        }
        if (any(
            (
                candidate is not None
                and (candidate not in menu or candidate[0] != method)
                for ((method, candidate)) in (selected.items())
            )
        )):
            raise ValueError("Unfrozen candidate redirection")
        ordered = [
            tuple(candidate) for (candidate) in (choices["development_order"][policy])
        ]
        if (len(set(ordered)) != len(ordered) or any(
            (candidate not in selected.values() for (candidate) in (ordered))
        )):
            raise ValueError("Frozen development choice does not match policy menu")
        winner = ordered[0] if (ordered) else None
        policy_menu = {
            candidate: means[candidate]
            for (candidate) in (selected.values())
            if (candidate is not None)
        }
        gap = None
        if (all(selected.values())):
            gap = CORE.paired_contrast(
                samples,
                [
                    (1, candidate_scores[selected["PersonaDrive"]]),
                    (-1, candidate_scores[selected["DawnRank"]]),
                ],
            )
            gap["interval"] = scalar(
                [row["hits"] for (row) in (gap["sample_contrasts"])], indices
            )
        result[policy] = {
            "method_candidates": selected,
            "development_order": ordered,
            "winner": winner,
            "method_means": {
                method: means.get(candidate)
                for ((method, candidate)) in (selected.items())
            },
            "PersonaDrive_minus_DawnRank": gap,
            "loss_same_two_candidate_policy": POLICIES.portfolio_loss(
                policy_menu, winner
            )
            if (len(policy_menu) == 2)
            else None,
            "loss_six_native": POLICIES.portfolio_loss(native, winner)
            if (policy != "RC")
            else None,
            "loss_enlarged_eight": POLICIES.portfolio_loss(means, winner)
            if (policy == "RC")
            else None,
            "E_outcomes_used_for_selection": False,
        }
    return result

def controls_panel(samples, baseline, real, controls, indices):
    if (tuple(controls) != CONTROL_SEEDS):
        raise ValueError("All five original ordered graph-control seeds required")
    effects = {
        seed: contrast(samples, rows, baseline, indices)
        for ((seed, rows)) in (controls.items())
    }
    real_effect = contrast(samples, real, baseline, indices)
    complete = all(
        (
            all((row["hits"] is not None for (row) in (rows.values())))
            for (rows) in ((baseline, real, *controls.values()))
        )
    )
    terms = [(1, real), (-1, baseline)] + [
        (coefficient, values)
        for (seed) in (CONTROL_SEEDS)
        for ((coefficient, values)) in (
            ((Fraction(-1, 5), controls[seed]), (Fraction(1, 5), baseline))
        )
    ]
    matched = CORE.paired_contrast(samples, terms)
    if (not complete and matched["value"] is not None):
        raise ValueError("Incomplete matched family was incorrectly collapsed")
    matched["interval"] = scalar(
        [row["hits"] for (row) in (matched["sample_contrasts"])], indices
    )
    absolute = []
    for (i, sample) in (enumerate(samples)):
        actual = real_effect["sample_contrasts"][i]["hits"]
        alternatives = [
            effects[seed]["sample_contrasts"][i]["hits"] for (seed) in (CONTROL_SEEDS)
        ]
        absolute.append(
            None
            if (actual is None or any((value is None for (value) in (alternatives))))
            else abs(actual) - sum((abs(value) for (value) in (alternatives))) / 5
        )
    return {
        "status": "COMPLETE_FIVE_CONTROL_FAMILY"
        if (complete)
        else "UNAVAILABLE_FULL_FIVE_CONTROL_FAMILY",
        "real": real_effect,
        "five_controls": effects,
        "signed_real_minus_mean_five_controls": matched,
        "absolute_real_minus_mean_five_controls": scalar(absolute, indices),
        "unavailable_controls_not_dropped": True,
        "control_seeds_are_not_people": True,
        "randomisation_p_value": None,
    }

def paired_resource(samples, e_values, string_values, indices, same_support_binding):
    if (
        same_support_binding.get("same_ordered_mutation_eligible_sets") is not True
        or same_support_binding.get(
            "STRING_values_reprojected_to_fixed_resource_intersection"
        )
        is not True
    ):
        raise ValueError(
            "Cross-resource magnitude requires actual identical support re-projection"
        )
    return contrast(samples, e_values, string_values, indices) | {
        "support_binding": same_support_binding,
        "independent_biological_cohort": False,
    }

def forecast_panel(rows, probabilities, samples, indices):
    if (set(rows) != set(samples) or set(probabilities) != set(samples)):
        raise ValueError("Prediction and outcome population mismatch")
    values = [rows[sample] for (sample) in (samples)]
    risk = [probabilities[sample] for (sample) in (samples)]
    if (any((value not in (None, 0, 1) for (value) in (values)))):
        raise ValueError("Frozen one-hit target must be binary or unavailable")
    if (any(
        (
            value is not None
            and (not math.isfinite(float(value)) or not 0 <= value <= 1)
            for (value) in (risk)
        )
    )):
        raise ValueError("Malformed frozen probability")
    complete = all((value is not None for (value) in (values + risk)))
    order = (
        sorted(
            samples,
            key = lambda sample: (-float(probabilities[sample]), sample.encode("utf-8")),
        )
        if (all((value is not None for (value) in (risk))))
        else None
    )
    allocated = None if (order is None) else order[: math.ceil(0.25 * len(samples))]
    brier = [
        None if (target is None or p is None) else (float(p) - target) ** 2
        for ((target, p)) in (zip(values, risk))
    ]
    capture = None
    total = None if (not complete) else sum(values)
    if (complete and total):
        capture = sum((rows[sample] for (sample) in (allocated))) / total
    calibration = []
    for (i) in (range(10)):
        selected = [
            sample
            for (sample) in (samples)
            if (
                probabilities[sample] is not None
                and probabilities[sample] >= i / 10
                and (
                    probabilities[sample] < (i + 1) / 10
                    if (i < 9)
                    else probabilities[sample] <= 1
                )
            )
        ]
        assessed = [sample for (sample) in (selected) if (rows[sample] is not None)]
        calibration.append(
            {
                "lower": i / 10,
                "upper": (i + 1) / 10,
                "selected_count": len(selected),
                "complete_count": len(assessed),
                "mean_probability": None
                if (not selected)
                else sum((probabilities[sample] for (sample) in (selected)))
                / len(selected),
                "event_fraction": None
                if (not selected or len(assessed) != len(selected))
                else sum((rows[sample] for (sample) in (selected))) / len(selected),
            }
        )
    return {
        "status": "COMPLETE" if (complete) else "UNAVAILABLE_FULL_POPULATION",
        "brier": scalar(brier, indices),
        "calibration": calibration,
        "budget_fraction": 0.25,
        "allocation": allocated,
        "capture": capture,
        "total_events": total,
        "E_refit": False,
        "tie_rule": "descending frozen risk then UTF-8 model identity",
    }

def evaluate_cohort(
    cohort, samples, getter, support, references, choices, bootstrap, stop = lambda: None
):
    indices = checked_indices(samples, bootstrap)
    if (set(references) != set(REFERENCES) or set(support) != set(samples)):
        raise ValueError("Complete inherited references and support required")
    panels = []
    cached_ranks = {}

    def ranks(method, network):
        key = (method, network)
        if (key not in cached_ranks):
            cached_ranks[key] = getter(method, network)
        return cached_ranks[key]

    for (support_name) in (SUPPORTS):
        eligibility = {
            sample: set(support[sample][support_name]) for (sample) in (samples)
        }
        candidates = {
            method: candidate_ranks(samples, method, ranks, eligibility)
            for (method) in (METHODS)
        }
        for (reference, k) in (itertools.product(REFERENCES, KS)):
            stop()
            positives = set(references[reference])
            scores = {
                (method, version): score_population(
                    samples, candidates[method][version], eligibility, positives, k
                )
                for (method) in (METHODS)
                for (version) in ((*POLICIES.VERSIONS, "RC"))
            }
            identity = {
                "cohort": cohort,
                "support": support_name,
                "reference": reference,
                "k": k,
                "primary": support_name == "XE_common"
                and reference == "NCG6_primary_all"
                and (k == 10),
                "cross_resource_projection": support_name
                == "fixed_resource_intersection",
            }
            panels.append(
                identity
                | {
                    "panel": "XA_TWO_METHOD_FIXED_POLICY",
                    "policies": selection_panel(samples, scores, choices, indices),
                    "candidate_summaries": {
                        method + "::" + version: summary(samples, values, indices)
                        for (((method, version), values)) in (scores.items())
                    },
                }
            )
            for (transition_index, transition) in (enumerate(TRANSITIONS)):
                for (family) in (("native", "persistent")):
                    endpoints = {}
                    for (method) in (METHODS):
                        if (family == "native"):
                            (old, new) = (
                                candidates[method][POLICIES.VERSIONS[position]]
                                for (position) in (
                                    (transition_index, transition_index + 1)
                                )
                            )
                        else:
                            (old, new) = (
                                ranks(
                                    method, "XE__" + transition + "__persistent_" + side
                                )
                                for (side) in (("old", "new"))
                            )
                        pair = [
                            score_population(samples, value, eligibility, positives, k)
                            for (value) in ((old, new))
                        ]
                        endpoints[method] = pair
                        panels.append(
                            identity
                            | {
                                "panel": family.upper() + "_RELEASE_CHANGE",
                                "transition": transition,
                                "method": method,
                                "effect": contrast(samples, pair[1], pair[0], indices),
                                "membership_and_coverage": membership_population(
                                    samples, old, new, eligibility, k, indices
                                ),
                            }
                        )
                    gaps = [
                        CORE.paired_contrast(
                            samples,
                            [
                                (1, endpoints["PersonaDrive"][side]),
                                (-1, endpoints["DawnRank"][side]),
                            ],
                        )
                        for (side) in ((0, 1))
                    ]
                    intervals = [
                        scalar(
                            [row["hits"] for (row) in (gap["sample_contrasts"])],
                            indices,
                        )
                        for (gap) in (gaps)
                    ]
                    panels.append(
                        identity
                        | {
                            "panel": family.upper() + "_METHOD_GAP_AND_REVERSAL",
                            "transition": transition,
                            "old_gap": gaps[0],
                            "new_gap": gaps[1],
                            "old_interval": intervals[0],
                            "new_interval": intervals[1],
                            "gap_change": CORE.paired_contrast(
                                samples,
                                [
                                    (1, endpoints["PersonaDrive"][1]),
                                    (-1, endpoints["DawnRank"][1]),
                                    (-1, endpoints["PersonaDrive"][0]),
                                    (1, endpoints["DawnRank"][0]),
                                ],
                            ),
                            "reversals": [
                                CORE.reversals(
                                    gaps[0]["value"],
                                    gaps[1]["value"],
                                    intervals[0],
                                    intervals[1],
                                    threshold,
                                )
                                for (threshold) in (DELTAS)
                            ],
                        }
                    )
            for (method, direction) in (itertools.product(METHODS, ("insert", "retract"))):
                prefix = "XE__" + TRANSITIONS[-1]
                baseline_rank = ranks(
                    method,
                    prefix
                    + "__persistent_"
                    + ("old" if (direction == "insert") else "new"),
                )
                actual_rank = ranks(
                    method, prefix + "__" + cohort + "__" + direction + "__PB"
                )
                control_ranks = {
                    seed: ranks(
                        method,
                        prefix
                        + "__"
                        + cohort
                        + "__"
                        + direction
                        + "__control_"
                        + str(seed)
                        + "__PB",
                    )
                    for (seed) in (CONTROL_SEEDS)
                }
                baseline = score_population(
                    samples, baseline_rank, eligibility, positives, k
                )
                actual = score_population(
                    samples, actual_rank, eligibility, positives, k
                )
                controls = {
                    seed: score_population(samples, rank, eligibility, positives, k)
                    for ((seed, rank)) in (control_ranks.items())
                }
                panels.append(
                    identity
                    | {
                        "panel": "PB_MATCHED_LOCALISATION",
                        "method": method,
                        "direction": direction,
                        "transition": TRANSITIONS[-1],
                        "effect": controls_panel(
                            samples, baseline, actual, controls, indices
                        ),
                        "real_membership_and_coverage": membership_population(
                            samples, baseline_rank, actual_rank, eligibility, k, indices
                        ),
                        "five_control_membership_and_coverage": {
                            seed: membership_population(
                                samples, baseline_rank, rank, eligibility, k, indices
                            )
                            for ((seed, rank)) in (control_ranks.items())
                        },
                        "scientific_control_states": {
                            seed: dict(
                                Counter(
                                    (rank[sample]["status"] for (sample) in (samples))
                                )
                            )
                            for ((seed, rank)) in (control_ranks.items())
                        },
                    }
                )
    return {
        "status": "CORE_ENDPOINT_PANELS",
        "cohort": cohort,
        "sample_order": samples,
        "projection_count": 48,
        "inherited_native_common_projections": 32,
        "additional_fixed_crossresource_projections": 16,
        "panels": panels,
    }
