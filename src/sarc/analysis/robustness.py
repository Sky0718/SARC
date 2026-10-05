import math
from fractions import Fraction

METHODS = ("DawnRank", "PRODIGY", "PersonaDrive")
SEEDS = (104729, 130363, 155921, 196613, 228017)
POLICIES = ("FIXED_V0", "FIXED_V2", "RC", "STABLE_ONLY_3")
CONTRASTS = (
    ("FIXED_V2", "FIXED_V0"),
    ("RC", "FIXED_V2"),
    ("STABLE_ONLY_3", "FIXED_V2"),
)
MASSES = ("returned", "assessed", "unassessed", "supported")

def seeds(method):
    if (method not in METHODS):
        raise ValueError("Unknown method")
    return SEEDS if (method == "PRODIGY") else (None,)

def mean(values):
    return (
        None
        if (not values or any((value is None for (value) in (values))))
        else sum(values) / len(values)
    )

def validate_split(split, production = True):
    models = split["models"]
    names = [model["sample_ID"] for (model) in (models)]
    donors = split["ordered_donors"]
    if (len(names) != len(set(names)) or len(donors) != len(set(donors))):
        raise ValueError("Duplicate model or donor identity")
    if (set(donors) != {model["study_donor_id"] for (model) in (models)}):
        raise ValueError("Donor population mismatch")
    for (donor) in (donors):
        states = {
            model["split"] for (model) in (models) if (model["study_donor_id"] == donor)
        }
        if (len(states) != 1 or not states <= {"development", "held_out"}):
            raise ValueError("Related models cross donor splits")
    counts = {
        label: sum(
            (
                any(
                    (
                        model["study_donor_id"] == donor and model["split"] == label
                        for (model) in (models)
                    )
                )
                for (donor) in (donors)
            )
        )
        for (label) in (("development", "held_out"))
    }
    if (production and (
        len(models),
        len(donors),
        counts["development"],
        counts["held_out"],
    ) != (
        85,
        83,
        25,
        58,
    )):
        raise ValueError("Original population changed")
    if (production and sum((model["split"] == "held_out" for (model) in (models))) != 59):
        raise ValueError("Original held-out model population changed")
    return {model["sample_ID"]: model for (model) in (models)}

def accumulate_ledger(rows, split, production = True):
    models = validate_split(split, production)
    groups = {}
    total_rows = 0
    selected_rows = 0
    for (row) in (rows):
        total_rows += 1
        if ((row["normal_reference_variant"], row["support"], row["offered_k"]) != (
            "all34",
            "native",
            10,
        ) or row["policy_id"] not in POLICIES):
            continue
        selected_rows += 1
        model = models.get(row["model_id"])
        if (model is None or (row["donor_id"], row["split"]) != (
            model["study_donor_id"],
            model["split"],
        )):
            raise ValueError("Candidate-ledger model/donor/split identity mismatch")
        method, seed, policy = (row["method"], row["seed"], row["policy_id"])
        if (seed not in seeds(method)):
            raise ValueError("Unexpected original algorithm seed")
        key = (row["model_id"], method, seed, policy)
        state = row["prediction_state"]
        if (state not in ("unavailable", "available_empty", "available_nonempty")):
            raise ValueError("Unknown prediction state")
        value = groups.setdefault(
            key,
            {
                "state": state,
                "genes": set(),
                "sentinel": False,
                "returned": Fraction(0),
                "assessed": Fraction(0),
                "supported": Fraction(0),
            },
        )
        if (value["state"] != state):
            raise ValueError("Mixed prediction states within one seed")
        gene = row["gene_id"]
        if (gene is None):
            if (
                state == "available_nonempty"
                or value["sentinel"]
                or value["genes"]
                or (row["weight_numerator"] is not None)
                or (row["weight_denominator"] is not None)
                or (row["assay_state"] != "not_applicable")
                or (row["binary_call"] is not None)
            ):
                raise ValueError("Invalid empty or unavailable sentinel")
            value["sentinel"] = True
            continue
        if (state != "available_nonempty" or value["sentinel"] or gene in value["genes"]):
            raise ValueError("Duplicate gene or incompatible prediction state")
        numerator, denominator = (row["weight_numerator"], row["weight_denominator"])
        if (
            type(numerator) is not int
            or type(denominator) is not int
            or denominator <= 0
        ):
            raise ValueError("Invalid exact fractional weight")
        weight = Fraction(numerator, denominator)
        if (weight <= 0 or weight > 1):
            raise ValueError("Fractional membership outside (0,1]")
        value["genes"].add(gene)
        value["returned"] += weight
        call = row["binary_call"]
        if (row["assay_state"] == "eligible_measured"):
            if (call not in (0, 1)):
                raise ValueError("Measured accepted binary call is not binary")
            value["assessed"] += weight
            value["supported"] += weight * int(call)
        elif (row["assay_state"] != "not_measured" or call is not None):
            raise ValueError("Unmeasured accepted call must remain missing")
    expected = {
        (sample, method, seed, policy)
        for (sample) in (models)
        for (method) in (METHODS)
        for (seed) in (seeds(method))
        for (policy) in (POLICIES)
    }
    if (set(groups) != expected):
        raise ValueError("Missing or extra model/method/seed/policy detail")
    for (value) in (groups.values()):
        if (value["returned"] > 10):
            raise ValueError("Returned mass exceeds original offered k10")
        if (value["state"] != "available_nonempty" and (not value["sentinel"])):
            raise ValueError("Missing empty/unavailable sentinel")
        value["unassessed"] = value["returned"] - value["assessed"]
        value.pop("genes")
        value.pop("sentinel")
        if (value["state"] == "unavailable"):
            for (field) in (MASSES):
                value[field] = None
    return (
        groups,
        {
            "total_accepted_ledger_rows": total_rows,
            "selected_derived_rows": selected_rows,
            "selected_seed_policy_groups": len(groups),
            "raw_screen_labels_opened": False,
            "seed_reduction_precedes_donor_reduction": True,
        },
    )

def thin_seed(latest, stable):
    if (latest["state"] == "unavailable" or stable["state"] == "unavailable"):
        return {
            "status": "UNAVAILABLE_PREDICTION",
            "q": None,
            **{
                field: None
                for (field) in (
                    (
                        "latest_returned",
                        "stable_returned",
                        "latest_assessed",
                        "stable_assessed",
                        "latest_unassessed",
                        "stable_unassessed",
                        "latest_vacancies",
                        "stable_vacancies",
                        "expected_returned",
                        "expected_assessed",
                        "expected_unassessed",
                        "expected_supported",
                        "stable_supported",
                        "stable_minus_expected_supported",
                        "stable_minus_latest_supported",
                        "complete_difference_lower",
                        "complete_difference_upper",
                    )
                )
            },
        }
    r_l, r_s = (latest["returned"], stable["returned"])
    if (not 0 <= r_s <= r_l <= 10):
        raise ValueError("Stable-only/latest returned mass relation impossible")
    q = r_s / r_l if (r_l) else Fraction(0)
    h_l, h_s, u_l, u_s = (
        latest["supported"],
        stable["supported"],
        latest["unassessed"],
        stable["unassessed"],
    )
    return {
        "status": "AVAILABLE",
        "q": q,
        "latest_returned": r_l,
        "stable_returned": r_s,
        "latest_assessed": latest["assessed"],
        "stable_assessed": stable["assessed"],
        "latest_unassessed": u_l,
        "stable_unassessed": u_s,
        "latest_vacancies": 10 - r_l,
        "stable_vacancies": 10 - r_s,
        "expected_returned": q * r_l,
        "expected_assessed": q * latest["assessed"],
        "expected_unassessed": q * u_l,
        "expected_supported": q * h_l,
        "stable_supported": h_s,
        "stable_minus_expected_supported": h_s - q * h_l,
        "stable_minus_latest_supported": h_s - h_l,
        "complete_difference_lower": h_s - q * (h_l + u_l),
        "complete_difference_upper": h_s + u_s - q * h_l,
    }

def donor_values(model_values, models, donor_order):
    return [
        mean(
            [
                model_values[model["sample_ID"]]
                for (model) in (models)
                if (model["study_donor_id"] == donor)
            ]
        )
        for (donor) in (donor_order)
    ]

def influence(values, donor_order):
    if (len(values) != len(donor_order) or len(values) < 2):
        raise ValueError("Invalid donor influence population")
    point = mean(values)
    records = [
        {"omitted_donor": donor, "mean": mean(values[:index] + values[index + 1 :])}
        for ((index, donor)) in (enumerate(donor_order))
    ]
    if (point is None):
        return {
            "status": "UNAVAILABLE_FULL_POPULATION",
            "full_mean": None,
            "full_donors": len(values),
            "available_donors": sum((value is not None for (value) in (values))),
            "leave_one_donor_out": records,
            "range": None,
            "largest_absolute_change": None,
            "sign_stable": None,
        }
    signs = lambda value: (value > 0) - (value < 0)
    return {
        "status": "AVAILABLE",
        "full_mean": point,
        "full_donors": len(values),
        "available_donors": len(values),
        "leave_one_donor_out": records,
        "range": [
            min((row["mean"] for (row) in (records))),
            max((row["mean"] for (row) in (records))),
        ],
        "largest_absolute_change": max(
            (abs(row["mean"] - point) for (row) in (records))
        ),
        "sign_stable": all((signs(row["mean"]) == signs(point) for (row) in (records))),
    }

def quantile_linear(values, probability):
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    low, high = (math.floor(position), math.ceil(position))
    return ordered[low] + (ordered[high] - ordered[low]) * (position - low)

def joint_sensitivity(columns, donor_order, bootstrap, expected_draws = 2000):
    unavailable = lambda reason: {
        "status": "UNAVAILABLE_ENTIRE_NINE_CONTRAST_FAMILY",
        "reason": reason,
        "contrasts_retained": list(columns),
        "dropped_contrasts": [],
        "new_draws_generated": False,
    }
    if (len(columns) != 9):
        return unavailable("Exactly nine original contrasts required")
    if (bootstrap is None):
        return unavailable("Accepted paired bootstrap carrier unavailable")
    try:
        shared = bootstrap["splits"]["held_out"]
        draws = shared["indices"]
        if (
            shared["donor_order"] != donor_order
            or len(draws) != expected_draws
            or any(
                (
                    len(draw) != len(donor_order)
                    or any(
                        (
                            type(index) is not int or not 0 <= index < len(donor_order)
                            for (index) in (draw)
                        )
                    )
                    for (draw) in (draws)
                )
            )
        ):
            return unavailable(
                "Accepted paired donor order or bootstrap index geometry invalid"
            )
    except (KeyError, TypeError):
        return unavailable("Accepted paired bootstrap schema unavailable")
    points, scales = ({}, {})
    n = len(donor_order)
    for (key, values) in (columns.items()):
        if (len(values) != n or any(
            (value is None or not math.isfinite(float(value)) for (value) in (values))
        )):
            return unavailable("Incomplete paired donor values: " + key)
        point = math.fsum((float(value) for (value) in (values))) / n
        scale = (
            math.sqrt(
                math.fsum(((float(value) - point) ** 2 for (value) in (values)))
                / (n * (n - 1))
            )
            if (n > 1)
            else 0
        )
        if (not math.isfinite(scale) or scale <= 0):
            return unavailable(
                "Zero or undefined original-sample donor standard error: " + key
            )
        points[key], scales[key] = (point, scale)
    maxima = []
    for (draw) in (draws):
        maxima.append(
            max(
                (
                    abs(
                        math.fsum((float(columns[key][index]) for (index) in (draw)))
                        / n
                        - points[key]
                    )
                    / scales[key]
                    for (key) in (columns)
                )
            )
        )
    critical = quantile_linear(maxima, 0.95)
    intervals = {
        key: {
            "mean": points[key],
            "original_sample_standard_error": scales[key],
            "approximate_simultaneous_interval": [
                points[key] - critical * scales[key],
                points[key] + critical * scales[key],
            ],
        }
        for (key) in (columns)
    }
    return {
        "status": "AVAILABLE_APPROXIMATE_POSTHOC_SENSITIVITY",
        "family_size": 9,
        "draws": len(draws),
        "quantile": 0.95,
        "quantile_interpolation": "linear",
        "scale": "original_sample_donor_standard_error_fixed_across_bootstrap_draws",
        "critical_maximum_absolute_standardised_deviation": critical,
        "intervals": intervals,
        "maxima_in_original_draw_order": maxima,
        "new_draws_generated": False,
        "formal_confirmatory_p_values_claimed": False,
        "exact_finite_sample_coverage_claimed": False,
    }

def analyse(groups, split, bootstrap):
    models = split["models"]
    held_models = [model for (model) in (models) if (model["split"] == "held_out")]
    held_set = {model["study_donor_id"] for (model) in (held_models)}
    held_donors = [
        donor for (donor) in (split["ordered_donors"]) if (donor in held_set)
    ]
    seed_thinning, model_thinning, thinning_summaries, contrasts, columns = (
        [],
        [],
        [],
        [],
        {},
    )
    for (method) in (METHODS):
        per_model = {}
        for (model) in (models):
            sample = model["sample_ID"]
            values = []
            for (seed) in (seeds(method)):
                value = thin_seed(
                    groups[sample, method, seed, "FIXED_V2"],
                    groups[sample, method, seed, "STABLE_ONLY_3"],
                )
                seed_thinning.append(
                    {
                        "model_id": sample,
                        "donor_id": model["study_donor_id"],
                        "split": model["split"],
                        "method": method,
                        "seed": seed,
                        **value,
                    }
                )
                values.append(value)
            fields = [key for (key) in (values[0]) if (key not in ("status", "q"))]
            averaged = {
                field: mean([value[field] for (value) in (values)])
                for (field) in (fields)
            }
            per_model[sample] = averaged
            model_thinning.append(
                {
                    "model_id": sample,
                    "donor_id": model["study_donor_id"],
                    "split": model["split"],
                    "method": method,
                    "seeds": list(seeds(method)),
                    "status": "AVAILABLE"
                    if (all((value["status"] == "AVAILABLE" for (value) in (values))))
                    else "UNAVAILABLE_PREDICTION",
                    **averaged,
                }
            )
        for (label) in (("held_out", "development", "all")):
            population = [
                model
                for (model) in (models)
                if (label == "all" or model["split"] == label)
            ]
            donors = [
                donor
                for (donor) in (split["ordered_donors"])
                if (any((model["study_donor_id"] == donor for (model) in (population))))
            ]
            summary = {
                field: mean(
                    donor_values(
                        {
                            sample: value[field]
                            for ((sample, value)) in (per_model.items())
                        },
                        population,
                        donors,
                    )
                )
                for (field) in (fields)
            }
            thinning_summaries.append(
                {
                    "method": method,
                    "stratum": label,
                    "full_models": len(population),
                    "full_donors": len(donors),
                    **summary,
                }
            )
        for (left, right) in (CONTRASTS):
            differences = {}
            for (model) in (held_models):
                sample = model["sample_ID"]
                left_value = mean(
                    [
                        groups[sample, method, seed, left]["supported"]
                        for (seed) in (seeds(method))
                    ]
                )
                right_value = mean(
                    [
                        groups[sample, method, seed, right]["supported"]
                        for (seed) in (seeds(method))
                    ]
                )
                differences[sample] = (
                    None
                    if (left_value is None or right_value is None)
                    else left_value - right_value
                )
            values = donor_values(differences, held_models, held_donors)
            key = method + "|" + left + "-" + right
            columns[key] = values
            contrasts.append(
                {
                    "contrast_id": key,
                    "method": method,
                    "left": left,
                    "right": right,
                    "donor_order": held_donors,
                    "donor_values": values,
                    **influence(values, held_donors),
                }
            )
    return {
        "status": "NEW_C_POSTHOC_DIAGNOSTICS_REQUIRE_FIRST_ACCEPTANCE",
        "scope": {
            "normal_reference": "all34",
            "support": "native",
            "offered_k": 10,
            "models_retained": len(models),
            "donors_retained": len(split["ordered_donors"]),
            "held_out_models": len(held_models),
            "held_out_donors": len(held_donors),
            "primary_contrasts_retained": 9,
            "other_original_normal_references_k_values_and_results_unchanged": True,
        },
        "thinning_estimand": "per_model_method_seed_outcome_independent_expected_mass_then_complete_seed_mean_then_equal_donor_mean",
        "thinning_is_not_realised_fixed_cardinality": True,
        "unknown_support_bounds": "conservative_marginal_not_shared_gene_sharp",
        "seed_thinning": seed_thinning,
        "model_thinning": model_thinning,
        "thinning_summaries": thinning_summaries,
        "donor_influence": contrasts,
        "joint_uncertainty": joint_sensitivity(columns, held_donors, bootstrap),
        "posthoc": True,
        "original_primary_analyses_overwritten": False,
        "new_ranks_or_fits": False,
    }
