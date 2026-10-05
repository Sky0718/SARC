import itertools
import json
import re
from collections import Counter
from fractions import Fraction
from . import robustness as HELPER

METHODS = ("DawnRank", "PRODIGY", "PersonaDrive")
BUDGETS = (1, 5, 10, 20)
SUBSETS = ("held_out", "development", "all")
POLICIES = ("FIXED_V0", "FIXED_V1", "FIXED_V2")
CONTRASTS = (
    ("FIXED_V1", "FIXED_V0"),
    ("FIXED_V2", "FIXED_V1"),
    ("FIXED_V2", "FIXED_V0"),
)
MASK = re.compile(
    '("continuous_LFC"\\s*:\\s*)(?:null|-?(?:0|[1-9]\\d*)(?:\\.\\d+)?(?:[eE][+-]?\\d+)?)(?=\\s*[,}])'
)

def unique_pairs(pairs):
    result = {}
    for (key, value) in (pairs):
        if (key in result):
            raise ValueError("Duplicate JSON key")
        result[key] = value
    return result

def masked_row(line):
    masked, count = MASK.subn(lambda match: match.group(1) + "null", line)
    if (count != 1):
        raise ValueError("Exactly one continuous token must be erased")
    row = json.loads(masked, object_pairs_hook = unique_pairs)
    if (row.pop("continuous_LFC") is not None):
        raise ValueError("Continuous value survived masking")
    row.pop("continuous_assay_state", None)
    return row

def selected(row):
    return (
        row["normal_reference_variant"] == "all34"
        and row["support"] == "native"
        and (row["offered_k"] in BUDGETS)
        and (row["policy_id"] in POLICIES)
    )

def accumulate(rows, split, production = True):
    models = HELPER.validate_split(split, production)
    groups, labels = ({}, {})
    total, kept = (0, 0)
    for (row) in (rows):
        total += 1
        if (not selected(row)):
            continue
        kept += 1
        sample = row["model_id"]
        if (sample not in models):
            raise ValueError("Ledger model_id is not a literal sample_ID")
        model = models[sample]
        if ((row["donor_id"], row["split"]) != (model["study_donor_id"], model["split"])):
            raise ValueError("Donor or split association changed")
        method, seed, policy, budget = (
            row["method"],
            row["seed"],
            row["policy_id"],
            row["offered_k"],
        )
        if (type(budget) is not int or seed not in HELPER.seeds(method)):
            raise ValueError("Budget or original seed invalid")
        key = (sample, method, seed, policy, budget)
        state = row["prediction_state"]
        if (state not in ("unavailable", "available_empty", "available_nonempty")):
            raise ValueError("Unknown prediction state")
        group = groups.setdefault(
            key, {"state": state, "weights": {}, "sentinel": False}
        )
        if (group["state"] != state):
            raise ValueError("Mixed prediction state")
        gene = row["gene_id"]
        if (gene is None):
            if (
                state == "available_nonempty"
                or group["sentinel"]
                or group["weights"]
                or (row["weight_numerator"] is not None)
                or (row["weight_denominator"] is not None)
                or (row["assay_state"] != "not_applicable")
                or (row["binary_call"] is not None)
            ):
                raise ValueError("Invalid prediction sentinel")
            group["sentinel"] = True
            continue
        if (
            not isinstance(gene, str)
            or not gene
            or state != "available_nonempty"
            or group["sentinel"]
            or (gene in group["weights"])
        ):
            raise ValueError("Duplicate gene or incompatible gene state")
        numerator, denominator = (row["weight_numerator"], row["weight_denominator"])
        if (
            type(numerator) is not int
            or type(denominator) is not int
            or denominator <= 0
        ):
            raise ValueError("Exact rational weight malformed")
        weight = Fraction(numerator, denominator)
        if (not 0 < weight <= 1):
            raise ValueError("Gene allocation outside original capacity")
        call = row["binary_call"]
        if (row["assay_state"] == "eligible_measured"):
            if (type(call) is not int or call not in (0, 1)):
                raise ValueError("Measured final Table5 call is not integer binary")
        elif (row["assay_state"] != "not_measured" or call is not None):
            raise ValueError("Missing Table5 state malformed")
        identity = (sample, gene)
        if (identity in labels and labels[identity] != call):
            raise ValueError("Repeated model/gene call or missingness conflict")
        labels[identity] = call
        group["weights"][gene] = weight
    expected = {
        (sample, method, seed, policy, budget)
        for (sample) in (models)
        for (method) in (METHODS)
        for (seed) in (HELPER.seeds(method))
        for (policy) in (POLICIES)
        for (budget) in (BUDGETS)
    }
    if (set(groups) != expected):
        raise ValueError("Incomplete original model/method/seed/release/budget family")
    for (key, group) in (groups.items()):
        if (sum(group["weights"].values(), Fraction(0)) > key[-1]):
            raise ValueError("Returned mass exceeds offered budget")
        if (group["state"] != "available_nonempty" and (not group["sentinel"])):
            raise ValueError("Missing prediction sentinel")
    return (
        groups,
        labels,
        {
            "rows_streamed": total,
            "selected_rows": kept,
            "seed_policy_groups": len(groups),
            "distinct_selected_model_gene_calls": len(labels),
        },
    )

def mean_allocations(states):
    if (any((state["state"] == "unavailable" for (state) in (states)))):
        return None
    genes = set().union(*(state["weights"] for (state) in (states)))
    return {
        gene: sum(
            (state["weights"].get(gene, Fraction(0)) for (state) in (states)),
            Fraction(0),
        )
        / len(states)
        for (gene) in (sorted(genes))
    }

def masses(weights, labels, budget):
    if (weights is None):
        return None
    returned = sum(weights.values(), Fraction(0))
    assessed = sum(
        (weight for ((gene, weight)) in (weights.items()) if (labels[gene] is not None)),
        Fraction(0),
    )
    supported = sum(
        (
            weight * labels[gene]
            for ((gene, weight)) in (weights.items())
            if (labels[gene] is not None)
        ),
        Fraction(0),
    )
    return {
        "returned": returned,
        "assessed": assessed,
        "unassessed": returned - assessed,
        "observed_supported": supported,
        "vacancy": budget - returned,
    }

def contrast(new, old, labels):
    if (new is None or old is None):
        return {
            "status": "UNAVAILABLE_PREDICTION",
            "observed_component": None,
            "lower": None,
            "upper": None,
            "coefficients": None,
        }
    coefficients = {
        gene: new.get(gene, Fraction(0)) - old.get(gene, Fraction(0))
        for (gene) in (sorted(set(new) | set(old)))
    }
    observed = sum(
        (
            value * labels[gene]
            for ((gene, value)) in (coefficients.items())
            if (labels[gene] is not None)
        ),
        Fraction(0),
    )
    lower = observed + sum(
        (
            min(value, Fraction(0))
            for ((gene, value)) in (coefficients.items())
            if (labels[gene] is None)
        ),
        Fraction(0),
    )
    upper = observed + sum(
        (
            max(value, Fraction(0))
            for ((gene, value)) in (coefficients.items())
            if (labels[gene] is None)
        ),
        Fraction(0),
    )
    if (lower > upper):
        raise ValueError("Invalid marginal extrema")
    status = (
        "DEFINITE_STRICT_LOSS"
        if (upper < 0)
        else "DEFINITE_STRICT_GAIN"
        if (lower > 0)
        else "GUARANTEED_ZERO"
        if (lower == upper == 0)
        else "SIGN_UNRESOLVED_OUTER_INTERVAL"
    )
    return {
        "status": status,
        "observed_component": observed,
        "lower": lower,
        "upper": upper,
        "coefficients": coefficients,
        "point_identified": lower == upper,
        "zero_attainability_claimed": lower == upper == 0,
        "labels": {gene: labels[gene] for (gene) in (coefficients)},
    }

def event_flags(lower, upper):
    return {
        "strict_loss": (upper < 0, lower < 0),
        "strict_gain": (lower > 0, upper > 0),
        "at_least_one_hit_loss": (upper <= -1, lower <= -1),
    }

def weights_for(models):
    counts = Counter((model["study_donor_id"] for (model) in (models)))
    return {
        model["sample_ID"]: Fraction(1, len(counts) * counts[model["study_donor_id"]])
        for (model) in (models)
    }

def summarise(records, models, budget):
    weights = weights_for(models)
    if (
        set(weights) != {row["sample_ID"] for (row) in (records)}
        or sum(weights.values(), Fraction(0)) != 1
    ):
        raise ValueError("Summary population or donor weights invalid")
    available = [row for (row) in (records) if (row["lower"] is not None)]
    unavailable = [row["sample_ID"] for (row) in (records) if (row["lower"] is None)]
    events = {}
    for (event) in (("strict_loss", "strict_gain", "at_least_one_hit_loss")):
        definite = [
            row
            for (row) in (available)
            if (event_flags(row["lower"], row["upper"])[event][0])
        ]
        possible = [
            row
            for (row) in (available)
            if (event_flags(row["lower"], row["upper"])[event][1])
        ]
        events[event] = {
            "definite_count": len(definite),
            "possible_count": len(possible),
            "available_contribution_lower": sum(
                (weights[row["sample_ID"]] for (row) in (definite)), Fraction(0)
            ),
            "available_contribution_upper": sum(
                (weights[row["sample_ID"]] for (row) in (possible)), Fraction(0)
            ),
        }
        events[event]["full_population_bounds"] = (
            None
            if (unavailable)
            else [
                events[event]["available_contribution_lower"],
                events[event]["available_contribution_upper"],
            ]
        )
    thresholds = sorted(
        {Fraction(-budget), Fraction(-1), Fraction(0), Fraction(1), Fraction(budget)}
        | {row[field] for (row) in (available) for (field) in (("lower", "upper"))}
    )
    cdf = []
    for (threshold) in (thresholds):
        low = sum(
            (
                weights[row["sample_ID"]]
                for (row) in (available)
                if (row["upper"] <= threshold)
            ),
            Fraction(0),
        )
        high = sum(
            (
                weights[row["sample_ID"]]
                for (row) in (available)
                if (row["lower"] <= threshold)
            ),
            Fraction(0),
        )
        cdf.append(
            {
                "threshold": threshold,
                "available_contribution_lower": low,
                "available_contribution_upper": high,
                "full_population_bounds": None if (unavailable) else [low, high],
            }
        )
    categories = {}
    for (status) in ((
        "DEFINITE_STRICT_LOSS",
        "DEFINITE_STRICT_GAIN",
        "GUARANTEED_ZERO",
        "SIGN_UNRESOLVED_OUTER_INTERVAL",
        "UNAVAILABLE_PREDICTION",
    )):
        matches = [row for (row) in (records) if (row["status"] == status)]
        categories[status] = {
            "models": len(matches),
            "donor_weighted_model_mass": sum(
                (weights[row["sample_ID"]] for (row) in (matches)), Fraction(0)
            ),
        }
    return {
        "status": "FULL_FINITE_POPULATION_AVAILABLE"
        if (not unavailable)
        else "FULL_FINITE_POPULATION_UNAVAILABLE",
        "models": len(models),
        "donors": len({model["study_donor_id"] for (model) in (models)}),
        "model_weights": weights,
        "unavailable_models": unavailable,
        "unavailable_mass": sum(
            (weights[sample] for (sample) in (unavailable)), Fraction(0)
        ),
        "categories": categories,
        "event_bounds": events,
        "cdf_definition": "F(t)=weighted model mass with Delta<=t",
        "cdf": cdf,
    }

def analyse(groups, labels, split):
    records, summaries = ([], [])
    for (model) in (split["models"]):
        sample = model["sample_ID"]
        calls = {
            gene: value
            for (((sample_key, gene), value)) in (labels.items())
            if (sample_key == sample)
        }
        for (method) in (METHODS):
            for (budget) in (BUDGETS):
                states = {
                    policy: [
                        groups[sample, method, seed, policy, budget]
                        for (seed) in (HELPER.seeds(method))
                    ]
                    for (policy) in (POLICIES)
                }
                allocations = {
                    policy: mean_allocations(value)
                    for ((policy, value)) in (states.items())
                }
                triplet = []
                for (new, old) in (CONTRASTS):
                    value = contrast(allocations[new], allocations[old], calls)
                    row = {
                        "sample_ID": sample,
                        "association_metadata": model,
                        "method": method,
                        "budget": budget,
                        "new": new,
                        "old": old,
                        "seeds": list(HELPER.seeds(method)),
                        "seed_states": {
                            policy: [state["state"] for (state) in (states[policy])]
                            for (policy) in ((old, new))
                        },
                        "old_mass": masses(allocations[old], calls, budget),
                        "new_mass": masses(allocations[new], calls, budget),
                        **value,
                    }
                    if (value["lower"] is not None and (
                        not -budget <= value["lower"] <= value["upper"] <= budget
                    )):
                        raise ValueError("Functional contrast exceeds offered capacity")
                    records.append(row)
                    triplet.append(value)
                if (all((row["coefficients"] is not None for (row) in (triplet)))):
                    genes = set().union(*(row["coefficients"] for (row) in (triplet)))
                    for (gene) in (genes):
                        if (triplet[0]["coefficients"].get(gene, 0) + triplet[1][
                            "coefficients"
                        ].get(gene, 0) != triplet[2]["coefficients"].get(gene, 0)):
                            raise ValueError(
                                "Release coefficient vector does not telescope"
                            )
                    if (
                        triplet[0]["observed_component"]
                        + triplet[1]["observed_component"]
                        != triplet[2]["observed_component"]
                    ):
                        raise ValueError("Known contribution does not telescope")
    for (subset) in (SUBSETS):
        models = [
            model
            for (model) in (split["models"])
            if (subset == "all" or model["split"] == subset)
        ]
        samples = {model["sample_ID"] for (model) in (models)}
        for (method, budget, (new, old)) in (itertools.product(
            METHODS, BUDGETS, CONTRASTS
        )):
            family = [
                row
                for (row) in (records)
                if (
                    row["sample_ID"] in samples
                    and (row["method"], row["budget"], row["new"], row["old"])
                    == (method, budget, new, old)
                )
            ]
            summaries.append(
                {
                    "subset": subset,
                    "method": method,
                    "budget": budget,
                    "new": new,
                    "old": old,
                    "primary_design_cell": subset == "held_out" and budget == 10,
                    **summarise(family, models, budget),
                }
            )
    return {
        "status": "NEW_FINITE_CASE_RELEASE_DIAGNOSTIC_REQUIRES_FIRST_INDEPENDENT_ACCEPTANCE",
        "exposed_exploratory": True,
        "numeric_carrier_opened": True,
        "old_science_replayed": False,
        "p_values": False,
        "population_equivalence_claimed": False,
        "shared_unknowns_jointly_extremised_across_groups": False,
        "three_pairs_are_two_adjacent_steps_plus_their_composition": True,
        "model_records": records,
        "summaries": summaries,
    }
