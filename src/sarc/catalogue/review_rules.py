import itertools
import math
import random
from fractions import Fraction

CONTEXTS = ("COAD_CCLE", "LUAD_CCLE")
METHODS = ("DawnRank", "PRODIGY", "PersonaDrive")
TRANSITIONS = ("11_0_to_11_5", "11_5_to_12_0")
SEEDS = (104729, 130363, 155921, 196613, 228017)
CONFIGURATIONS = (Fraction(0), Fraction(1, 2), Fraction(1))
BUDGETS = (Fraction(0), Fraction(1, 10), Fraction(1, 4), Fraction(1, 2), Fraction(1))
DELTAS = (Fraction(1), Fraction(1, 2), Fraction(2))
FAMILIES = ("strong", "candidate")
DEFAULT = Fraction(1, 2)

def exact(value):
    if (isinstance(value, bool)):
        raise ValueError("Boolean is not a numerical value")
    if (isinstance(value, Fraction)):
        return value
    if (isinstance(value, int)):
        return Fraction(value)
    if (isinstance(value, dict)):
        numerator = value.get("numerator")
        denominator = value.get("denominator")
        if (
            (type(numerator) is not int)
            or (type(denominator) is not int)
            or (denominator <= 0)
        ):
            raise ValueError("Invalid exact carrier")
        result = Fraction(numerator, denominator)
        if (("value" in value) and (
            (not math.isfinite(value["value"])) or (float(result) != value["value"])
        )):
            raise ValueError("Inconsistent exact carrier")
        return result
    if (isinstance(value, str)):
        return Fraction(value)
    if (isinstance(value, float) and math.isfinite(value)):
        return Fraction(str(value))
    raise ValueError("Unsupported numerical value")

def to_json(value):
    if (isinstance(value, Fraction)):
        return {
            "numerator": value.numerator,
            "denominator": value.denominator,
            "value": float(value),
        }
    if (isinstance(value, dict)):
        if (any(not isinstance(key, str) for (key) in (value))):
            raise ValueError("JSON keys must be strings")
        return {key: to_json(child) for ((key, child)) in (value.items())}
    if (isinstance(value, (tuple, list))):
        return [to_json(child) for (child) in (value)]
    if (isinstance(value, (set, frozenset))):
        return [to_json(child) for (child) in (sorted(value, key = utf8))]
    if (isinstance(value, float) and (not math.isfinite(value))):
        raise ValueError("Non-finite JSON value")
    return value

def utf8(value):
    if ((not isinstance(value, str)) or (not value)):
        raise ValueError("Identity must be a non-empty string")
    return value.encode("utf-8")

def identities(values):
    values = tuple(values)
    if (len(values) != len(set(values))):
        raise ValueError("Duplicate identities")
    for (value) in (values):
        utf8(value)
    return frozenset(values)

def ranked_membership(output, eligible, k):
    if ((type(k) is not int) or (k < 1)):
        raise ValueError("Invalid leading-list scale")
    state = output["method_status"]
    genes = output["genes"]
    scores = output["scores"]
    if (state in ("FAILED", "UNAVAILABLE")):
        if ((genes is not None) or (scores is not None)):
            raise ValueError("Failed output cannot expose a ranking")
        return None, None
    if (state not in ("SUCCESS", "SUCCESS_EMPTY")):
        raise ValueError("Unknown method state")
    if (genes is None):
        raise ValueError("Successful output lacks its full roster")
    identities(genes)
    if ((state == "SUCCESS") != bool(genes)):
        raise ValueError("Method state and native roster disagree")
    if (scores is not None):
        if ((len(scores) != len(genes)) or any(
            isinstance(value, bool)
            or (not isinstance(value, (int, float, Fraction)))
            or (not math.isfinite(value))
            for (value) in (scores)
        )):
            raise ValueError("Invalid full score array")
        if (any(left < right for ((left, right)) in (zip(scores, scores[1:])))):
            raise ValueError("Native score ordering changed")
    retained = [
        (gene, index if (scores is None) else scores[index])
        for ((index, gene)) in (enumerate(genes))
        if (gene in eligible)
    ]
    membership = {}
    occupied = 0
    for (score, iterator) in (itertools.groupby(retained, key = lambda item: item[1])):
        group = list(iterator)
        available = min(len(group), max(0, k - occupied))
        if (available):
            membership.update(
                {gene: Fraction(available, len(group)) for ((gene, unused)) in (group)}
            )
        occupied += len(group)
        if (occupied >= k):
            break
    return frozenset(gene for ((gene, unused)) in (retained)), membership

def weighted_jaccard(old, new):
    universe = old.keys() | new.keys()
    denominator = sum(
        (max(old.get(gene, 0), new.get(gene, 0)) for (gene) in (universe)), Fraction(0)
    )
    if (denominator == 0):
        return None
    return (
        sum(
            (min(old.get(gene, 0), new.get(gene, 0)) for (gene) in (universe)),
            Fraction(0),
        )
        / denominator
    )

def pair_features(old, new, eligible, k = 10):
    eligible = identities(eligible)
    old_roster, old_weights = ranked_membership(old, eligible, k)
    new_roster, new_weights = ranked_membership(new, eligible, k)
    result = {
        "actionable": (old_weights is not None) and (new_weights is not None),
        "k": k,
        "old_method_status": old["method_status"],
        "new_method_status": new["method_status"],
        "_eligible": eligible,
        "_old_weights": old_weights,
        "_new_weights": new_weights,
        "_intersection": None,
        "_roster_difference": None,
        "old_roster_count": None if (old_roster is None) else len(old_roster),
        "new_roster_count": None if (new_roster is None) else len(new_roster),
        "features": None,
    }
    if (not result["actionable"]):
        return result
    common = old_roster & new_roster
    changed = old_roster ^ new_roster
    mass = sum(old_weights.values(), Fraction(0)) + sum(
        new_weights.values(), Fraction(0)
    )
    persistent_mass = sum(
        (
            abs(new_weights.get(gene, 0) - old_weights.get(gene, 0))
            for (gene) in (common)
        ),
        Fraction(0),
    )
    roster_mass = sum(
        (
            abs(new_weights.get(gene, 0) - old_weights.get(gene, 0))
            for (gene) in (changed)
        ),
        Fraction(0),
    )
    total_mass = sum(
        (
            abs(new_weights.get(gene, 0) - old_weights.get(gene, 0))
            for (gene) in (old_roster | new_roster)
        ),
        Fraction(0),
    )
    if ((persistent_mass + roster_mass) != total_mass):
        raise ValueError("Movement partition does not close")
    native_jaccard = weighted_jaccard(old_weights, new_weights)
    native_distance = Fraction(0) if (native_jaccard is None) else (1 - native_jaccard)
    unused, fixed_old = ranked_membership(old, common, k)
    unused, fixed_new = ranked_membership(new, common, k)
    fixed_jaccard = weighted_jaccard(fixed_old, fixed_new)
    fixed_distance = None if (fixed_jaccard is None) else (1 - fixed_jaccard)
    result.update(
        {
            "_intersection": common,
            "_roster_difference": changed,
            "intersection_count": len(common),
            "roster_difference_count": len(changed),
            "features": {
                "p": persistent_mass / mass if (mass) else Fraction(0),
                "r": roster_mass / mass if (mass) else Fraction(0),
                "dN": native_distance,
                "dF": fixed_distance,
                "dF_or_fallback": native_distance
                if (fixed_distance is None)
                else fixed_distance,
                "native_fractional_jaccard": native_jaccard,
                "fixed_fractional_jaccard": fixed_jaccard,
                "membership_mass": mass,
                "persistent_movement_mass": persistent_mass,
                "roster_movement_mass": roster_mass,
                "empty_pair_convention": mass == 0,
                "fixed_fallback": fixed_distance is None,
            },
        }
    )
    return result

def pair_outcome(features, positives):
    positives = identities(positives) & features["_eligible"]
    upper = Fraction(min(features["k"], len(positives)))
    hits = []
    for (name) in (("_old_weights", "_new_weights")):
        weights = features[name]
        hits.append(
            None
            if (weights is None)
            else sum(
                (value for ((gene, value)) in (weights.items()) if (gene in positives)),
                Fraction(0),
            )
        )
    old, new = hits
    lower = (Fraction(0) if (new is None) else new) - (upper if (old is None) else old)
    bound = (upper if (new is None) else new) - (Fraction(0) if (old is None) else old)
    delta = None if ((old is None) or (new is None)) else new - old
    result = {
        "old_hits": old,
        "new_hits": new,
        "delta_hits": delta,
        "delta_lower": lower,
        "delta_upper": bound,
        "positive_count": len(positives),
        "signed_persistent": None,
        "signed_roster": None,
    }
    if (delta is not None):
        old_weights = features["_old_weights"]
        new_weights = features["_new_weights"]
        for (field, genes) in ((
            ("signed_persistent", features["_intersection"]),
            ("signed_roster", features["_roster_difference"]),
        )):
            result[field] = sum(
                (
                    new_weights.get(gene, 0) - old_weights.get(gene, 0)
                    for (gene) in (genes)
                    if (gene in positives)
                ),
                Fraction(0),
            )
        if ((result["signed_persistent"] + result["signed_roster"]) != delta):
            raise ValueError("Signed reference-hit identity does not close")
    return result

def reduce_seed_pairs(sample_id, method, pairs):
    utf8(sample_id)
    if (method not in METHODS):
        raise ValueError("Unknown method")
    expected = SEEDS if (method == "PRODIGY") else (None,)
    observed = [pair["master_seed"] for (pair) in (pairs)]
    if ((len(observed) != len(expected)) or (set(observed) != set(expected))):
        raise ValueError("Exact complete seed set is required")
    ordered = sorted(pairs, key = lambda pair: expected.index(pair["master_seed"]))
    complete = all(pair["features"]["actionable"] for (pair) in (ordered))
    outcomes = [pair["outcome"] for (pair) in (ordered)]
    result = {
        "sample_id": sample_id,
        "actionable": complete,
        "expected_seed_count": len(expected),
        "complete_seed_count": sum(
            pair["features"]["actionable"] for (pair) in (ordered)
        ),
        "master_seeds": list(expected),
        "features": None,
        "outcome": {},
        "seed_records": [],
    }
    for (field) in ((
        "old_hits",
        "new_hits",
        "delta_hits",
        "delta_lower",
        "delta_upper",
        "signed_persistent",
        "signed_roster",
    )):
        values = [outcome[field] for (outcome) in (outcomes)]
        result["outcome"][field] = (
            None
            if (any(value is None for (value) in (values)))
            else sum(values, Fraction(0)) / len(expected)
        )
    if (len({outcome["positive_count"] for (outcome) in (outcomes)}) != 1):
        raise ValueError("Seed pairs have different frozen references")
    result["outcome"]["positive_count"] = outcomes[0]["positive_count"]
    if (complete):
        feature_rows = [pair["features"]["features"] for (pair) in (ordered)]
        result["features"] = {}
        for (field) in ((
            "p",
            "r",
            "dN",
            "dF",
            "dF_or_fallback",
            "native_fractional_jaccard",
            "fixed_fractional_jaccard",
            "membership_mass",
            "persistent_movement_mass",
            "roster_movement_mass",
        )):
            values = [row[field] for (row) in (feature_rows)]
            result["features"][field] = (
                None
                if (any(value is None for (value) in (values)))
                else sum(values, Fraction(0)) / len(expected)
            )
        result["features"]["empty_pair_seed_count"] = sum(
            row["empty_pair_convention"] for (row) in (feature_rows)
        )
        result["features"]["fixed_fallback_seed_count"] = sum(
            row["fixed_fallback"] for (row) in (feature_rows)
        )
        if ((
            result["outcome"]["signed_persistent"] + result["outcome"]["signed_roster"]
        ) != result["outcome"]["delta_hits"]):
            raise ValueError("Seed-averaged signed identity does not close")
    for (pair) in (ordered):
        result["seed_records"].append(
            {
                "master_seed": pair["master_seed"],
                "features": {
                    key: value
                    for ((key, value)) in (pair["features"].items())
                    if (not key.startswith("_"))
                },
                "outcome": pair["outcome"],
            }
        )
    return result

def review_priority(row, rule, parameter = DEFAULT):
    if (not row["actionable"]):
        return None
    if (rule == "direct_exact"):
        value = row["outcome"]["delta_hits"]
        if (value is None):
            raise ValueError("Actionable row lacks the reference endpoint")
        return abs(value)
    parameter = exact(parameter)
    if ((rule not in FAMILIES) or (parameter not in CONFIGURATIONS)):
        raise ValueError("Unknown review family or configuration")
    feature = row["features"]
    if ((feature is None) or any(
        feature[name] is None for (name) in (("p", "r", "dN", "dF_or_fallback"))
    )):
        raise ValueError("Actionable row lacks decision features")
    if (rule == "strong"):
        return (1 - parameter) * feature["dN"] + parameter * max(
            feature["dF_or_fallback"], feature["r"]
        )
    return (1 - parameter) * feature["p"] + parameter * feature["r"]

def ordered_allocation(rows, rule, parameter = DEFAULT):
    if ((rule not in (*FAMILIES, "direct_exact", "random", "full_inspection")) or (
        (rule in FAMILIES) and (exact(parameter) not in CONFIGURATIONS)
    )):
        raise ValueError("Unknown review family or configuration")
    identities(row["sample_id"] for (row) in (rows))
    actionable = [row for (row) in (rows) if (row["actionable"])]
    canonical = sorted((row["sample_id"] for (row) in (actionable)), key = utf8)
    if (rule == "random"):
        random.Random(20260926).shuffle(canonical)
        return {
            "order": canonical,
            "priority_tie_groups": [],
            "rule": rule,
            "parameter": None,
        }
    if (rule == "full_inspection"):
        return {
            "order": canonical,
            "priority_tie_groups": [],
            "rule": rule,
            "parameter": None,
        }
    scores = {
        row["sample_id"]: review_priority(row, rule, parameter)
        for (row) in (actionable)
    }
    ordered = sorted(canonical, key = lambda sample: (-scores[sample], utf8(sample)))
    ties = []
    for (priority, group) in (itertools.groupby(
        ordered, key = lambda sample: scores[sample]
    )):
        group = list(group)
        if (len(group) > 1):
            ties.append({"priority": priority, "sample_ids": group})
    return {
        "order": ordered,
        "priority_tie_groups": ties,
        "rule": rule,
        "parameter": exact(parameter) if (rule in FAMILIES) else None,
    }

def allocate(rows, rule, budget, parameter = DEFAULT):
    budget = exact(budget)
    if (budget not in BUDGETS):
        raise ValueError("Unfrozen budget")
    if ((rule == "full_inspection") and (budget != 1)):
        raise ValueError("Full inspection is only a 100 percent comparator")
    ordered = ordered_allocation(rows, rule, parameter)
    count = min(len(ordered["order"]), math.ceil(budget * len(ordered["order"])))
    return ordered | {
        "budget_fraction": budget,
        "actionable_count": len(ordered["order"]),
        "allocated": ordered["order"][:count],
        "allocated_count": count,
    }

def event_state(row, delta):
    delta = exact(delta)
    if (delta not in DELTAS):
        raise ValueError("Unfrozen event threshold")
    outcome = row["outcome"]
    value = outcome["delta_hits"]
    lower = outcome["delta_lower"]
    upper = outcome["delta_upper"]
    if (lower > upper):
        raise ValueError("Reversed endpoint bounds")
    if ((value is not None) and ((value < lower) or (value > upper))):
        raise ValueError("Endpoint outside attainable bounds")
    if (row["actionable"] and (value is None)):
        raise ValueError("Actionable outcome is unavailable")
    minimum = Fraction(0) if (lower <= 0 <= upper) else min(abs(lower), abs(upper))
    maximum = max(abs(lower), abs(upper))
    return {
        "event": None if (value is None) else abs(value) >= delta,
        "event_lower": int(minimum >= delta),
        "event_upper": int(maximum >= delta),
    }

def ratio(numerator, denominator):
    return None if (denominator == 0) else Fraction(numerator, denominator)

def budget_endpoints(rows, allocation, delta = 1, sample_weights = None):
    population = identities(row["sample_id"] for (row) in (rows))
    selected = identities(allocation["allocated"])
    by_id = {row["sample_id"]: row for (row) in (rows)}
    if ((not selected.issubset(population)) or any(
        not by_id[sample]["actionable"] for (sample) in (selected)
    )):
        raise ValueError("Allocation contains unavailable or foreign samples")
    weights = (
        {sample: 1 for (sample) in (population)}
        if (sample_weights is None)
        else dict(sample_weights)
    )
    if ((set(weights) != population) or any(
        (type(value) is not int) or (value < 0) for (value) in (weights.values())
    )):
        raise ValueError("Invalid frozen-population reweighting")
    full = sum(weights.values())
    actionable = sum(
        weights[row["sample_id"]] for (row) in (rows) if (row["actionable"])
    )
    inspected = sum(weights[sample] for (sample) in (selected))
    counts = {
        "event_known": 0,
        "events_known": 0,
        "non_events_known": 0,
        "events_unavailable": 0,
        "captured_known": 0,
        "missed_known": 0,
        "non_event_inspections": 0,
        "unavailable_inspections": 0,
        "events_lower": 0,
        "events_upper": 0,
        "captured_lower": 0,
        "captured_upper": 0,
        "missed_lower": 0,
        "missed_upper": 0,
    }
    for (row) in (rows):
        sample = row["sample_id"]
        weight = weights[sample]
        state = event_state(row, delta)
        inspected_here = sample in selected
        event = state["event"]
        counts["event_known" if (event is not None) else "events_unavailable"] += weight
        if (event is not None):
            counts["events_known" if (event) else "non_events_known"] += weight
            if (event):
                counts["captured_known" if (inspected_here) else "missed_known"] += (
                    weight
                )
            elif (inspected_here):
                counts["non_event_inspections"] += weight
        elif (inspected_here):
            counts["unavailable_inspections"] += weight
        for (side) in (("lower", "upper")):
            counts["events_" + side] += weight * state["event_" + side]
            counts[("captured_" if (inspected_here) else "missed_") + side] += (
                weight * state["event_" + side]
            )
    known = counts["events_unavailable"] == 0
    return {
        "full_population_count": full,
        "actionable_count": actionable,
        "allocated_sample_ids": allocation["allocated"],
        "allocated_original_count": len(selected),
        "inspected_weighted_count": inspected,
        "delta": exact(delta),
        "counts": counts,
        "captured_yield": None
        if (counts["unavailable_inspections"])
        else ratio(counts["captured_known"], full),
        "captured_yield_lower": ratio(counts["captured_lower"], full),
        "captured_yield_upper": ratio(counts["captured_upper"], full),
        "precision": None
        if (counts["unavailable_inspections"])
        else ratio(counts["captured_known"], inspected),
        "known_event_capture_rate": ratio(
            counts["captured_known"], counts["events_known"]
        ),
        "full_event_capture_rate": ratio(
            counts["captured_known"], counts["events_known"]
        )
        if (known)
        else None,
        "full_event_capture_lower": ratio(
            counts["captured_lower"], counts["captured_lower"] + counts["missed_upper"]
        ),
        "full_event_capture_upper": ratio(
            counts["captured_upper"], counts["captured_upper"] + counts["missed_lower"]
        ),
        "full_event_rate": ratio(counts["events_known"], full) if (known) else None,
        "event_status": "COMPLETE" if (known) else "UNAVAILABLE_EVENTS_RETAINED",
        "zero_event_denominator": counts["events_upper"] == 0,
        "fixed_allocation_reweighted": sample_weights is not None,
    }

def primary_strata(strata, contexts = CONTEXTS):
    contexts = tuple(contexts)
    if (
        (not contexts)
        or (len(contexts) != len(set(contexts)))
        or (not set(contexts).issubset(CONTEXTS))
    ):
        raise ValueError("Unknown training contexts")
    selected = [
        stratum
        for (stratum) in (strata)
        if (
            (stratum["context"] in contexts)
            and (stratum["support"] == "common")
            and (stratum["reference"] == "NCG6_primary_all")
            and (stratum["k"] == 10)
        )
    ]
    keys = [(row["context"], row["method"], row["transition"]) for (row) in (selected)]
    expected = set(itertools.product(contexts, METHODS, TRANSITIONS))
    if ((len(keys) != len(expected)) or (set(keys) != expected)):
        raise ValueError("Incomplete primary training-stratum population")
    return selected

def fit_review(strata, family, contexts = CONTEXTS):
    if (family not in FAMILIES):
        raise ValueError("Unknown fitted review family")
    training = primary_strata(strata, contexts)
    trials = []
    for (parameter) in (CONFIGURATIONS):
        values = []
        records = []
        for (stratum) in (training):
            allocation = allocate(stratum["rows"], family, Fraction(1, 4), parameter)
            endpoints = budget_endpoints(stratum["rows"], allocation, 1)
            values.append(endpoints["captured_yield"])
            records.append(
                {
                    "stratum_id": stratum["stratum_id"],
                    "allocation": allocation,
                    "endpoints": endpoints,
                }
            )
        objective = (
            None
            if (any(value is None for (value) in (values)))
            else sum(values, Fraction(0)) / len(values)
        )
        trials.append(
            {"parameter": parameter, "objective": objective, "records": records}
        )
    complete = all(trial["objective"] is not None for (trial) in (trials))
    selected = (
        sorted(
            trials,
            key = lambda trial: (
                -trial["objective"],
                trial["parameter"] != DEFAULT,
                str(trial["parameter"]).encode("utf-8"),
            ),
        )[0]["parameter"]
        if (complete)
        else DEFAULT
    )
    return {
        "family": family,
        "parameter": selected,
        "status": "FITTED_PRIMARY_DEVELOPMENT" if (complete) else "DEFAULT_NOT_TUNED",
        "training_contexts": list(contexts),
        "training_stratum_count": len(training),
        "trials": trials,
    }

def curves(stratum, fitted, deltas = DELTAS):
    if (set(fitted) != set(FAMILIES)):
        raise ValueError("Both equally budgeted review families are required")
    result = []
    for (rule) in ((*FAMILIES, "direct_exact", "random", "full_inspection")):
        parameter = fitted[rule]["parameter"] if (rule in FAMILIES) else DEFAULT
        for (budget) in ((Fraction(1),) if (rule == "full_inspection") else BUDGETS):
            allocation = allocate(stratum["rows"], rule, budget, parameter)
            for (delta) in (deltas):
                result.append(
                    {
                        "rule": rule,
                        "parameter": parameter if (rule in FAMILIES) else None,
                        "allocation": allocation,
                        "endpoints": budget_endpoints(
                            stratum["rows"], allocation, delta
                        ),
                    }
                )
    return {key: value for ((key, value)) in (stratum.items()) if (key != "rows")} | {
        "curves": result
    }

def loco_review(strata):
    primary_strata(strata)
    folds = []
    for (held) in (CONTEXTS):
        training = tuple(context for (context) in (CONTEXTS) if (context != held))
        fitted = {
            family: fit_review(strata, family, training) for (family) in (FAMILIES)
        }
        held_strata = [
            stratum for (stratum) in (strata) if (stratum["context"] == held)
        ]
        folds.append(
            {
                "held_context": held,
                "training_contexts": list(training),
                "fitted": fitted,
                "evaluation": [curves(stratum, fitted) for (stratum) in (held_strata)],
            }
        )
    return {
        "status": "TWO_CONTEXT_DEVELOPMENT_DIAGNOSTIC_NOT_PROTECTED_VALIDATION",
        "folds": folds,
    }

def evaluate_review(strata):
    fitted = {family: fit_review(strata, family) for (family) in (FAMILIES)}
    return {
        "status": "DEVELOPMENT_REVIEW_COMPARISON",
        "fitted": fitted,
        "strata": [curves(stratum, fitted) for (stratum) in (strata)],
    }

def reweight_allocation(rows, allocation, indices, delta = 1):
    sample_order = [row["sample_id"] for (row) in (rows)]
    identities(sample_order)
    n = len(sample_order)
    if ((len(indices) != n) or any(
        (type(index) is not int) or (index < 0) or (index >= n) for (index) in (indices)
    )):
        raise ValueError("Resampling must use the frozen complete population")
    weights = {sample: 0 for (sample) in (sample_order)}
    for (index) in (indices):
        weights[sample_order[index]] += 1
    return budget_endpoints(rows, allocation, delta, weights)
