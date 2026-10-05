import math
import time
from collections import Counter, defaultdict
from fractions import Fraction as F

SEEDS = (104729, 130363, 155921, 196613, 228017)
METHODS = ("PRODIGY", "DawnRank")
CELLS = ("G0H0", "G1H0", "G0H1", "G1H1")

def combine(terms):
    output = defaultdict(F)
    for (scale, vector) in (terms):
        for (key, value) in (vector.items()):
            output[key] += F(scale) * value
    return {key: value for (key, value) in (output.items()) if (value)}

def direction(lower, upper):
    if (upper < 0):
        return "negative"
    if (lower > 0):
        return "positive"
    return "zero" if (lower == upper == 0) else "unresolved"

def bound(vector, calls):
    known, low, high = F(), F(), F()
    unknown = 0
    for (key, value) in (vector.items()):
        state, label = calls[key]
        if (state == "MEASURED"):
            known += value * label
        else:
            low += min(value, 0)
            high += max(value, 0)
            unknown += 1
    return {
        "known": known,
        "lower": known + low,
        "upper": known + high,
        "sign": direction(known + low, known + high),
        "unknown_active_identities": unknown,
    }

def sample_vector(sample, weights):
    return {(sample, gene): value for (gene, value) in (weights.items()) if (value)}

def summary(models, vectors, calls):
    groups = defaultdict(list)
    for (model) in (models):
        groups[model["donor_id"]].append(model)
    model_rows, donor_rows, donor_vectors = [], [], {}
    for (model) in (models):
        vector = sample_vector(model["sample_id"], vectors[model["model_id"]])
        model_rows.append(
            {key: model[key] for (key) in (("model_id", "sample_id", "donor_id"))}
            | bound(vector, calls)
        )
    for (donor, members) in (sorted(groups.items())):
        vector = combine(
            [
                (
                    F(1, len(members)),
                    sample_vector(model["sample_id"], vectors[model["model_id"]]),
                )
                for (model) in (members)
            ]
        )
        donor_vectors[donor] = vector
        donor_rows.append({"donor_id": donor, "models": len(members), **bound(vector, calls)})
    population = bound(
        combine([(F(1, len(groups)), vector) for (vector) in (donor_vectors.values())]), calls
    )
    population["donor_directions"] = dict(Counter(row["sign"] for (row) in (donor_rows)))
    return {"population": population, "donors": donor_rows, "models": model_rows}, donor_vectors

def is_opposed(first, second):
    return (first, second) in (("negative", "positive"), ("positive", "negative"))

def make_lattice(first, second, calls):
    keys = sorted(set(first) | set(second))
    known = [F(), F()]
    unknown = []
    for (key) in (keys):
        pair = (first.get(key, F()), second.get(key, F()))
        state, value = calls[key]
        if (state == "MEASURED"):
            known = [known[index] + pair[index] * value for (index) in (range(2))]
        else:
            unknown.append((key, pair))
    scales = [
        math.lcm(known[index].denominator, *(pair[index].denominator for (key, pair) in (unknown)))
        for (index) in (range(2))
    ]
    integer_known = tuple(int(known[index] * scales[index]) for (index) in (range(2)))
    integer_terms = [
        (key, tuple(int(pair[index] * scales[index]) for (index) in (range(2))))
        for (key, pair) in (unknown)
    ]
    return scales, integer_known, integer_terms

def reachable(scales, known, terms, limit = 1000000, seconds = 30):
    started = time.monotonic()
    states = {known: 0}
    maximum = 1
    for (index, (key, pair)) in (enumerate(terms)):
        newer = states.copy()
        for (coordinates, witness) in (states.items()):
            target = (coordinates[0] + pair[0], coordinates[1] + pair[1])
            if (target not in newer):
                newer[target] = witness | (1 << index)
        states = newer
        maximum = max(maximum, len(states))
        if (maximum > limit or time.monotonic() - started > seconds):
            return {
                "status": "UNAVAILABLE_RESOURCE_LIMIT",
                "states": len(states),
                "terms_processed": index + 1,
                "unknown_terms": len(terms),
            }, None
    witnesses = {}
    for (point, mask) in (states.items()):
        category = "opposed" if (point[0] * point[1] < 0) else "not_opposed"
        if (category not in witnesses):
            assignment = [
                {"sample_id": key[0], "gene": key[1], "value": (mask >> index) & 1}
                for (index, (key, pair)) in (enumerate(terms))
            ]
            witnesses[category] = {
                "PRODIGY": F(point[0], scales[0]),
                "DawnRank": F(point[1], scales[1]),
                "assignment": assignment,
            }
    output = {
        "status": "EXACT_COMPLETE",
        "unknown_terms": len(terms),
        "distinct_reachable_pairs": len(states),
        "maximum_states": maximum,
        "minimum_opposition": int("not_opposed" not in witnesses),
        "maximum_opposition": int("opposed" in witnesses),
        "witnesses": witnesses,
        "seconds": time.monotonic() - started,
    }
    return output, set(states)

def analyse(models, allocations, means, degrees, calls, opposition = True):
    updates, donor_forms = {}, {}
    degree_update = {
        row["model_id"]: combine(
            [(1, degrees[row["model_id"]]["G1"]), (-1, degrees[row["model_id"]]["G0"])]
        )
        for (row) in (models)
    }
    updates["degree_G1_minus_G0"], _ = summary(models, degree_update, calls)
    for (method) in (METHODS):
        joint = {
            row["model_id"]: combine(
                [
                    (1, means[(method, row["model_id"], "G1H1")]),
                    (-1, means[(method, row["model_id"], "G0H0")]),
                ]
            )
            for (row) in (models)
        }
        updates[method + "_joint"], donor_forms[method] = summary(models, joint, calls)
        residual = {
            mid: combine([(1, vector), (-1, degree_update[mid])])
            for (mid, vector) in (joint.items())
        }
        updates[method + "_joint_minus_degree_update"], _ = summary(models, residual, calls)
    output = {
        "models": len(models),
        "donors": len({row["donor_id"] for (row) in (models)}),
        "budget": 10,
        "update_comparator": updates,
    }
    mean_directions = {
        method: {row["donor_id"]: row["sign"] for (row) in (updates[method + "_joint"]["donors"])}
        for (method) in (METHODS)
    }
    mean_opposed = sorted(
        donor
        for (donor) in (mean_directions["PRODIGY"])
        if (is_opposed(mean_directions["PRODIGY"][donor], mean_directions["DawnRank"][donor]))
    )
    seed_rows, seed_sets = [], {}
    for (seed) in (SEEDS):
        vectors = {
            row["model_id"]: combine(
                [
                    (1, allocations[("PRODIGY", row["model_id"], "G1H1", str(seed))]),
                    (-1, allocations[("PRODIGY", row["model_id"], "G0H0", str(seed))]),
                ]
            )
            for (row) in (models)
        }
        reduced, _ = summary(models, vectors, calls)
        selected = []
        for (row) in (reduced["donors"]):
            donor = row["donor_id"]
            dawn = next(
                record
                for (record) in (updates["DawnRank_joint"]["donors"])
                if (record["donor_id"] == donor)
            )
            opposed = is_opposed(row["sign"], dawn["sign"])
            if (opposed):
                selected.append(donor)
            seed_rows.append(
                {
                    "seed": seed,
                    "donor_id": donor,
                    "prodigy_lower": row["lower"],
                    "prodigy_upper": row["upper"],
                    "prodigy_sign": row["sign"],
                    "dawn_lower": dawn["lower"],
                    "dawn_upper": dawn["upper"],
                    "dawn_sign": dawn["sign"],
                    "opposed": opposed,
                    "seed_mean_opposed": donor in mean_opposed,
                }
            )
        seed_sets[str(seed)] = selected
    frequencies = {
        donor: sum(donor in selected for (selected) in (seed_sets.values()))
        for (donor) in (sorted(mean_directions["PRODIGY"]))
    }
    output["seed_opposition"] = {
        "seed_mean_opposed": mean_opposed,
        "sets": seed_sets,
        "intersection": sorted(set.intersection(*(set(value) for (value) in (seed_sets.values())))),
        "union": sorted(set.union(*(set(value) for (value) in (seed_sets.values())))),
        "frequency_by_donor": frequencies,
        "pairwise": [
            {
                "seed_a": first,
                "seed_b": second,
                "intersection": sorted(
                    set(seed_sets[str(first)]) & set(seed_sets[str(second)])
                ),
            }
            for (position, first) in (enumerate(SEEDS))
            for (second) in (SEEDS[position + 1 :])
        ],
        "all_donor_seed_rows": len(seed_rows),
        "frequency_histogram": {
            str(count): sum(value == count for (value) in (frequencies.values()))
            for (count) in (range(6))
        },
        "seed_mean_opposed_frequencies": {donor: frequencies[donor] for (donor) in (mean_opposed)},
    }
    output["seed_rows"] = seed_rows
    output["evidence_zero"], output["evidence_rows"] = evidence_diagnostics(
        models, means, calls
    )
    if (opposition):
        output["joint_opposition"] = opposition_bounds(donor_forms, calls, mean_opposed)
    return output

def evidence_diagnostics(models, means, calls):
    evidence_rows, evidence_summary = [], []
    for (method) in (METHODS):
        for (graph) in (("G0", "G1")):
            vectors = {
                row["model_id"]: combine(
                    [
                        (1, means[(method, row["model_id"], graph + "H1")]),
                        (-1, means[(method, row["model_id"], graph + "H0")]),
                    ]
                )
                for (row) in (models)
            }
            reduced, forms = summary(models, vectors, calls)
            for (unit) in (("models", "donors")):
                selected = []
                for (row) in (reduced[unit]):
                    vector = (
                        vectors[row["model_id"]]
                        if (unit == "models")
                        else forms[row["donor_id"]]
                    )
                    l1 = sum((abs(value) for (value) in (vector.values())), F())
                    record = {
                        "method": method,
                        "graph": graph,
                        "unit": unit[:-1],
                        **row,
                        "l1": l1,
                        "changed": bool(vector),
                        "known_zero": row["known"] == 0,
                        "zero_identified": row["lower"] == row["upper"] == 0,
                    }
                    evidence_rows.append(record)
                    selected.append(record)
                subgroup = [row for (row) in (selected) if (row["changed"] and row["known_zero"])]
                evidence_summary.append(
                    {
                        "method": method,
                        "graph": graph,
                        "unit": unit[:-1],
                        "total": len(selected),
                        "changed": sum(row["changed"] for (row) in (selected)),
                        "changed_known_zero": len(subgroup),
                        "changed_known_zero_identified_zero": sum(
                            row["zero_identified"] for (row) in (subgroup)
                        ),
                        "changed_known_zero_unresolved": sum(
                            row["sign"] == "unresolved" for (row) in (subgroup)
                        ),
                    }
                )
    return evidence_summary, evidence_rows

def opposition_bounds(donor_forms, calls, mean_opposed):
    optional, started = [], time.monotonic()
    for (donor) in (sorted(donor_forms["PRODIGY"])):
        if (time.monotonic() - started > 600):
            optional.append({"donor_id": donor, "status": "UNAVAILABLE_TOTAL_TIME_LIMIT"})
            continue
        scales, known, terms = make_lattice(
            donor_forms["PRODIGY"][donor], donor_forms["DawnRank"][donor], calls
        )
        result, _ = reachable(scales, known, terms)
        optional.append(
            {
                "donor_id": donor,
                "scales": scales,
                "integer_known": known,
                "integer_terms": [
                    {"sample_id": key[0], "gene": key[1], "pair": pair} for (key, pair) in (terms)
                ],
                **result,
            }
        )
    complete = all(row["status"] == "EXACT_COMPLETE" for (row) in (optional))
    return {
        "status": "EXACT_COMPLETE" if (complete) else "PARTIAL_NO_GLOBAL_SHARP_BOUND",
        "donors": optional,
        "minimum_opposing_donors": sum(row["minimum_opposition"] for (row) in (optional))
        if (complete)
        else None,
        "maximum_opposing_donors": sum(row["maximum_opposition"] for (row) in (optional))
        if (complete)
        else None,
        "separately_identified_opposing_donors": mean_opposed,
    }
