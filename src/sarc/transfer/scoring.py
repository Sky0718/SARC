import itertools
import math
import random
from fractions import Fraction

def score_list(genes, scores, eligible, positives, k = 10):
    eligible = frozenset(eligible)
    positive_ids = frozenset(positives) & eligible
    if (not isinstance(k, int) or isinstance(k, bool) or k < 1):
        raise ValueError("Invalid review scale")
    output = {
        "eligible_count": len(eligible),
        "positive_count": len(positive_ids),
        "k": k,
        "lower": Fraction(0),
        "upper": Fraction(min(k, len(positive_ids))),
    }
    if (genes is None):
        if (scores is not None):
            raise ValueError("Failed list cannot expose scores")
        return output | {
            "status": "UNAVAILABLE",
            "hits": None,
            "recall": None,
            "native_output_count": None,
            "common_output_count": None,
            "positive_output_count": None,
        }
    if (len(genes) != len(set(genes)) or any(
        (not isinstance(gene, str) or not gene for (gene) in (genes))
    )):
        raise ValueError("Invalid ranked identities")
    if (scores is not None and (
        len(scores) != len(genes)
        or any((not math.isfinite(value) for (value) in (scores)))
        or any((a < b for ((a, b)) in (zip(scores, scores[1:]))))
    )):
        raise ValueError("Invalid ranked scores")
    retained = [
        (gene, position if (scores is None) else scores[position])
        for ((position, gene)) in (enumerate(genes))
        if (gene in eligible)
    ]
    hits = Fraction(0)
    used = 0
    for (value, entries) in (itertools.groupby(retained, key = lambda item: item[1])):
        group = list(entries)
        slots = min(len(group), max(0, k - used))
        hits += Fraction(
            slots * sum((gene in positive_ids for ((gene, unused)) in (group))),
            len(group),
        )
        used += len(group)
        if (used >= k):
            break
    return output | {
        "status": "SUCCESS" if (positive_ids) else "SUCCESS_NO_REFERENCE_POSITIVES",
        "hits": hits,
        "recall": hits / len(positive_ids) if (positive_ids) else None,
        "lower": hits,
        "upper": hits,
        "native_output_count": len(genes),
        "common_output_count": len(retained),
        "positive_output_count": sum(
            (gene in positive_ids for ((gene, value)) in (retained))
        ),
    }

def mean_records(records):
    if (not records):
        raise ValueError("Empty planned population")
    n = len(records)
    complete = [row for (row) in (records) if (row["hits"] is not None)]
    return {
        "population_count": n,
        "complete_count": len(complete),
        "value": sum((row["hits"] for (row) in (complete))) / n
        if (len(complete) == n)
        else None,
        "lower": sum((row["lower"] for (row) in (records))) / n,
        "upper": sum((row["upper"] for (row) in (records))) / n,
        "complete_subset_value": sum((row["hits"] for (row) in (complete)))
        / len(complete)
        if (complete)
        else None,
    }

def repeat_summary(records):
    if (len(records) != 5):
        raise ValueError("All five planned stochastic repeats are required")
    result = mean_records(records)
    result["hits"] = result.pop("value")
    available = [row["hits"] for (row) in (records) if (row["hits"] is not None)]
    result["repeat_minimum"] = (
        min(available) if (len(available) == len(records)) else None
    )
    result["repeat_maximum"] = (
        max(available) if (len(available) == len(records)) else None
    )
    return result

def paired_contrast(population, terms):
    if (
        not population
        or len(population) != len(set(population))
        or any((set(rows) != set(population) for ((coefficient, rows)) in (terms)))
    ):
        raise ValueError("Contrast population mismatch")
    results = []
    for (sample) in (population):
        complete = all(
            (rows[sample]["hits"] is not None for ((coefficient, rows)) in (terms))
        )
        hits = (
            sum((coefficient * rows[sample]["hits"] for ((coefficient, rows)) in (terms)))
            if (complete)
            else None
        )
        lower = sum(
            (
                coefficient * rows[sample]["lower" if (coefficient > 0) else "upper"]
                for ((coefficient, rows)) in (terms)
            )
        )
        upper = sum(
            (
                coefficient * rows[sample]["upper" if (coefficient > 0) else "lower"]
                for ((coefficient, rows)) in (terms)
            )
        )
        results.append(
            {"sample_id": sample, "hits": hits, "lower": lower, "upper": upper}
        )
    return mean_records(results) | {"sample_contrasts": results}

def bootstrap_indices(n, draws = 2000, seed = 20260926):
    if (n < 1 or draws < 2):
        raise ValueError("Invalid bootstrap geometry")
    rng = random.Random(seed)
    return [[rng.randrange(n) for (unused) in (range(n))] for (draw) in (range(draws))]

def percentile(sorted_values, probability):
    index = (len(sorted_values) - 1) * probability
    low = math.floor(index)
    high = math.ceil(index)
    return sorted_values[low] + (index - low) * (
        sorted_values[high] - sorted_values[low]
    )

def paired_interval(result, indices):
    values = [row["hits"] for (row) in (result["sample_contrasts"])]
    if (any((value is None for (value) in (values)))):
        return {
            "lower_95": None,
            "upper_95": None,
            "status": "UNAVAILABLE_FULL_POPULATION_CONTRAST",
        }
    n = len(values)
    vector = [float(value) for (value) in (values)]
    if (any(
        (
            len(draw) != n
            or any((position < 0 or position >= n for (position) in (draw)))
            for (draw) in (indices)
        )
    )):
        raise ValueError("Bootstrap population mismatch")
    means = sorted(
        (
            math.fsum((vector[position] for (position) in (draw))) / n
            for (draw) in (indices)
        )
    )
    return {
        "lower_95": percentile(means, 0.025),
        "upper_95": percentile(means, 0.975),
        "status": "CONDITIONAL_DESCRIPTIVE_PERCENTILE",
        "draws": len(indices),
    }

def reversals(old_gap, new_gap, old_interval, new_interval, delta = 1):
    if (old_gap is None or new_gap is None):
        return {
            "nominal": None,
            "meaningful": None,
            "uncertainty_supported": None,
            "delta": delta,
        }
    nominal = old_gap * new_gap < 0
    meaningful = (
        old_gap > delta and new_gap < -delta or (old_gap < -delta and new_gap > delta)
    )
    intervals_complete = all(
        (
            item.get(key) is not None
            for (item) in ((old_interval, new_interval))
            for (key) in (("lower_95", "upper_95"))
        )
    )
    supported = (
        meaningful
        and (
            old_interval["lower_95"] > delta
            and new_interval["upper_95"] < -delta
            or (old_interval["upper_95"] < -delta and new_interval["lower_95"] > delta)
        )
        if (intervals_complete)
        else None
    )
    return {
        "nominal": bool(nominal),
        "meaningful": bool(meaningful),
        "uncertainty_supported": supported,
        "delta": delta,
    }

def leading_membership(genes, scores, eligible, k):
    if (genes is None):
        return None
    score_list(genes, scores, eligible, set(), k)
    rows = [
        (gene, index if (scores is None) else scores[index])
        for ((index, gene)) in (enumerate(genes))
        if (gene in eligible)
    ]
    membership = {}
    used = 0
    for (value, entries) in (itertools.groupby(rows, key = lambda item: item[1])):
        entries = list(entries)
        slots = min(len(entries), max(0, k - used))
        if (slots):
            membership.update(
                {gene: Fraction(slots, len(entries)) for ((gene, score)) in (entries)}
            )
        used += len(entries)
        if (used >= k):
            break
    return membership

def leading_stability(old, new):
    if (old is None or new is None):
        return {
            "status": "UNAVAILABLE",
            "jaccard": None,
            "fractional_membership_jaccard": None,
        }
    union = old.keys() | new.keys()
    if (not union):
        return {
            "status": "BOTH_SUCCESSFUL_EMPTY",
            "jaccard": None,
            "fractional_membership_jaccard": None,
        }
    weighted_union = sum(
        (max(old.get(gene, 0), new.get(gene, 0)) for (gene) in (union))
    )
    return {
        "status": "SUCCESS",
        "jaccard": Fraction(len(old.keys() & new.keys()), len(union)),
        "fractional_membership_jaccard": sum(
            (min(old.get(gene, 0), new.get(gene, 0)) for (gene) in (union))
        )
        / weighted_union,
    }

def json_number(value):
    if (isinstance(value, Fraction)):
        return {
            "numerator": value.numerator,
            "denominator": value.denominator,
            "value": float(value),
        }
    if (isinstance(value, dict)):
        return {key: json_number(child) for ((key, child)) in (value.items())}
    if (isinstance(value, (list, tuple))):
        return [json_number(child) for (child) in (value)]
    return value
