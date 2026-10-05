import math
from collections import defaultdict
from fractions import Fraction

CELLS = ("G0H0", "G1H0", "G0H1", "G1H1")

SEEDS = (104729, 130363, 155921, 196613, 228017)

BUDGETS = (1, 5, 10, 20)

METHODS = ("PRODIGY", "DawnRank")

CONTRASTS = {
    "graph_at_old_evidence": {"G1H0": 1, "G0H0": -1},
    "graph_at_new_evidence": {"G1H1": 1, "G0H1": -1},
    "evidence_at_old_graph": {"G0H1": 1, "G0H0": -1},
    "evidence_at_new_graph": {"G1H1": 1, "G1H0": -1},
    "interaction": {"G1H1": 1, "G1H0": -1, "G0H1": -1, "G0H0": 1},
    "joint": {"G1H1": 1, "G0H0": -1},
}

def require(condition, message):
    if (not condition):
        raise ValueError(message)

def allocate(scores, budget):
    require(type(budget) is int and budget > 0, "Positive integer budget required")
    groups = defaultdict(list)
    for (gene, score) in (scores.items()):
        require(math.isfinite(score), "Nonfinite rank score")
        groups[score].append(gene)
    weights, left = {}, Fraction(budget)
    for (score) in (sorted(groups, reverse = True)):
        if (not left):
            break
        genes = groups[score]
        amount = min(Fraction(len(genes)), left)
        weights.update((gene, amount / len(genes)) for (gene) in (genes))
        left -= amount
    require(
        sum(weights.values(), Fraction()) == min(budget, len(scores)),
        "Allocation capacity identity",
    )
    return weights

def mean_weights(rows):
    require(rows, "Empty averaging set")
    output = defaultdict(Fraction)
    for (row) in (rows):
        for (gene, value) in (row.items()):
            output[gene] += value / len(rows)
    return dict(output)

def sign(lower, upper):
    return (
        "negative"
        if (upper < 0)
        else "positive"
        if (lower > 0)
        else "zero"
        if (lower == upper == 0)
        else "unresolved"
    )

def bound(coefficients, calls):
    known, lower, upper = Fraction(), Fraction(), Fraction()
    for (key, value) in (coefficients.items()):
        state, outcome = calls[key]
        if (state == "MEASURED"):
            known += value * outcome
        else:
            lower += min(value, 0)
            upper += max(value, 0)
    lower, upper = lower + known, upper + known
    return {
        "known": known,
        "lower": lower,
        "upper": upper,
        "width": upper - lower,
        "sign": sign(lower, upper),
    }

def distance(first, second):
    genes = set(first) | set(second)
    values = {
        gene: second.get(gene, Fraction()) - first.get(gene, Fraction()) for (gene) in (genes)
    }
    addition = sum((max(value, 0) for (value) in (values.values())), Fraction())
    removal = sum((max(-value, 0) for (value) in (values.values())), Fraction())
    return {
        "l1": addition + removal,
        "addition": addition,
        "removal": removal,
        "changed_genes": sum(value != 0 for (value) in (values.values())),
    }
