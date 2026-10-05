from collections import defaultdict
from fractions import Fraction
from itertools import product

from ..precision.io import require

def state_name(state):
    a, b, h = state
    return f"P{a}_D{b}_H{h}"

def state_names():
    return tuple(state_name(state) for (state) in (product((0, 1), repeat = 3)))

def contrasts():
    result = {}
    for (axis, label) in (enumerate(("P", "D", "H"))):
        other = [position for (position) in (range(3)) if (position != axis)]
        for (values) in (product((0, 1), repeat = 2)):
            low = [0, 0, 0]
            for (position, value) in (zip(other, values)):
                low[position] = value
            high = low.copy()
            high[axis] = 1
            suffix = "_".join(
                f"{('P', 'D', 'H')[position]}{value}" for (position, value) in (zip(other, values))
            )
            result[f"{label}_{suffix}"] = {state_name(high): 1, state_name(low): -1}
    for (axes, label) in ((((0, 1), "PD"), ((0, 2), "PH"), ((1, 2), "DH"))):
        fixed = next(position for (position) in (range(3)) if (position not in axes))
        for (value) in ((0, 1)):
            multipliers = {}
            for (values) in (product((0, 1), repeat = 2)):
                state = [0, 0, 0]
                state[fixed] = value
                for (position, bit) in (zip(axes, values)):
                    state[position] = bit
                multipliers[state_name(state)] = (-1) ** (2 - sum(values))
            result[f"{label}_{('P', 'D', 'H')[fixed]}{value}"] = multipliers
    result["PDH"] = {
        state_name(state): (-1) ** (3 - sum(state)) for (state) in (product((0, 1), repeat = 3))
    }
    for (h) in ((0, 1)):
        result[f"GRAPH_H{h}"] = {state_name((1, 1, h)): 1, state_name((0, 0, h)): -1}
    result["JOINT"] = {state_name((1, 1, 1)): 1, state_name((0, 0, 0)): -1}
    require(len(result) == 22, "Incomplete component contrast family")
    return result

def categories():
    return ("ALL", "BOTH_CORE", "ORGANOID_ONLY", "CELL_ONLY", "NEITHER_CORE", "NOT_BOTH_CORE")

def core_sets(document):
    organoid = set(document["sets"]["organoid"])
    cell = set(document["sets"]["cell_line"])
    require(
        len(organoid) == 751 and len(cell) == 1121, "Source-fixed core category sizes differ"
    )
    return organoid, cell

def category_map(genes, document):
    organoid, cell = core_sets(document)
    result = {}
    for (gene) in (genes):
        if (gene in organoid and gene in cell):
            result[gene] = "BOTH_CORE"
        elif (gene in organoid):
            result[gene] = "ORGANOID_ONLY"
        elif (gene in cell):
            result[gene] = "CELL_ONLY"
        else:
            result[gene] = "NEITHER_CORE"
    return result

def belongs(gene, selected, mapping):
    return (
        selected == "ALL"
        or (selected == "NOT_BOTH_CORE" and mapping[gene] != "BOTH_CORE")
        or selected == mapping[gene]
    )

def signed_coefficients(allocations, multipliers, samples, alpha):
    result = defaultdict(Fraction)
    for (sample) in (samples):
        for (state, multiplier) in (multipliers.items()):
            for (gene, amount) in (allocations[(state, sample)].items()):
                result[(sample, gene)] += alpha[sample] * multiplier * amount
    return {key: value for (key, value) in (result.items()) if (value != 0)}

def interval(coefficients, labels):
    known, lower, upper = Fraction(0), Fraction(0), Fraction(0)
    for (key, coefficient) in (coefficients.items()):
        state, value = labels[key]
        if (state == "MEASURED"):
            known += coefficient * value
        else:
            lower += min(coefficient, Fraction(0))
            upper += max(coefficient, Fraction(0))
    lower += known
    upper += known
    sign = (
        "NEGATIVE"
        if (upper < 0)
        else "POSITIVE"
        if (lower > 0)
        else "ZERO"
        if (lower == upper == 0)
        else "UNRESOLVED"
    )
    return {
        "R": sum(coefficients.values(), Fraction(0)),
        "H": known,
        "lower": lower,
        "upper": upper,
        "sign": sign,
    }
