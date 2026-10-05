import math
from collections import Counter, defaultdict
from fractions import Fraction as F
from . import influence as D

CELLS = ("G0H0", "G1H0", "G0H1", "G1H1")

METHODS = ("PRODIGY", "DawnRank")

ALL_METHODS = METHODS + ("DawnRank_strict",)

def need(ok, message):
    if (not ok):
        raise ValueError(message)

def fweights(items):
    return {gene: F(weight) for (gene, weight) in (items) if (F(weight))}

def mean_weights(items):
    result = defaultdict(F)
    for (item) in (items):
        for (gene, weight) in (item.items()):
            result[gene] += weight / len(items)
    return {key: value for (key, value) in (result.items()) if (value)}

def diff(first, second):
    return {
        gene: first.get(gene, F()) - second.get(gene, F())
        for (gene) in (sorted(set(first) | set(second)))
        if (first.get(gene, F()) != second.get(gene, F()))
    }

def add(first, second):
    result = defaultdict(F, first)
    for (key, value) in (second.items()):
        result[key] += value
    return {key: value for (key, value) in (result.items()) if (value)}

def category_gene(gene, sets, roles):
    org, cell = gene in sets["organoid"], gene in sets["cell_line"]
    core = (
        "BOTH_CORE"
        if (org and cell)
        else "ORGANOID_ONLY"
        if (org)
        else "CELL_LINE_ONLY"
        if (cell)
        else "NOT_LISTED_IN_EITHER"
    )
    row = roles.get(gene)
    if (row is None):
        role = "UNLISTED"
    else:
        onco, tsg = "oncogene" in row["cosmic_tokens"], "TSG" in row["cosmic_tokens"]
        role = (
            "DUAL_ONCOGENE_TSG"
            if (onco and tsg)
            else "ONCOGENE"
            if (onco)
            else "TSG"
            if (tsg)
            else "OTHER_OR_UNSPECIFIED"
        )
    return {
        "core": core,
        "role": role,
        "auc_common_essential": gene in sets["auc_common"],
        "organoid_exclusive_source_flag": gene in sets["organoid_exclusive"],
        "source_role": row,
    }

def summarise_contrast(models, local, calls):
    population, donorcoeff = {}, defaultdict(dict)
    permodel = {}
    for (row) in (models):
        coeff = {
            (row["sample_id"], gene): value
            for (gene, value) in (local[row["model_id"]].items())
            if (value)
        }
        permodel[row["model_id"]] = D.bounds(coeff, calls)
        population.update({key: row["alpha"] * value for (key, value) in (coeff.items())})
        donorcoeff[row["donor_id"]].update(
            {key: row["alpha"] * 83 * value for (key, value) in (coeff.items())}
        )
    donors = {donor: D.bounds(coeff, calls) for (donor, coeff) in (sorted(donorcoeff.items()))}
    aggregate = D.bounds(population, calls)
    for (key) in (("known", "unknown_negative", "unknown_positive", "lower", "upper", "width")):
        need(
            aggregate[key] == sum((value[key] for (value) in (donors.values())), F()) / 83,
            "Donor aggregation differs",
        )
    return {
        "population": aggregate,
        "models": permodel,
        "donors": donors,
        "donor_sign_counts": dict(Counter(row["sign"] for (row) in (donors.values()))),
        "nonzero_identities": len(population),
    }, population

def level_summary(models, allocations, calls):
    values = {
        row["model_id"]: D.mass(allocations[row["model_id"]], row["sample_id"], calls)[0]
        for (row) in (models)
    }
    donors = {
        donor: D.weighted(
            [
                (row["alpha"] * 83, values[row["model_id"]])
                for (row) in (models)
                if (row["donor_id"] == donor)
            ]
        )
        for (donor) in (sorted({row["donor_id"] for (row) in (models)}))
    }
    return {
        "models": values,
        "donors": donors,
        "population": D.weighted([(F(1, 83), value) for (value) in (donors.values())]),
    }

def ranks(values):
    groups = defaultdict(list)
    for (index, value) in (enumerate(values)):
        groups[value].append(index)
    result, count = [None] * len(values), 0
    for (value) in (sorted(groups)):
        group = groups[value]
        rank = F(2 * count + len(group) + 1, 2)
        for (index) in (group):
            result[index] = rank
        count += len(group)
    return result

def spearman(first, second):
    if (len(first) < 2):
        return None
    x, y = ranks(first), ranks(second)
    mx, my = sum(x, F()) / len(x), sum(y, F()) / len(y)
    numerator = sum(((a - mx) * (b - my) for (a, b) in (zip(x, y))), F())
    sx, sy = sum(((a - mx) ** 2 for (a) in (x)), F()), sum(((b - my) ** 2 for (b) in (y)), F())
    return float(numerator) / math.sqrt(float(sx * sy)) if (sx and sy) else None

def quantile(values, p):
    values = sorted(values)
    if (not values):
        return None
    index = (len(values) - 1) * p
    lower = int(index)
    return values[lower] + (values[min(lower + 1, len(values) - 1)] - values[lower]) * (
        index - lower
    )
