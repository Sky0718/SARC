from collections import Counter
from fractions import Fraction

from .allocation import BUDGETS, CELLS, SEEDS, require
from .influence import labels_from

def fraction(value):
    if (isinstance(value, dict)):
        require(set(value) == {"numerator", "denominator"}, "Invalid rational object")
        return Fraction(value["numerator"], value["denominator"])
    require(not isinstance(value, (float, bool)), "Use exact rational strings")
    return Fraction(value)

def weights(value, budget = 10):
    rows = list(value.items()) if (isinstance(value, dict)) else value
    result = {gene: fraction(amount) for (gene, amount) in (rows)}
    require(len(rows) == len(result), "Duplicate allocated gene")
    require(all(0 <= amount <= 1 for (amount) in (result.values())), "Invalid allocation")
    require(sum(result.values(), Fraction()) <= budget, "Allocation exceeds budget")
    return {gene: amount for (gene, amount) in (result.items()) if (amount)}

def roster(document, expected_models = 85, expected_donors = 83):
    rows = document["models"] if (isinstance(document, dict)) else document
    require(len(rows) == expected_models, "Incomplete model roster")
    require(len({row["model_id"] for (row) in (rows)}) == len(rows), "Duplicate model")
    require(len({row["sample_id"] for (row) in (rows)}) == len(rows), "Duplicate sample")
    counts = Counter(str(row["donor_id"]) for (row) in (rows))
    require(len(counts) == expected_donors, "Incomplete donor roster")
    result = []
    for (row) in (rows):
        model = dict(row)
        model["donor_id"] = str(model["donor_id"])
        model["alpha"] = Fraction(1, expected_donors * counts[model["donor_id"]])
        if ("alpha" in row):
            require(fraction(row["alpha"]) == model["alpha"], "Donor weight mismatch")
        if ("methods" in row):
            model["methods"] = {
                method: {cell.replace("_", ""): weights(value) for (cell, value) in (arms.items())}
                for (method, arms) in (row["methods"].items())
            }
            require(
                all(set(arms) == set(CELLS) for (arms) in (model["methods"].values())),
                "Every method requires four complete states",
            )
        result.append(model)
    return result

def calls(documents):
    if (isinstance(documents, dict)):
        documents = [documents]
    return labels_from([(str(index), item) for (index, item) in (enumerate(documents))])

def rank_allocations(models, documents):
    model_ids = {row["model_id"] for (row) in (models)}
    result = {}
    for (document) in (documents):
        method = document["method"]
        require(method in ("PRODIGY", "DawnRank"), "Unknown ranking method")
        for (row) in (document["allocations"]):
            seed = row.get("seed", row.get("master_seed", "deterministic"))
            if (method == "DawnRank"):
                seed = "deterministic"
            cell = row["cell"].replace("_", "")
            budget = row["budget"]
            if (method == "DawnRank" and budget == 10):
                require(
                    row.get("numerically_certified") is True,
                    "Primary DawnRank cutoff must be certified",
                )
            require(row["model_id"] in model_ids, "Unknown allocation model")
            require(cell in CELLS and budget in BUDGETS, "Unknown state or budget")
            require(
                seed in SEEDS if (method == "PRODIGY") else seed == "deterministic",
                "Algorithm seed differs",
            )
            key = (method, row["model_id"], cell, seed, budget)
            value = weights(row["weights"], budget)
            if (key in result):
                require(method == "DawnRank" and result[key] == value, "Duplicate state")
            result[key] = value
    expected = {
        (method, model, cell, seed, budget)
        for (method) in (("PRODIGY", "DawnRank"))
        for (model) in (model_ids)
        for (cell) in (CELLS)
        for (seed) in (SEEDS if (method == "PRODIGY") else ("deterministic",))
        for (budget) in (BUDGETS)
    }
    require(set(result) == expected, "Incomplete method/model/state/seed/budget matrix")
    return result

def cover_labels(models, allocations, selected):
    samples = {row["model_id"]: row["sample_id"] for (row) in (models)}
    requested = {
        (samples[model], gene)
        for ((method, model, cell, seed, budget), values) in (allocations.items())
        for (gene) in (values)
    }
    require(requested <= set(selected), "Every selected identity needs an explicit label state")
