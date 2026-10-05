from collections import defaultdict
from .contextual_order import allocation, compatibility

def build(ranks, eligible, models):
    by_model = {m["sample_ID"]: m for (m) in (models)}
    if (len(by_model) != 85 or len({m["study_donor_id"] for (m) in (models)}) != 83):
        raise ValueError("Incomplete original population")
    if (any(
        not m["study_donor_id"] or m["split"] not in ("development", "held_out")
        for (m) in (models)
    )):
        raise ValueError("Missing donor or invalid split")
    for (donor) in ({m["study_donor_id"] for (m) in (models)}):
        if (len({m["split"] for (m) in (models) if (m["study_donor_id"] == donor)}) != 1):
            raise ValueError("Donor split leakage")
    for (name, n, d) in ((("development", 26, 25), ("held_out", 59, 58))):
        selected = [m for (m) in (models) if (m["split"] == name)]
        if (len(selected) != n or len({m["study_donor_id"] for (m) in (selected)}) != d):
            raise ValueError("Invalid original subset")
    budgets = (1, 5, 10, 20)
    groups = defaultdict(dict)
    for (row) in (ranks):
        method, pool, release, seed, model = row["key"]
        if (pool != "all34"):
            continue
        key = (method, release, seed)
        if (model in groups[key]):
            raise ValueError("Duplicate ranking")
        menu = set(
            eligible["pools"][pool][model]["native_eligible_literal_gene_labels"]
        )
        groups[key][model] = {
            k: allocation(row["genes"], row["scores"], menu, k) for (k) in (budgets)
        }
    axes = {
        (m, r, s)
        for (m) in (("DawnRank", "PersonaDrive", "PRODIGY"))
        for (r) in (("native_11_0", "native_11_5", "native_12_0"))
        for (s) in (
            (104729, 130363, 155921, 196613, 228017) if (m == "PRODIGY") else (None,)
        )
    }
    if (set(groups) != axes or any(set(x) != set(by_model) for (x) in (groups.values()))):
        raise ValueError("Incomplete methods, versions, seeds or models")
    for (key, rows) in (sorted(groups.items(), key = lambda x: str(x[0]))):
        for (population) in (("all85", "development", "held_out")):
            chosen = {
                m: w
                for ((m, w)) in (rows.items())
                if (population == "all85" or by_model[m]["split"] == population)
            }
            for (budget) in ((1, 5, 20, "joint_1_5_10_20")):
                observations = {
                    m + "|" + str(k): w[k]
                    for ((m, w)) in (chosen.items())
                    for (k) in (budgets if (isinstance(budget, str)) else (budget,))
                }
                yield key, population, budget, chosen, observations

def analyse(ranks, eligible, models):
    output = []
    for (key, population, budget, chosen, observations) in (build(ranks, eligible, models)):
        result, unused = compatibility(observations)
        result.update(
            method = key[0],
            release = key[1],
            seed = key[2],
            population = population,
            budget = budget,
            model_count = len(chosen),
            allocation_count = len(observations),
            available_allocations = sum(w is not None for (w) in (observations.values())),
            all_planned_observations_available = all(
                w is not None for (w) in (observations.values())
            ),
            unavailable_allocations = [
                m for ((m, w)) in (observations.items()) if (w is None)
            ],
            empty_allocations = [m for ((m, w)) in (observations.items()) if (w == {})],
            all_returned_allocations = [
                m
                for ((m, w)) in (observations.items())
                if (w and all(v == 1 for (v) in (w.values())))
            ],
        )
        output.append(result)
    if (len(output) != 252):
        raise ValueError("Missing fixed analysis group")
    return {"groups": output}
