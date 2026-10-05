from .allocation import CELLS, SEEDS, mean_weights, require
from .inputs import fraction, roster, weights

def native_menus(pro, dawn):
    models = roster(pro)
    dawn_models = {row["model_id"]: row for (row) in (roster(dawn))}
    require(set(dawn_models) == {row["model_id"] for (row) in (models)}, "Method rosters differ")
    result = list(pro["native"])
    for (row) in (dawn["allocations"]):
        if (row["budget"] != 10):
            continue
        model = dawn_models[row["model_id"]]
        require(row["numerically_certified"] is True, "Uncertified DawnRank allocation")
        result.append({**row, "method": "DawnRank", "menu": model["query"]})
    return result

def records(pro, dawn, certificates):
    models = roster(pro)
    indexed = {}
    for (row) in (native_menus(pro, dawn)):
        key = (row["method"], row["model_id"], row["cell"].replace("_", ""), row["seed"])
        require(key not in indexed, "Duplicate native menu")
        indexed[key] = row
    cuts = {}
    for (row) in (certificates):
        if (int(row["budget"]) != 10):
            continue
        key = (row["sample_id"], row["state"].replace("_", ""))
        require(key not in cuts, "Duplicate cutoff record")
        cuts[key] = row
    output = []
    consumed = set()
    for (model) in (models):
        query = model.get("query", model.get("eligible_genes"))
        current = {
            "model": {**model, "eligible_genes": query},
            "pro": {},
            "dawn": {},
            "cutoff": [],
        }
        for (cell) in (CELLS):
            label = cell[:2] + "_" + cell[2:]
            for (method, name, seeds) in ((
                ("PRODIGY", "pro", SEEDS),
                ("DawnRank", "dawn", ("deterministic",)),
            )):
                rows = []
                for (seed) in (seeds):
                    key = (method, model["model_id"], cell, seed)
                    require(key in indexed, "Incomplete native matrix")
                    consumed.add(key)
                    rows.append(indexed[key])
                menus = [set(row["menu"]) for (row) in (rows)]
                require(all(menu <= set(query) for (menu) in (menus)), "Menu outside query")
                current[name][label] = {
                    "weights": mean_weights([weights(row["weights"]) for (row) in (rows)]),
                    "menu": sorted(set().union(*menus)),
                    "seed_menu_counts": [len(menu) for (menu) in (menus)],
                    "seed_menus_identical": all(menu == menus[0] for (menu) in (menus)),
                    "statuses": [row.get("native_status", "SUCCESS") for (row) in (rows)],
                }
            cut = cuts[(model["sample_id"], cell)]
            gap = indexed[("DawnRank", model["model_id"], cell, "deterministic")][
                "cutoff_gap_exact"
            ]
            current["cutoff"].append(
                {
                    "cell": label,
                    "cutoff_gap": None if (gap is None) else fraction(gap),
                    "boundary_tie": gap is not None and fraction(gap) == 0,
                    "iterations": int(cut["iterations"]),
                    "cutoff_certified": cut["certified"] == "TRUE",
                }
            )
        output.append(current)
    require(consumed == set(indexed), "Extra native states")
    require(
        set(cuts) == {(model["sample_id"], cell) for (model) in (models) for (cell) in (CELLS)},
        "Incomplete or extra cutoff states",
    )
    return output
