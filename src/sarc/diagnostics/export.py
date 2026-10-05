import struct
from fractions import Fraction

from ..analysis import prodigy_factorial
from .allocation import BUDGETS, allocate, require

def serial(value):
    if (isinstance(value, Fraction)):
        return str(value)
    if (isinstance(value, dict)):
        return {key: serial(item) for (key, item) in (value.items())}
    if (isinstance(value, (list, tuple))):
        return [serial(item) for (item) in (value)]
    if (isinstance(value, set)):
        return [serial(item) for (item) in (sorted(value))]
    return value

def prodigy_allocations(scope, carriers):
    frozen = prodigy_factorial.freeze(scope, carriers)
    require(frozen["scientifically_complete"], "Complete PRODIGY states required")
    by_model = {row["model_id"]: row for (row) in (scope["models"])}
    records, menus = [], []
    for (cell) in (carriers["cells"]):
        model = by_model[cell["model_id"]]
        query = set(model["eligible_genes"])
        name = "G" + str(scope["graph_states"].index(cell["graph"])) + cell["evidence"]
        for (state) in (cell["seed_states"]):
            checked = prodigy_factorial.allocate(state, model)
            require(checked is not None, "Unavailable PRODIGY state")
            scores = {
                gene: struct.unpack("<d", bytes.fromhex(code))[0]
                for (gene, code) in (zip(state["genes"], state["scores_binary64_le"]))
                if (gene in query)
            }
            for (budget) in (BUDGETS):
                weights = allocate(scores, budget)
                if (budget == 10):
                    require(
                        {gene: value for (gene, value) in (checked.items()) if (value)} == weights,
                        "Native allocation convention differs",
                    )
                row = {
                    "method": "PRODIGY",
                    "model_id": model["model_id"],
                    "sample_id": model["sample_id"],
                    "donor_id": model["donor_id"],
                    "cell": name,
                    "seed": state["master_seed"],
                    "budget": budget,
                    "status": "AVAILABLE",
                    "native_status": state["status"],
                    "capacity": min(budget, len(scores)),
                    "menu": sorted(scores),
                    "weights": weights,
                }
                records.append(row)
                if (budget == 10):
                    menus.append({**row, "menu": sorted(scores)})
    models = [
        {key: model[key] for (key) in (("model_id", "sample_id", "donor_id"))}
        | {"query": model["eligible_genes"], "eligible_genes": model["eligible_genes"]}
        for (model) in (scope["models"])
    ]
    return {"method": "PRODIGY", "models": models, "allocations": records, "native": menus}

def packaged(operation):
    def invoke(**arguments):
        return serial(operation(**arguments))

    return invoke
