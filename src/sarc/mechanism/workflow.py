from pathlib import Path

from ..precision.io import read_json, require, resolve
from . import functional, score_reference

def run(config_path):
    config_path = Path(config_path).resolve()
    config = read_json(config_path)
    operation = config["operation"]
    require(
        operation in ("freeze", "evaluate", "score_reference"),
        "Unknown component analysis operation",
    )
    output = resolve(config_path.parent, config["output"])
    output.mkdir(parents = True, exist_ok = False)
    if (operation == "freeze"):
        document = read_json(resolve(config_path.parent, config["allocations"]))
        return functional.freeze(document, output)
    categories = read_json(resolve(config_path.parent, config["categories"]))
    if (operation == "evaluate"):
        frozen = resolve(config_path.parent, config["frozen_directory"])
        document = read_json(frozen / "allocations.json")
        requested = read_json(frozen / "label_requests.json")["requests"]
        labels = read_json(resolve(config_path.parent, config["labels"]))
        require(
            {(row["sample_id"], row["gene"]) for (row) in (requested)}
            == {(row["sample_id"], row["gene"]) for (row) in (labels["calls"])},
            "Frozen selected-label request support differs",
        )
        return functional.evaluate(document, labels, categories, output)
    return score_reference.analyse(
        resolve(config_path.parent, config["scores_directory"]), categories, output
    )
