import csv
import gzip
import json
from pathlib import Path

from ..analysis.run import read_input
from . import workflow
from .review_rules import to_json

OPERATIONS = {
    "reference_join": workflow.reference_join,
    "cohort_features": workflow.cohort_features,
    "development": workflow.development,
    "frozen_choices": workflow.frozen_choices,
    "protected": workflow.protected,
    "selection_optimism": workflow.selection_optimism,
}

def run(config_path):
    config_path = Path(config_path).resolve()
    base = config_path.parent
    config = read_input(base, config_path.name)
    if (config["operation"] not in OPERATIONS):
        raise ValueError("Unknown catalogue operation")
    output = (base / config["output"]).resolve()
    if (output.exists()):
        raise FileExistsError(
            "Use a distinct output path; existing results are retained"
        )
    arguments = {}
    for (name, specification) in (config.get("inputs", {}).items()):
        if (isinstance(specification, dict) and specification.get("format") == "tsv"):
            path = (base / specification["path"]).resolve()
            opener = gzip.open if (path.suffix == ".gz") else open
            with opener(path, "rt", encoding = "utf-8-sig", newline = "") as stream:
                arguments[name] = list(csv.DictReader(stream, delimiter = "\t"))
        else:
            arguments[name] = read_input(base, specification)
        if (isinstance(specification, dict) and "select" in specification):
            if (not isinstance(specification["select"], list) or any(
                not isinstance(key, str) for (key) in (specification["select"])
            )):
                raise ValueError(
                    "A nested input selection must list explicit object keys"
                )
            for (key) in (specification["select"]):
                arguments[name] = arguments[name][key]
    options = config.get("options", {})
    if (set(arguments) & set(options)):
        raise ValueError("An argument cannot be both an input and an option")
    arguments.update(options)
    result = OPERATIONS[config["operation"]](**arguments)
    output.parent.mkdir(parents = True, exist_ok = True)
    with output.open("x", encoding = "utf-8", newline = "\n") as stream:
        json.dump(
            to_json(result), stream, ensure_ascii = False, allow_nan = False, indent = 2
        )
        stream.write("\n")
    return {"operation": config["operation"], "output": str(output)}
