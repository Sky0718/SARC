from fractions import Fraction
from pathlib import Path

from ..resources.io import check, read_json, write_json
from . import panels
from . import intervals
from . import policies

def restore_exact(value):
    if (isinstance(value, list)):
        return [restore_exact(item) for (item) in (value)]
    if (isinstance(value, dict)):
        if (
            set(value) == {"numerator", "denominator"}
            and value["denominator"] is not None
        ):
            return Fraction(value["numerator"], value["denominator"])
        return {key: restore_exact(item) for ((key, item)) in (value.items())}
    return value

def load_config(config_path):
    return read_json(config_path), Path(config_path).resolve().parent

def load_input(root, value):
    path = Path(value)
    return restore_exact(read_json(path if (path.is_absolute()) else root / path))

def candidate_measurements(value):
    check(
        isinstance(value, dict),
        "Candidate measurements require a method::version keyed object",
    )
    result = {}
    for (key, rows) in (value.items()):
        parts = key.split("::") if (isinstance(key, str)) else ()
        check(
            len(parts) == 2 and all(parts),
            "Candidate measurement keys must be method::version",
        )
        candidate = tuple(parts)
        check(candidate not in result, "Duplicate candidate measurements")
        result[candidate] = rows
    return result

def evaluate_transfer(config_path):
    config, root = load_config(config_path)
    ranks = load_input(root, config["ranks"])
    support = load_input(root, config["support"])
    references = load_input(root, config["references"])
    choices = load_input(root, config["choices"])
    bootstrap = load_input(root, config["bootstrap"])

    def get_ranks(method, network):
        return ranks[method][network]

    result = panels.evaluate_cohort(
        config["cohort"],
        config["samples"],
        get_ranks,
        support,
        references,
        choices,
        bootstrap,
    )
    output = Path(config["output"])
    write_json(output if (output.is_absolute()) else root / output, result)

def select_portfolio(config_path):
    config, root = load_config(config_path)
    rows = load_input(root, config["development"])
    scores = {(row["method"], row["version"]): row["score"] for (row) in (rows)}
    check(len(scores) == len(rows), "Duplicate development candidate")
    result = policies.select_policies(
        scores, methods = tuple(config.get("methods", policies.METHODS))
    )
    output = Path(config["output"])
    write_json(output if (output.is_absolute()) else root / output, result)

def calculate_intervals(config_path):
    config, root = load_config(config_path)
    functions = {
        "selection": intervals.selection_intervals,
        "transfer": intervals.transfer_intervals,
        "capture": intervals.capture_interval,
    }
    arguments = load_input(root, config["arguments"])
    if (config["kind"] == "selection"):
        arguments["measurements"] = candidate_measurements(arguments["measurements"])
    if (config["kind"] == "transfer"):
        for (key) in (("string_measurements", "e_measurements")):
            arguments[key] = candidate_measurements(arguments[key])
    result = functions[config["kind"]](**arguments)
    output = Path(config["output"])
    write_json(output if (output.is_absolute()) else root / output, result)
