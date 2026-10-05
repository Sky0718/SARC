import os
import subprocess
import sys
from pathlib import Path

from .io import load_json, resolve_path, save_json, unique_axis, validate_seeds

PATH_FIELDS = {
    "counts",
    "mutations",
    "tumour",
    "normal",
    "mutation",
    "network",
    "normalisation_names",
    "normalise_source",
    "pathways",
    "explicit_evidence",
    "common_genes",
    "inherited_H0",
    "inherited_H1",
    "native_genes_H0",
    "native_genes_H1",
    "edges_G0",
    "edges_G1",
    "loader_source",
    "ranking_source",
    "outliers",
}

def resolved_inputs(value, base):
    if (isinstance(value, list)):
        return [resolved_inputs(item, base) for (item) in (value)]
    if (not isinstance(value, dict)):
        return value
    result = {}
    for (key, item) in (value.items()):
        if (key in PATH_FIELDS and isinstance(item, str)):
            result[key] = resolve_path(item, base)
        elif (key in ("original_sources", "library_paths")):
            result[key] = [resolve_path(path, base) for (path) in (item)]
        elif (
            key in ("rscript", "python_executable")
            and isinstance(item, str)
            and ("/" in item or "\\" in item)
        ):
            result[key] = resolve_path(item, base)
        else:
            result[key] = resolved_inputs(item, base)
    return result

def validate_prodigy_tasks(config):
    tasks = config["tasks"]
    cells = unique_axis(config["cells"], "cell")
    seeds = validate_seeds(config)
    secondary = config.get("parameter_scope", "primary") == "catalogue_secondary"
    if (secondary and len(cells) != 1):
        raise ValueError("A secondary job must represent one original native network")
    if (config["mode"] == "factorial" and cells != ["G0_H0", "G1_H0", "G0_H1", "G1_H1"]):
        raise ValueError("All four graph/evidence cells are required in their declared order")
    expected = {
        (cell, sample, seed)
        for (cell) in (cells)
        for (sample) in (config["models"])
        for (seed) in (seeds)
    }
    actual = [(task["cell"], task["sample"], task["master_seed"]) for (task) in (tasks)]
    if (len(actual) != len(set(actual)) or set(actual) != expected):
        raise ValueError("PRODIGY tasks must contain every cell, model and seed exactly once")
    for (task) in (tasks):
        if (any(
            task.get(field, config.get(field)) != config.get(field)
            for (field) in (("method", "mode", "parameter_scope", "context"))
        )):
            raise ValueError("Task overrides cannot change the scientific parameter scope")
        if (task.get("alpha", config["alpha"]) != config["alpha"]):
            raise ValueError("Task alpha differs from its declared parameter scope")
        if (secondary and task["network_id"] != config["network_id"]):
            raise ValueError("Secondary tasks must retain the declared native network")
        if (type(task["patient_seed"]) is not int or not 1 <= task["patient_seed"] <= 2147483646):
            raise ValueError("Supply the recorded integer patient_seed for every task")
        if (
            config["mode"] == "factorial"
            and "explicit_evidence" not in task
            and task.get("input_status") != "UNAVAILABLE"
        ):
            raise ValueError("Each factorial task requires explicit inherited evidence")
        if (config["mode"] == "native" and "explicit_evidence" in task):
            raise ValueError("Native evidence construction must be declared as native")

def run_job(config_path):
    config_path = Path(config_path).resolve()
    config = resolved_inputs(load_json(config_path), config_path.parent)
    validate_seeds(config)
    unique_axis(config["models"], "model")
    if (len(config["models"]) != config["expected_models"]):
        raise ValueError("The complete declared model population is required")
    if (config["mode"] not in ("native", "factorial")):
        raise ValueError("Unknown scientific input mode")
    if (config["method"] == "PRODIGY"):
        validate_prodigy_tasks(config)
    elif (config["method"] == "DawnRank"):
        pair = (config.get("epsilon"), config.get("max_iterations"))
        allowed = {(1e-4, 100), (1e-6, 10000), (1e-8, 10000), (1e-12, 10000)}
        if (pair not in allowed or type(config.get("max_iterations")) is not int):
            raise ValueError(
                "DawnRank requires explicit epsilon and max_iterations for its declared precision setting"
            )
    elif (config["method"] == "PersonaDrive" and config["mode"] != "native"):
        raise ValueError("PersonaDrive requires its full native cohort construction")
    elif (config["method"] not in ("DawnRank", "PersonaDrive")):
        raise ValueError("Unknown method")
    output = Path(resolve_path(config["output"], config_path.parent))
    output.mkdir(parents = True, exist_ok = False)
    config["output"] = str(output)
    serialised = output / "job.json"
    save_json(serialised, config)
    environment = dict(os.environ)
    environment.update(
        {
            "PYTHONDONTWRITEBYTECODE": "1",
            "PYTHONHASHSEED": "0",
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }
    )
    if (config.get("library_paths")):
        environment["R_LIBS_USER"] = os.pathsep.join(config["library_paths"])
    temporary = output / "tmp"
    temporary.mkdir()
    environment.update(
        {"TMPDIR": str(temporary), "TEMP": str(temporary), "TMP": str(temporary)}
    )
    if (config["method"] == "PersonaDrive"):
        command = [
            config.get("python_executable", sys.executable),
            str(Path(__file__).with_name("persona.py")),
            str(serialised),
        ]
    else:
        command = [
            config["rscript"],
            "--vanilla",
            str(Path(__file__).parent / "r" / "entry.R"),
            str(serialised),
        ]
    with (output / "native.log").open("x", encoding = "utf-8") as stream:
        process = subprocess.run(
            command,
            env = environment,
            stdout = stream,
            stderr = subprocess.STDOUT,
            check = False,
        )
    result = {
        "method": config["method"],
        "mode": config["mode"],
        "return_code": process.returncode,
        "output": str(output),
    }
    save_json(output / "execution.json", result)
    if (process.returncode != 0):
        raise RuntimeError("Native method returned an error; its output and log are preserved")
    return result
