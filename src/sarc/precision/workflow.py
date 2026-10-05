import os
import subprocess
from pathlib import Path

import numpy as np

from .io import (
    allocate,
    binary_matrix,
    exact_score,
    read_axes,
    read_csv,
    read_json,
    require,
    resolve,
    roster_rows,
    write_json,
)

def export_allocations(directory, roster, budgets):
    directory = Path(directory)
    manifest = read_json(directory / "scores.json")
    genes, samples = read_axes(directory)
    query = binary_matrix(directory / "mutation.bin", len(genes), len(samples))
    require(np.isin(query, [0, 1]).all(), "Mutation query is not binary")
    certificate_rows = read_csv(directory / "cutoff_certificates.csv")
    lookup = {
        (row["state"], row["sample_id"], int(row["budget"])): row for (row) in (certificate_rows)
    }
    require(len(lookup) == len(certificate_rows), "Duplicate numerical certificate")
    expected = {
        (state, sample, budget)
        for (state) in (manifest["states"])
        for (sample) in (samples)
        for (budget) in (budgets)
    }
    require(set(lookup) == expected, "Incomplete numerical certificates")
    require([row["sample_id"] for (row) in (roster)] == samples, "Roster order differs")
    roster = [
        {
            **row,
            "query": [genes[position] for (position) in (np.flatnonzero(query[:, index]))],
            "eligible_genes": [genes[position] for (position) in (np.flatnonzero(query[:, index]))],
        }
        for (index, row) in (enumerate(roster))
    ]
    allocations = []
    for (state) in (manifest["states"]):
        scores = binary_matrix(directory / (state + "_scores.bin"), len(genes), len(samples))
        require((scores >= 0).all(), "Negative native scores")
        for (index, model) in (enumerate(roster)):
            values = {
                genes[position]: exact_score(scores[position, index])
                for (position) in (np.flatnonzero(query[:, index]))
            }
            for (budget) in (budgets):
                weights, gap, threshold = allocate(values, budget)
                record = lookup[(state, model["sample_id"], budget)]
                allocations.append(
                    {
                        **model,
                        "state": state,
                        "cell": state,
                        "seed": "deterministic",
                        "budget": budget,
                        "weights": weights,
                        "numerically_certified": record["certified"] == "TRUE",
                        "stopped_by_threshold": record["stopped_by_threshold"] == "TRUE",
                        "cutoff_gap_exact": gap,
                        "cutoff_score_exact": threshold,
                    }
                )
    document = {
        "schema": "dawn_precision_allocations_v1",
        "method": "DawnRank",
        "models": roster,
        "genes": genes,
        "budgets": budgets,
        "master_seeds": [104729, 130363, 155921, 196613, 228017],
        "seed_semantics": "DETERMINISTIC_ALIASES",
        "allocations": allocations,
    }
    write_json(directory / "allocations.json", document)
    if (all((directory / (graph + "_degree.bin")).is_file() for (graph) in (("G0", "G1")))):
        degree_rows = []
        for (graph) in (("G0", "G1")):
            degrees = np.fromfile(directory / (graph + "_degree.bin"), dtype = "<f8")
            require(
                len(degrees) == len(genes)
                and np.isfinite(degrees).all()
                and (degrees >= 0).all(),
                "Invalid graph degrees",
            )
            for (index, model) in (enumerate(roster)):
                values = {
                    genes[position]: exact_score(degrees[position])
                    for (position) in (np.flatnonzero(query[:, index]))
                }
                for (budget) in (budgets):
                    weights = allocate(values, budget)[0]
                    degree_rows.append(
                        {**model, "state": graph, "budget": budget, "weights": weights}
                    )
        write_json(
            directory / "degree_allocations.json",
            {
                "models": roster,
                "allocations": degree_rows,
                "degree_definition": "Graph-wide degree restricted to the unchanged full molecular query, including zero-degree genes",
            },
        )
    return document

def run(config_path):
    config_path = Path(config_path).resolve()
    config = read_json(config_path)
    require(config["mode"] in ("factorial", "components"), "Unknown precision mode")
    require(config["epsilon"] in (1e-4, 1e-6, 1e-8, 1e-12), "Undeclared precision threshold")
    expected_limit = 100 if (config["epsilon"] == 1e-4) else 10000
    require(
        config["max_iterations"] == expected_limit,
        "Iteration limit differs from specified threshold study",
    )
    require(
        config["budgets"] and len(config["budgets"]) == len(set(config["budgets"])),
        "Repeated or missing budgets",
    )
    require(
        all(type(value) is int and value in (1, 5, 10, 20) for (value) in (config["budgets"])),
        "Unsupported allocation budget",
    )
    roster_path = resolve(config_path.parent, config["roster"])
    roster_document = read_json(roster_path)
    initial_roster = (
        roster_document["models"] if (isinstance(roster_document, dict)) else roster_document
    )
    planned_samples = [row["sample_id"] for (row) in (initial_roster)]
    roster_rows(
        roster_path, planned_samples, config["expected_models"], config["expected_donors"]
    )
    output = resolve(config_path.parent, config["output"])
    output.mkdir(parents = True, exist_ok = False)
    prepared = resolve(config_path.parent, config["prepared"])
    specification = {
        **config,
        "prepared": str(prepared),
        "output": str(output),
        "expected_samples": planned_samples,
    }
    if (config["mode"] == "components"):
        require(
            config["epsilon"] == 1e-12, "Component intervention requires the strict reference"
        )
        specification["reuse_directory"] = str(
            resolve(config_path.parent, config["reuse_directory"])
        )
    job = output / "job.json"
    write_json(job, specification)
    environment = dict(os.environ)
    environment.update(
        {
            key: "1"
            for (key) in ((
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS",
            ))
        }
    )
    temporary = output / "tmp"
    temporary.mkdir()
    environment.update({key: str(temporary) for (key) in (("TMPDIR", "TEMP", "TMP"))})
    if (config.get("library_paths")):
        environment["R_LIBS_USER"] = os.pathsep.join(
            str(resolve(config_path.parent, value)) for (value) in (config["library_paths"])
        )
    module = Path(__file__).resolve().parent
    executable = config["rscript"]
    if ("/" in executable or "\\" in executable):
        executable = str(resolve(config_path.parent, executable))
    with (output / "numerical.log").open("x", encoding = "utf-8") as log:
        subprocess.run(
            [executable, "--vanilla", str(module / "entry.R"), str(job), str(module)],
            stdout = log,
            stderr = subprocess.STDOUT,
            env = environment,
            check = True,
        )
    genes, samples = read_axes(output)
    roster = roster_rows(
        roster_path, samples, config["expected_models"], config["expected_donors"]
    )
    require(len(genes) == config["expected_genes"], "Output gene count differs")
    export_allocations(output, roster, config["budgets"])
    return {"output": str(output), "allocations": str(output / "allocations.json")}
