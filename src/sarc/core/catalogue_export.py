import math
from pathlib import Path

from .export import STATES, native_rank_rows
from .io import load_json, resolve_path, save_json

def nonempty_text(value, name):
    if (not isinstance(value, str) or not value):
        raise ValueError("A non-empty " + name + " is required")
    return value

def parameter_matches(method, parameter, specification):
    if (method == "PersonaDrive"):
        return parameter == "original"
    field, prefix, default = (
        ("mu", "mu_", 3) if (method == "DawnRank") else ("alpha", "alpha_", None)
    )
    value = specification.get(field, default)
    if (
        not parameter.startswith(prefix)
        or isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
    ):
        return False
    return value == float(parameter[len(prefix) :])

def check_job(job, cells):
    first = cells[0]
    factors = ("context", "method", "network_id", "parameter_id")
    if (any(any(cell[key] != first[key] for (key) in (factors)) for (cell) in (cells))):
        raise ValueError("One native job must map to one exact catalogue factor group")
    population = first["sample_order"]
    if (
        not population
        or len(population) != len(set(population))
        or job["models"] != population
        or any(cell["sample_order"] != population for (cell) in (cells))
        or job["expected_models"] != len(population)
        or job["context"] != first["context"]
        or job["method"] != first["method"]
        or job["mode"] != "native"
    ):
        raise ValueError("Native job and catalogue population or identity differ")
    if (any(not isinstance(sample, str) or not sample for (sample) in (population))):
        raise ValueError("The complete literal sample population is required")
    for (field) in (("parameter_id", "network_id")):
        if (field in job and job[field] != first[field]):
            raise ValueError(
                "Native job cannot be relabelled as another catalogue cell"
            )
    if (not parameter_matches(first["method"], first["parameter_id"], job)):
        raise ValueError("Native numerical parameter differs from the catalogue cell")
    seeds = [cell["master_seed"] for (cell) in (cells)]
    if (len(seeds) != len(set(seeds))):
        raise ValueError("Repeated catalogue seed mapping")
    if (first["method"] == "PRODIGY"):
        tasks = job["tasks"]
        task_keys = [(task["master_seed"], task["sample"]) for (task) in (tasks)]
        expected = {(seed, sample) for (seed) in (seeds) for (sample) in (population)}
        if (len(task_keys) != len(set(task_keys)) or set(task_keys) != expected):
            raise ValueError(
                "Every native task must map to exactly one catalogue state"
            )
        if (len(job["cells"]) != 1 or any(
            task["cell"] != job["cells"][0] for (task) in (tasks)
        )):
            raise ValueError("A native catalogue export requires one graph per job")
        for (task) in (tasks):
            effective = job | task
            if (
                effective["network_id"] != first["network_id"]
                or effective["context"] != first["context"]
                or not parameter_matches("PRODIGY", first["parameter_id"], effective)
            ):
                raise ValueError(
                    "Native task parameter or network differs from its cell"
                )
    elif (seeds != [None] or job["network_id"] != first["network_id"]):
        raise ValueError("Deterministic native jobs require one matching unseeded cell")

def check_prodigy_results(directory, job):
    tasks = {(task["master_seed"], task["sample"]): task for (task) in (job["tasks"])}
    seen = set()
    for (item) in (load_json(directory / "task_results.json")):
        value = load_json(directory / item["result"])
        key = (value["master_seed"], value["sample_id"])
        if (key not in tasks or key in seen):
            raise ValueError("Repeated or foreign native task output")
        seen.add(key)
        task = tasks[key]
        effective = job | task
        for (field) in (("cell", "sample_id", "master_seed", "status")):
            if (item[field] != value[field]):
                raise ValueError("Native result and task index identities differ")
        if (
            value["cell"] != task["cell"]
            or value["method"] != "PRODIGY"
            or value["network_id"] != effective["network_id"]
            or ("context" in value and value["context"] != effective["context"])
            or ("alpha" in value and value["alpha"] != effective["alpha"])
        ):
            raise ValueError("Native output provenance differs from the explicit job")
    if (seen != set(tasks)):
        raise ValueError("A declared native task output is missing")

def sample_record(row):
    state = row["status"]
    genes, scores = row["genes"], row["scores"]
    if (state not in STATES):
        raise ValueError("Unknown native catalogue status")
    available = state in ("SUCCESS", "SUCCESS_EMPTY")
    if (available):
        if (
            not isinstance(genes, list)
            or not isinstance(scores, list)
            or len(genes) != len(scores)
            or len(genes) != len(set(genes))
            or bool(genes) != (state == "SUCCESS")
            or any(not isinstance(gene, str) or not gene for (gene) in (genes))
            or any(
                isinstance(score, bool)
                or not isinstance(score, (int, float))
                or not math.isfinite(score)
                for (score) in (scores)
            )
        ):
            raise ValueError("Incomplete or invalid native catalogue ranking")
    elif (genes is not None or scores is not None):
        raise ValueError("Unavailable native ranks must remain null")
    return {
        "method_status": state,
        "genes": genes,
        "scores": scores,
        "reason": row.get("reason"),
        "source_id": row["source_id"],
    }

def catalogue_cells(config, base):
    from ..catalogue.matrix import validate_registry
    from ..catalogue.workflow import registry_index

    development = config["development"]
    if (type(development) is not bool):
        raise ValueError("The development/protected catalogue scope must be explicit")
    registry = registry_index(
        load_json(resolve_path(config["registry"], base)), development
    )
    if (not development):
        validate_registry(registry, next(iter(registry.values()))["sample_order"])
    if (any(
        len(cell["sample_order"]) != (36 if (development) else 396)
        for (cell) in (registry.values())
    )):
        raise ValueError(
            "Every catalogue cell must retain its complete sample population"
        )
    completed = {}
    directories = set()
    for (entry) in (config["jobs"]):
        directory = Path(resolve_path(entry["directory"], base))
        source_id = nonempty_text(entry["source_id"], "native source ID")
        identifiers = entry["cell_ids"]
        if (
            not isinstance(identifiers, list)
            or not identifiers
            or len(identifiers) != len(set(identifiers))
            or any(
                identifier not in registry or identifier in completed
                for (identifier) in (identifiers)
            )
            or directory in directories
        ):
            raise ValueError(
                "Native jobs and catalogue cell assignments must be unique"
            )
        directories.add(directory)
        cells = [registry[identifier] for (identifier) in (identifiers)]
        job = load_json(directory / "job.json")
        check_job(job, cells)
        if (job["method"] == "PRODIGY"):
            check_prodigy_results(directory, job)
        first = cells[0]
        records = native_rank_rows(
            directory, first["context"], first["network_id"], source_id
        )
        lookup = {}
        for (record) in (records):
            key = tuple(record["key"])
            if (
                len(key) != 5
                or key[:3] != (first["method"], first["context"], first["network_id"])
                or key[3:] in lookup
                or record["source_id"] != source_id
            ):
                raise ValueError("Native rank identity is duplicated or mismatched")
            lookup[key[3:]] = sample_record(record)
        expected = {
            (cell["master_seed"], sample)
            for (cell) in (cells)
            for (sample) in (cell["sample_order"])
        }
        if (set(lookup) != expected):
            raise ValueError("Native model-seed states cannot be omitted or mixed")
        for (cell) in (cells):
            completed[cell["cell_id"]] = {
                "cell_id": cell["cell_id"],
                "samples": {
                    sample: lookup[(cell["master_seed"], sample)]
                    for (sample) in (cell["sample_order"])
                },
            }
    if (set(completed) != set(registry)):
        raise ValueError(
            "Every registered catalogue cell must have a complete native output"
        )
    return [completed[identifier] for (identifier) in (registry)]

def export_catalogue(config_path):
    config_path = Path(config_path).resolve()
    config = load_json(config_path)
    if (config["operation"] != "catalogue_cells"):
        raise ValueError("The catalogue_cells export operation is required")
    result = catalogue_cells(config, config_path.parent)
    destination = Path(resolve_path(config["output"], config_path.parent))
    destination.parent.mkdir(parents = True, exist_ok = True)
    save_json(destination, result)
    return {
        "operation": "catalogue_cells",
        "output": str(destination),
        "cells": len(result),
    }
