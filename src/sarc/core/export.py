import csv
import math
from pathlib import Path

from .io import SEEDS, load_json, resolve_path, save_json, validate_seeds

CELLS = ("G0_H0", "G1_H0", "G0_H1", "G1_H1")
STATES = ("SUCCESS", "SUCCESS_EMPTY", "FAILED", "UNAVAILABLE")

def read_csv(path):
    with Path(path).open(newline = "", encoding = "utf-8-sig") as stream:
        return list(csv.DictReader(stream))

def number(token, integer = False):
    if (token in ("", "NA", "NaN")):
        return None
    value = int(token) if (integer) else float(token)
    if (not math.isfinite(value)):
        raise ValueError("Non-finite native diagnostic")
    return value

def native_rows(directory):
    directory = Path(directory)
    scores = {}
    for (row) in (read_csv(directory / "ranked_scores.csv")):
        key = (row.get("cell", "G0_H0"), row["sample_id"])
        scores.setdefault(key, []).append(row)
    output = {}
    for (state) in (read_csv(directory / "sample_status.csv")):
        key = (state.get("cell", "G0_H0"), state["sample_id"])
        if (key in output or state["status"] not in STATES):
            raise ValueError("Repeated native state or unsupported status")
        rows = scores.pop(key, [])
        if ([int(row["rank_position"]) for (row) in (rows)] != list(
            range(1, len(rows) + 1)
        )):
            raise ValueError("Native rank order is incomplete")
        if (state["status"] in ("FAILED", "UNAVAILABLE") and rows):
            raise ValueError("Unavailable state contains native ranks")
        if (state["status"] in ("SUCCESS", "SUCCESS_EMPTY") and (
            state["status"] == "SUCCESS"
        ) != bool(rows)):
            raise ValueError("Native empty/nonempty status differs from rank rows")
        if (int(state["candidate_count"]) != len(rows)):
            raise ValueError("Native candidate count differs from complete rank rows")
        output[key] = (state, rows)
    if (scores):
        raise ValueError("Native ranks without a corresponding state")
    return output

def dawn_cells(directory):
    result = {cell: {} for (cell) in (CELLS)}
    job = load_json(Path(directory) / "job.json")
    if (job["method"] != "DawnRank" or job["mode"] != "factorial"):
        raise ValueError("DawnRank factorial output is required")
    for ((cell, sample), (state, rows)) in (native_rows(directory).items()):
        if (cell not in result):
            raise ValueError("Native cell lies outside the four-cell design")
        result[cell][sample] = {
            "status": state["status"],
            "reason": state["reason"],
            "scores": [
                [row["gene"], float(row["score"]), float(row["native_percent_rank"])]
                for (row) in (rows)
            ],
            "iterations": number(state["iterations"], True),
            "terminal_update_norm": number(state["terminal_update_norm"]),
            "termination_reason": state["termination_reason"],
            "fixed_point_error_certified": state["fixed_point_error_certified"].upper()
            == "TRUE",
        }
    if (any(set(rows) != set(job["models"]) for (rows) in (result.values()))):
        raise ValueError("Full model-cell population is required")
    return result

def prodigy_carriers(directory, scope, cell_ids, source_id):
    directory = Path(directory)
    job = load_json(directory / "job.json")
    if (job["method"] != "PRODIGY" or job["mode"] != "factorial" or not source_id):
        raise ValueError(
            "PRODIGY factorial output and an explicit source ID are required"
        )
    models = {model["sample_id"]: model for (model) in (scope["models"])}
    if (job["models"] != list(models)):
        raise ValueError("Scientific scope and native ordered model population differ")
    identifiers = {
        (row["model_id"], row["graph"], row["evidence"]): row["signature_id"]
        for (row) in (cell_ids)
    }
    if (len(identifiers) != len(cell_ids) or len(set(identifiers.values())) != len(
        cell_ids
    )):
        raise ValueError("Cell IDs must be explicit and unique")
    groups = {}
    for (item) in (load_json(directory / "task_results.json")):
        value = load_json(directory / item["result"])
        graph_key, evidence = item["cell"].split("_")
        graph = scope["graph_states"][("G0", "G1").index(graph_key)]
        model = models[item["sample_id"]]
        key = (model["model_id"], graph, evidence)
        for (field) in (("cell", "sample_id", "master_seed", "status")):
            if (value[field] != item[field]):
                raise ValueError("Task index and native result identities differ")
        if (value["status"] not in STATES):
            raise ValueError("Unknown native PRODIGY status")
        unavailable = value["status"] in ("FAILED", "UNAVAILABLE")
        genes = None if (unavailable) else value["genes"]
        bits = None if (unavailable) else value["scores_little_endian_float64_hex"]
        if (not unavailable and len(bits) != 16 * len(genes)):
            raise ValueError("Lossless score encoding does not match gene count")
        state = {
            "master_seed": item["master_seed"],
            "status": value["status"],
            "reason": value.get("reason"),
            "genes": genes,
            "scores_binary64_le": None
            if (unavailable)
            else [bits[index : index + 16] for (index) in (range(0, len(bits), 16))],
            "evidence_id": model["evidence_id"][evidence],
            "native_binding": {"source_id": source_id, "task": item["task"]},
        }
        if (not state["reason"] and value["status"] != "SUCCESS"):
            raise ValueError("Empty or unavailable output requires its native reason")
        groups.setdefault(key, []).append(state)
    expected = {
        (model["model_id"], graph, evidence)
        for (model) in (models.values())
        for (graph) in (scope["graph_states"])
        for (evidence) in (scope["evidence_states"])
    }
    if (set(groups) != expected or set(identifiers) != expected):
        raise ValueError("Complete model-cell identities are required")
    cells = []
    for (model) in (models.values()):
        for (graph) in (scope["graph_states"]):
            for (evidence) in (scope["evidence_states"]):
                key = (model["model_id"], graph, evidence)
                states = sorted(groups[key], key = lambda state: state["master_seed"])
                if ([state["master_seed"] for (state) in (states)] != list(SEEDS)):
                    raise ValueError("Exactly five complete seed states are required")
                cells.append(
                    {
                        "model_id": model["model_id"],
                        "sample_id": model["sample_id"],
                        "donor_id": model["donor_id"],
                        "root_axis_id": model["root_axis_id"],
                        "graph": graph,
                        "evidence": evidence,
                        "signature_id": identifiers[key],
                        "common_domain_id": scope["common_domain_id"],
                        "pathway_dictionary_id": scope["pathway_dictionary_id"],
                        "graph_id": scope["graph_id"][graph],
                        "evidence_id": model["evidence_id"][evidence],
                        "seed_states": states,
                    }
                )
    return {
        "schema": "common_domain_factorial_native_carriers_v1",
        "data_class": scope["data_class"],
        "declared_complete": True,
        "cells": cells,
    }

def native_rank_rows(directory, pool, release, source_id):
    directory = Path(directory)
    job = load_json(directory / "job.json")
    if (job["mode"] != "native" or not source_id):
        raise ValueError("Native results and an explicit source ID are required")
    seeds = validate_seeds(job)
    records = []
    if (job["method"] == "PRODIGY"):
        if (len(job["cells"]) != 1):
            raise ValueError("Export one native release job at a time")
        for (item) in (load_json(directory / "task_results.json")):
            value = load_json(directory / item["result"])
            unavailable = value["status"] in ("FAILED", "UNAVAILABLE")
            records.append(
                {
                    "key": [
                        "PRODIGY",
                        pool,
                        release,
                        value["master_seed"],
                        value["sample_id"],
                    ],
                    "status": value["status"],
                    "reason": value.get("reason"),
                    "genes": None if (unavailable) else value["genes"],
                    "scores": None if (unavailable) else value["scores"],
                    "source_id": source_id,
                }
            )
    else:
        for ((cell, sample), (state, rows)) in (native_rows(directory).items()):
            unavailable = state["status"] in ("FAILED", "UNAVAILABLE")
            records.append(
                {
                    "key": [job["method"], pool, release, None, sample],
                    "status": state["status"],
                    "reason": state.get("reason"),
                    "genes": None
                    if (unavailable)
                    else [row["gene"] for (row) in (rows)],
                    "scores": None
                    if (unavailable)
                    else [float(row["score"]) for (row) in (rows)],
                    "source_id": source_id,
                }
            )
    expected = {
        (seed, sample)
        for (seed) in (seeds if (job["method"] == "PRODIGY") else (None,))
        for (sample) in (job["models"])
    }
    keys = [(row["key"][3], row["key"][4]) for (row) in (records)]
    if (len(keys) != len(set(keys)) or set(keys) != expected):
        raise ValueError("Full native model-seed population is required")
    return records

def export_results(config_path):
    config_path = Path(config_path).resolve()
    config = load_json(config_path)
    base = config_path.parent
    operation = config["operation"]
    if (operation == "dawn_cells"):
        result = dawn_cells(resolve_path(config["directory"], base))
    elif (operation == "prodigy_carriers"):
        result = prodigy_carriers(
            resolve_path(config["directory"], base),
            load_json(resolve_path(config["scope"], base)),
            load_json(resolve_path(config["cell_ids"], base)),
            config["source_id"],
        )
    elif (operation == "native_ranks"):
        result = []
        for (job) in (config["jobs"]):
            result.extend(
                native_rank_rows(
                    resolve_path(job["directory"], base),
                    job["pool"],
                    job["release"],
                    job["source_id"],
                )
            )
        keys = [tuple(row["key"]) for (row) in (result)]
        if (len(keys) != len(set(keys))):
            raise ValueError("Duplicate native rank group across jobs")
    elif (operation in ("continuous_selections", "eight_model_eligibility")):
        from .selections import continuous_selections, eight_model_eligibility

        result = (
            continuous_selections
            if (operation == "continuous_selections")
            else eight_model_eligibility
        )(config, base)
    else:
        raise ValueError("Unsupported native result export")
    destination = Path(resolve_path(config["output"], base))
    destination.parent.mkdir(parents = True, exist_ok = True)
    save_json(destination, result)
    return {"operation": operation, "output": str(destination)}
