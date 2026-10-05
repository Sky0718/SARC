import csv
from pathlib import Path

from .io import load_json, resolve_path, save_json

def write_matrix(path, matrix):
    with Path(path).open("x", newline = "", encoding = "utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["gene"] + list(matrix.samples))
        for (gene, values) in (zip(matrix.genes, matrix.values)):
            writer.writerow(
                [gene] + [format(float(value), ".17g") for (value) in (values)]
            )

def network_nodes(path):
    nodes = set()
    edges = set()
    with Path(path).open(newline = "", encoding = "utf-8-sig") as stream:
        reader = csv.DictReader(stream, delimiter = "\t")
        if (reader.fieldnames != ["gene1", "gene2", "combined_score"]):
            raise ValueError("Expected the original retained-network columns")
        for (row) in (reader):
            first, second = row["gene1"], row["gene2"]
            score = float(row["combined_score"])
            pair = tuple(sorted((first, second)))
            if (
                not first
                or not second
                or first == second
                or pair in edges
                or not 700 < score <= 1000
                or score != int(score)
            ):
                raise ValueError("Invalid retained network edge")
            edges.add(pair)
            nodes.update(pair)
    return nodes

def prepare_molecular(config_path):
    from . import molecular

    if (molecular.np.__version__ != "1.26.4" or molecular.pd.__version__ != "2.2.3"):
        raise ValueError("Molecular preparation requires NumPy 1.26.4 and pandas 2.2.3")
    config_path = Path(config_path).resolve()
    config = load_json(config_path)
    base = config_path.parent
    contract = load_json(resolve_path(config["contract"], base))
    observed = load_json(resolve_path(config["observed_genes"], base))
    if (
        not isinstance(observed, list)
        or len(observed) != len(set(observed))
        or any(not isinstance(gene, str) or not gene for (gene) in (observed))
    ):
        raise ValueError("An explicit unique observed-gene universe is required")
    quantitative = molecular.quantitative_inputs(
        resolve_path(config["source_root"], base), contract, set(observed)
    )
    mutation = molecular.load_matrix(
        resolve_path(config["mutation"], base), config["mutation_sample_map"]
    )
    mutation, excluded = molecular.remove_unresolved(mutation, binary = True)
    models = tuple(contract["original_tumour_order"])
    if (len(models) != config["expected_models"] or set(mutation.samples) != set(models)):
        raise ValueError("The complete mutation and tumour populations must agree")
    mutation = mutation.select(mutation.genes, models)
    nodes = network_nodes(resolve_path(config["network"], base))
    dawn = molecular.dawn_projection(quantitative, mutation, nodes)
    prodigy = molecular.prodigy_projection(quantitative, mutation, nodes)
    output = Path(resolve_path(config["output"], base))
    output.mkdir(parents = True, exist_ok = False)
    for (name, values) in ((("dawn", dawn), ("prodigy", prodigy))):
        directory = output / name
        directory.mkdir()
        for (field) in (
            ("tumour", "normal", "mutation")
            if (name == "dawn")
            else ("counts", "mutation")
        ):
            write_matrix(directory / (field + ".csv"), values[field])
        save_json(directory / "genes.json", values["genes"])
    with (output / "dawn" / "normalisation_names.csv").open(
        "x", newline = "", encoding = "utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(["role", "canonical", "original"])
        for (role) in (("tumour", "normal")):
            writer.writerows(
                (role, canonical, original)
                for ((canonical, original)) in (
                    zip(dawn[role].samples, dawn["normalisation_names"][role])
                )
            )
    save_json(
        output / "prodigy" / "sample_origins.json", list(prodigy["sample_origins"])
    )
    record = {
        "models": list(models),
        "dawn_genes": len(dawn["genes"]),
        "prodigy_genes": len(prodigy["genes"]),
        "normal_samples": list(dawn["normal"].samples),
        "count_samples": list(prodigy["counts"].samples),
        "excluded": {**quantitative["excluded"], "mutation": excluded},
    }
    save_json(output / "prepared.json", record)
    return record
