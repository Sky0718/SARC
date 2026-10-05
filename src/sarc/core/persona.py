import csv
import importlib.util
import json
import math
import os
import shutil
import sys
from pathlib import Path

def load_module(path, name):
    specification = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(value)
    return value

def write_json(path, value):
    with Path(path).open("x", encoding = "utf-8") as stream:
        json.dump(value, stream, indent = 2, allow_nan = False)
        stream.write("\n")

def write_csv(path, fields, rows):
    with Path(path).open("x", newline = "", encoding = "utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames = fields)
        writer.writeheader()
        writer.writerows(rows)

def run_persona(config):
    import networkx as nx
    import numpy as np
    import pandas as pd
    import tqdm

    versions = {
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "networkx": nx.__version__,
        "tqdm": tqdm.__version__,
    }
    expected = {
        "numpy": "1.26.4",
        "pandas": "2.2.3",
        "networkx": "2.8.8",
        "tqdm": "4.66.5",
    }
    if (
        versions != expected
        or os.environ.get("PYTHONHASHSEED") != "0"
        or pd.options.mode.copy_on_write
    ):
        raise ValueError(
            "PersonaDrive requires the declared package and iteration-order settings"
        )
    loader = load_module(config["loader_source"], "sarc_persona_loader")
    ranking = load_module(config["ranking_source"], "sarc_persona_rank")
    graph = load_module(
        Path(__file__).with_name("persona_graph.py"), "sarc_persona_graph"
    )
    edges = pd.read_csv(
        config["network"],
        sep = "\t",
        dtype = {"gene1": str, "gene2": str},
        keep_default_na = False,
    )
    if (list(edges.columns) != ["gene1", "gene2", "combined_score"]):
        raise ValueError("Expected the three-column retained network")
    ppi = {}
    for (first, second) in (edges[["gene1", "gene2"]].itertuples(index = False, name = None)):
        ppi.setdefault(first, set()).add(second)
        ppi.setdefault(second, set()).add(first)
    mutations = loader.load_mutations(config["mutations"])
    outliers = loader.load_DEGs(config["outliers"])
    mutations = mutations.drop(
        index = [gene for (gene) in (mutations.index) if (gene not in ppi)]
    )
    outliers = outliers.drop(
        index = [gene for (gene) in (outliers.index) if (gene not in ppi)]
    )
    samples = mutations.columns.tolist()
    if (samples != config["pipeline_sample_order"] or set(outliers.columns) != set(
        samples
    )):
        raise ValueError("The complete cohort or its source order differs")
    codec = config["pipeline_to_canonical"]
    if ([codec[sample] for (sample) in (samples)] != config["models"] or len(
        set(codec.values())
    ) != config["expected_models"]):
        raise ValueError("The sample codec must be the complete frozen bijection")
    if (any("_" in sample or "." in sample for (sample) in (samples))):
        raise ValueError(
            "Pipeline sample identifiers must preserve the original underscore parser"
        )
    output = Path(config["output"])
    os.chdir(output)
    (output / "data").mkdir()
    shutil.copyfile(config["pathways"], output / "data" / "kegg_pathways_v1.txt")
    cancer, dataset, network = config["cancer"], config["dataset"], config["network_id"]
    results = output / "results" / dataset / (cancer + "_" + network)
    results.mkdir(parents = True, exist_ok = False)
    graph.construct_pbns(mutations, outliers, ppi, cancer, dataset, network)
    ranking.calculate_pps("graphs/", dataset, cancer, network)
    similarity = pd.read_csv(
        results / "pps_matrix.csv", index_col = 0, float_precision = "legacy"
    )
    captured = {}
    original = ranking.calculate_infl_scores

    def observe(graph_value, genes, sample, pps, pathways):
        values = original(graph_value, genes, sample, pps, pathways)
        captured[sample] = values
        return values

    ranking.calculate_infl_scores = observe
    ranking.PersonaDrive(
        "graphs/", similarity, cancer, dataset, ranking.load_pathways(), network
    )
    rows = []
    statuses = []
    aliases = []
    for (sample) in (samples):
        ordered = sorted(
            captured.get(sample, {}).items(), key = lambda item: (-item[1], item[0])
        )
        if (any(not math.isfinite(float(score)) for ((gene, score)) in (ordered))):
            raise ValueError("PersonaDrive returned non-finite scores")
        for (index, (gene, score)) in (enumerate(ordered, 1)):
            rows.append(
                {
                    "context": config["context"],
                    "sample_id": codec[sample],
                    "method": "PersonaDrive",
                    "network_id": network,
                    "gene": gene,
                    "score": format(float(score), ".17g"),
                    "rank_position": index,
                }
            )
        state = "SUCCESS" if (ordered) else "SUCCESS_EMPTY"
        statuses.append(
            {
                "sample_id": codec[sample],
                "status": state,
                "candidate_count": len(ordered),
            }
        )
        aliases.extend(
            {
                "sample_id": codec[sample],
                "master_seed": seed,
                "status": state,
                "seed_interpretation": "DETERMINISTIC_REFERENCE",
            }
            for (seed) in (config["seeds"])
        )
    write_csv(
        output / "ranked_scores.csv",
        [
            "context",
            "sample_id",
            "method",
            "network_id",
            "gene",
            "score",
            "rank_position",
        ],
        rows,
    )
    write_csv(
        output / "sample_status.csv",
        ["sample_id", "status", "candidate_count"],
        statuses,
    )
    write_csv(
        output / "seed_aliases.csv",
        ["sample_id", "master_seed", "status", "seed_interpretation"],
        aliases,
    )
    write_json(
        output / "completion.json",
        {
            "method": "PersonaDrive",
            "models": len(samples),
            "network_id": network,
            "packages": versions,
        },
    )

def main():
    if (len(sys.argv) != 2):
        raise ValueError("One configuration path is required")
    config = json.loads(Path(sys.argv[1]).read_text(encoding = "utf-8"))
    try:
        run_persona(config)
    except Exception as error:
        write_json(
            Path(config["output"]) / "failure.json",
            {
                "status": "ERROR",
                "error_type": type(error).__name__,
                "reason": str(error),
            },
        )
        raise

if (__name__ == "__main__"):
    main()
