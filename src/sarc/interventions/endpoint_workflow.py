from itertools import product

from ..analysis.run import run_operation
from . import endpoints
from . import features
from .original_prediction import MODELS, fit_family

def require(value, message):
    if (not value):
        raise ValueError(message)

def indexed_ranks(rows):
    result = {}
    for (row) in (rows):
        key = (row["network_id"], row["method"], row["sample_id"], row["master_seed"])
        require(key not in result, "Duplicate native rank identity")
        result[key] = {name: row[name] for (name) in (("status", "genes", "scores"))}
    return result

def evaluate_slice(document):
    cohort, population = document["cohort"], document["population"]
    require(cohort in endpoints.COHORTS, "Unknown cohort")
    require(
        document["transition"] in ("11_0_to_11_5", "11_5_to_12_0")
        and document["direction"] in ("insert", "retract"),
        "Unknown transition or direction",
    )
    require(
        len(population) == len(set(population)) == endpoints.COHORTS[cohort],
        "Complete ordered cohort is required",
    )
    graphs = document["graphs"]
    expected = {("real", label, None) for (label) in (endpoints.CLASSES)}
    expected.update(
        ("matched_control", label, seed)
        for ((label, seed)) in (product(endpoints.CLASSES[:-1], endpoints.CONTROLS))
    )
    actual = [(row["kind"], row["class_id"], row["control_seed"]) for (row) in (graphs)]
    require(
        len(actual) == len(set(actual)) and set(actual) == expected,
        "All 21 real and compact control slots are required",
    )
    require(
        len({row["network_id"] for (row) in (graphs)}) == len(graphs),
        "Duplicate graph identity",
    )
    reference = document["references"]
    require(
        set(reference) == set(endpoints.METHODS),
        "Complete method reference mapping is required",
    )
    require(
        all(
            set(reference[method]) == set(population)
            for (method) in (endpoints.METHODS)
        ),
        "Complete biological reference roster is required",
    )
    ranks = indexed_ranks(document["ranks"])
    baseline_id = document["baseline_network_id"]
    require(
        baseline_id not in {row["network_id"] for (row) in (graphs)},
        "Baseline and edited graph identities must differ",
    )
    used, records = set(), {}
    for (graph) in (graphs):
        require(
            graph["input_status"] in ("READY", "FAILED_MATCHING"),
            "Explicit graph availability is required",
        )
        require(
            graph["kind"] != "real" or graph["input_status"] == "READY",
            "A real intervention cannot be a failed matching control",
        )
        require(
            graph["kind"] != "matched_control" or bool(graph.get("matching_evidence")),
            "Control matching diagnostics are required",
        )
        for (method, sample) in (product(endpoints.METHODS, population)):
            baseline, intervention = {}, {}
            for (seed) in (endpoints.seeds(method)):
                before = (baseline_id, method, sample, seed)
                require(before in ranks, "Missing baseline rank state")
                baseline[seed] = ranks[before]
                used.add(before)
                after = (graph["network_id"], method, sample, seed)
                if (graph["input_status"] == "READY"):
                    require(after in ranks, "Missing intervention rank state")
                    intervention[seed] = ranks[after]
                    used.add(after)
                else:
                    require(
                        after not in ranks,
                        "Failed matching cannot contain an inferred method output",
                    )
                    intervention[seed] = {
                        "status": "UNAVAILABLE",
                        "genes": None,
                        "scores": None,
                    }
            records[(graph["network_id"], method, sample)] = endpoints.compare_rankings(
                method, baseline, intervention, reference[method][sample]
            )
    require(used == set(ranks), "Unmatched native rank records remain")
    counts = endpoints.validate_draws(document["bootstrap_indices"], len(population))
    panels = endpoints.panel_summaries(population, graphs, records, counts)
    return {
        "cohort": cohort,
        "transition": document["transition"],
        "direction": document["direction"],
        "population": population,
        "graphs": graphs,
        "panels": panels,
        "individual_records": [
            {
                "network_id": key[0],
                "method": key[1],
                "sample_id": key[2],
                "result": value,
            }
            for ((key, value)) in (records.items())
        ],
    }

def fit_original_predictors(rows):
    return {name: fit_family(rows, name) for (name) in (MODELS)}

def input_features(baseline, changed, modules, mutations, degs, method):
    graphs = []
    for (rows) in ((baseline, changed)):
        graph = {}
        for (row) in (rows):
            edge = (row["gene1"], row["gene2"])
            require(edge not in graph, "Repeated canonical input edge")
            graph[edge] = row["combined_score"]
        graphs.append(graph)
    return features.graph_features(*graphs, modules, set(mutations), set(degs), method)

def run(config_path):
    return run_operation(
        config_path,
        {
            "edit_slice": evaluate_slice,
            "original_predictors": fit_original_predictors,
            "input_features": input_features,
        },
    )
