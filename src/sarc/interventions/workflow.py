import math
import random
from pathlib import Path

import networkx as nx

from ..resources.io import check, read_json, read_rows, write_json
from . import matching as controls
from . import balance
from . import graph_edits as real

POPULATIONS = {"COAD_CCLE": 36, "LUAD_CCLE": 36, "COAD_TCGA": 396}
TRANSITIONS = ("11_0_to_11_5", "11_5_to_12_0")

def path_from(config_path, value):
    path = Path(value)
    return path if (path.is_absolute()) else Path(config_path).resolve().parent / path

def save_graph(output, name, rows):
    destination = output / "graphs" / (name + ".tsv")
    destination.parent.mkdir(parents = True, exist_ok = True)
    destination.write_bytes(real.graph_bytes(rows))
    return {"network_id": name, "edges": str(destination), "edge_count": len(rows)}

def prepare_graphs(config_path):
    config = read_json(config_path)
    check(
        nx.__version__ == "2.8.8", "NetworkX 2.8.8 is required for the fixed partition"
    )
    check(
        set(config["cohorts"]) == set(POPULATIONS),
        "All three complete cohorts are required",
    )
    check(
        tuple(row["id"] for (row) in (config["transitions"])) == TRANSITIONS,
        "Both consecutive transitions are required",
    )
    output = path_from(config_path, config["output"])
    output.mkdir(parents = True, exist_ok = False)
    mutations = {}
    for (cohort, record) in (config["cohorts"].items()):
        samples = record["samples"]
        check(
            len(samples) == len(set(samples)) == POPULATIONS[cohort],
            "The sample roster is incomplete",
        )
        mutations[cohort] = real.mutation_summary(
            path_from(config_path, record["mutation"]), samples
        )
    records = []
    for (transition) in (config["transitions"]):
        old = real.edges(path_from(config_path, transition["baseline"]))
        new = real.edges(path_from(config_path, transition["followup"]))
        nodes = {
            row["preferred_name"]
            for (row) in (
                read_rows(path_from(config_path, transition["persistent_nodes"]))
            )
        }
        check(
            {gene for (row) in (old + new) for (gene) in (row[:2])} <= nodes,
            "Graph endpoints exceed persistent identities",
        )
        graph = nx.Graph()
        graph.add_nodes_from(sorted(nodes))
        graph.add_edges_from(sorted(row[:2] for (row) in (old)))
        groups = nx.community.louvain_communities(
            graph, weight = None, resolution = 1, threshold = 1e-7, seed = 104729
        )
        groups = sorted(tuple(sorted(group)) for (group) in (groups))
        modules = {
            gene: index for ((index, group)) in (enumerate(groups)) for (gene) in (group)
        }
        prefix = transition["id"]
        write_json(output / (prefix + "__modules.json"), {"modules": groups})
        changed = {row[:2] for (row) in (old)} ^ {row[:2] for (row) in (new)}
        for (cohort, (counts, _)) in (mutations.items()):
            threshold = math.ceil(POPULATIONS[cohort] * 0.05)
            anchors = {gene for (gene) in (nodes) if (counts.get(gene, 0) >= threshold)}
            classes = {
                edge: real.classify(edge, anchors, modules) for (edge) in (changed)
            }
            write_json(
                output / (prefix + "__" + cohort + "__classes.json"),
                {
                    "anchors": sorted(anchors),
                    "anchor_threshold": threshold,
                    "edits": [
                        {"gene1": a, "gene2": b, "class": classes[(a, b)]}
                        for ((a, b)) in (sorted(classes))
                    ],
                },
            )
            for (direction, baseline, other) in ((
                ("insert", old, new),
                ("retract", new, old),
            )):
                for (class_id, selected) in (real.CLASS_SETS.items()):
                    rows, edits = real.intervene(baseline, other, classes, selected)
                    name = "XB__" + "__".join((prefix, cohort, direction, class_id))
                    records.append(save_graph(output, name, rows) | {"edits": edits})
    write_json(output / "graphs.json", records)

def prepare_matched_controls(config_path):
    config = read_json(config_path)
    check(config["transition"] in TRANSITIONS, "Unknown transition")
    check(config["cohort"] in POPULATIONS, "Unknown cohort")
    check(config["direction"] in ("insert", "retract"), "Unknown direction")
    seed_ids = [row["seed"] for (row) in (config["seeds"])]
    check(
        seed_ids == [271828, 314159, 173205],
        "The compact control population contains three fixed seeds",
    )
    output = path_from(config_path, config["output"])
    output.mkdir(parents = True, exist_ok = False)
    baseline_rows = real.edges(path_from(config_path, config["baseline"]))
    baseline = {row[:2]: row[2] for (row) in (baseline_rows)}
    observed = {
        row[:2]: row[2]
        for (row) in (real.edges(path_from(config_path, config["followup"])))
    }
    degrees = {
        row["preferred_name"]: 0
        for (row) in (read_rows(path_from(config_path, config["persistent_nodes"])))
    }
    for (edge) in (baseline):
        for (node) in (edge):
            degrees[node] += 1
    assignments = read_json(path_from(config_path, config["classes"]))
    classes = {
        (row["gene1"], row["gene2"]): row["class"] for (row) in (assignments["edits"])
    }
    removed = {
        label: {
            edge
            for (edge) in (baseline.keys() - observed.keys())
            if (classes[edge] == label)
        }
        for (label) in (controls.CLASSES)
    }
    added = {
        label: {
            edge
            for (edge) in (observed.keys() - baseline.keys())
            if (classes[edge] == label)
        }
        for (label) in (controls.CLASSES)
    }
    records = []
    for (seed) in (config["seeds"]):
        check(
            isinstance(seed["random_seed"], int),
            "An exact recorded generator seed is required",
        )
        rng = random.Random(seed["random_seed"])
        removal, removal_info = controls.jointly_match_removals(
            baseline, removed, degrees, rng
        )
        addition, addition_info = controls.jointly_match_additions(
            baseline, observed, added, degrees, rng
        )
        edits = {
            "removed": {key: sorted(value) for ((key, value)) in (removal.items())},
            "added": {
                key: [(*edge, score) for ((edge, score)) in (sorted(value.items()))]
                for ((key, value)) in (addition.items())
            },
        }
        write_json(output / (str(seed["seed"]) + "__edits.json"), edits)
        try:
            result = balance.verify_matches(
                baseline, observed, removed, added, removal, addition, degrees
            )
        except ValueError as error:
            write_json(
                output / (str(seed["seed"]) + "__failure.json"),
                {
                    "error": str(error),
                    "seed": seed["seed"],
                    "removal": removal_info,
                    "addition": addition_info,
                },
            )
            raise
        for (label) in ((*controls.CLASSES, "PB_PW")):
            selected = ("PB", "PW") if (label == "PB_PW") else (label,)
            remove = set().union(*(removal[key] for (key) in (selected)))
            add = {
                edge: score
                for (key) in (selected)
                for ((edge, score)) in (addition[key].items())
            }
            rows = [row for (row) in (baseline_rows) if (row[:2] not in remove)] + [
                (*edge, add[edge]) for (edge) in (sorted(add))
            ]
            check(
                len(rows) == len(baseline) - len(remove) + len(add),
                "Control graph multiplicity differs",
            )
            name = "XB__" + "__".join(
                (
                    config["transition"],
                    config["cohort"],
                    config["direction"],
                    "control_" + str(seed["seed"]),
                    label,
                )
            )
            records.append(save_graph(output, name, rows) | {"balance": result[label]})
        write_json(output / (str(seed["seed"]) + "__balance.json"), result)
    write_json(output / "graphs.json", records)
