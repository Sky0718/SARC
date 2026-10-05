import csv
from collections import Counter
from pathlib import Path
from ..core.io import load_json, resolve_path
from ..evaluation import localisation
from ..evaluation.prediction import require, validate_rows
from . import features_r2, features_r7, mechanism, specs_r1

def edges(path):
    with Path(path).open(encoding = "utf-8-sig", newline = "") as stream:
        records = list(csv.DictReader(stream, delimiter = "\t"))
    result = {
        (row["gene1"], row["gene2"]): int(row["combined_score"]) for (row) in (records)
    }
    require(len(result) == len(records), "Duplicate canonical graph edge")
    mechanism.canonical_edges(result)
    return result

def extract(config, base):
    round_id = config["round"]
    require(
        round_id in ("R2", "R5", "R7"),
        "Only the recorded localisation, mechanism and invariant feature rounds are supported",
    )
    slots = load_json(resolve_path(config["slots"], base))
    molecular = load_json(resolve_path(config["molecular"], base))
    require(
        set(molecular) == set(specs_r1.COHORTS),
        "Both full development cohorts required",
    )
    require(
        len(slots)
        == len({(slot["cohort"], slot["network_id"]) for (slot) in (slots)})
        == 128,
        "Complete unique 128 graph slots required",
    )
    require(
        Counter(slot["cohort"] for (slot) in (slots))
        == dict.fromkeys(specs_r1.COHORTS, 64),
        "Full graph population differs",
    )
    for (cohort, record) in (molecular.items()):
        population = record["population"]
        require(
            len(population) == len(set(population)) == 36
            and set(record["native_mutations"])
            == set(record["degs"])
            == set(population),
            "Complete aligned 36-sample molecular inputs required",
        )
        require(
            record["literal_normalisation"]
            == "ORIGINAL_NATIVE_INVALID_LABELS_EXCLUDED",
            "Original native molecular-label semantics required",
        )
    graphs, baselines, cache, rows = {}, {}, {}, []
    pathways = (
        mechanism.parse_pathways(
            Path(resolve_path(config["pathways"], base)).read_text(encoding = "utf-8")
        )
        if (round_id in ("R5", "R7"))
        else None
    )
    for (slot) in (slots):
        cohort = slot["cohort"]
        require(
            slot["transition"] in specs_r1.TRANSITIONS
            and slot["control_seed"] in (None, 271828, 314159),
            "Protected or undeclared graph slot",
        )
        require(
            slot["baseline_network_id"]
            in {slot["transition"] + "__persistent_" + end for (end) in (("old", "new"))},
            "Original persistent baseline identity required",
        )
        record = molecular[cohort]
        population = record["population"]
        mutation = {
            sample: set(genes) for (sample, genes) in (record["native_mutations"].items())
        }
        degs = {sample: set(genes) for (sample, genes) in (record["degs"].items())}
        graph_paths = [
            resolve_path(slot[name], base)
            for (name) in (("baseline_edges", "changed_edges"))
        ]
        for (path) in (graph_paths):
            if (path not in graphs):
                graphs[path] = edges(path)
        old, new = [graphs[path] for (path) in (graph_paths)]
        if (round_id == "R2"):
            prepared = features_r2.prepare_graph(old, new)
            blocks = {
                sample: features_r2.strict_block(
                    prepared, mutation[sample], degs[sample]
                )
                for (sample) in (population)
            }
            ranking_path = resolve_path(slot["old_rankings"], base)
            if (ranking_path not in baselines):
                baselines[ranking_path] = load_json(ranking_path)
            ranking = baselines[ranking_path]
            require(
                ranking["cohort"] == cohort
                and ranking["network_id"] == slot["baseline_network_id"]
                and ranking["population"] == population,
                "Old native ranking identity differs",
            )
        else:
            key = (cohort, graph_paths[0])
            if (key not in cache):
                cache[key] = mechanism.prepare_baseline(
                    old, population, mutation, degs, pathways
                )
            blocks = (
                mechanism.mechanism_features(cache[key], new)
                if (round_id == "R5")
                else features_r7.extract_features(cache[key], new)
            )
        for (method) in (specs_r1.METHODS):
            masters = list(specs_r1.MASTERS) if (method == "PRODIGY") else [None]
            for (sample) in (population):
                row = {
                    "row_id": "::".join((cohort, slot["network_id"], method, sample)),
                    "cohort": cohort,
                    "network_id": slot["network_id"],
                    "method": method,
                    "sample_id": sample,
                    "transition": slot["transition"],
                    "control_seed": slot["control_seed"],
                    "master_seeds": masters,
                }
                if (round_id == "R2"):
                    seeded = ranking["methods"][method]
                    require(
                        [entry["master_seed"] for (entry) in (seeded)] == masters
                        and all(
                            set(entry["samples"]) == set(population) for (entry) in (seeded)
                        ),
                        "Complete old baseline algorithm seeds and samples required",
                    )
                    row.update(
                        strict_input = blocks[sample],
                        old_output_assisted = features_r2.assisted_block(
                            prepared,
                            mutation[sample],
                            {
                                entry["master_seed"]: entry["samples"][sample]
                                for (entry) in (seeded)
                            },
                            method,
                        ),
                    )
                else:
                    row[
                        "mechanism_input" if (round_id == "R5") else "invariant_input"
                    ] = blocks[sample]
                rows.append(row)
    require(
        len(rows) == len({row["row_id"] for (row) in (rows)}) == 13824,
        "Complete unique feature row population required",
    )
    return features_r2.plain(rows)

def join(original, features, round_id):
    if (round_id == "R2"):
        return localisation.join_rows(original, features)[0]
    require(round_id in ("R5", "R7"), "Undeclared feature join")
    validate_rows(original)
    indexed = {row["row_id"]: row for (row) in (features)}
    require(
        len(indexed) == len(features) == len(original)
        and set(indexed) == {row["row_id"] for (row) in (original)},
        "Exact complete feature population required",
    )
    names = mechanism.FEATURE_NAMES if (round_id == "R5") else features_r7.FEATURE_NAMES
    block_name = "mechanism_input" if (round_id == "R5") else "invariant_input"
    result, shared = [], {}
    for (row) in (original):
        extra = indexed[row["row_id"]]
        require(
            all(row[name] == extra[name] for (name) in (localisation.IDENTITIES)),
            "Immutable row metadata differs",
        )
        block = extra[block_name]
        require(
            block["status"] == "AVAILABLE"
            and block["reason"] is None
            and set(block["features"]) == set(names),
            "Complete finite input-only feature block required",
        )
        values = {name: specs_r1.number(block["features"][name]) for (name) in (names)}
        require(
            round_id != "R5" or all(value >= 0 for (value) in (values.values())),
            "Negative absolute mechanism feature",
        )
        key = (row["cohort"], row["network_id"], row["sample_id"])
        require(
            key not in shared or shared[key] == block,
            "Input-only features differ between methods",
        )
        shared[key] = block
        require(
            not set(names).intersection(row["features"]),
            "Original feature overwrite forbidden",
        )
        result.append({**row, "features": {**row["features"], **values}})
    validate_rows(result)
    return result
