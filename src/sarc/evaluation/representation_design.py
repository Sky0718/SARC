import itertools
import math
from fractions import Fraction

COHORTS = ("COAD_CCLE", "LUAD_CCLE")
TRANSITIONS = ("11_0_to_11_5", "11_5_to_12_0")
DIRECTIONS = ("insert", "retract")
SEEDS = (104729, 130363, 155921, 196613, 228017)
REPRESENTATIONS = ("ORIGINAL", "PERMUTED", "CANONICAL")
REPEATS = (1, 2)

def binding(value):
    if (
        not isinstance(value, dict)
        or not isinstance(value.get("source_id"), str)
        or (not value["source_id"])
    ):
        raise ValueError("Explicit scientific source identity required")
    return {"source_id": value["source_id"]}

def choose_samples(records, selected_samples):
    if (set(selected_samples) != set(COHORTS)):
        raise ValueError("Both archived cohort selections are required")
    selected = {}
    for (cohort) in (COHORTS):
        matches = [record for (record) in (records) if (record["context"] == cohort)]
        if (not matches):
            raise ValueError("Missing complete cohort metadata")
        roster = matches[0]["tumours"]
        if (
            len(roster) != 36
            or len(set(roster)) != 36
            or any((record["tumours"] != roster for (record) in (matches)))
        ):
            raise ValueError("Complete unchanged 36-sample roster required")
        sample = selected_samples[cohort]
        if (sample not in roster):
            raise ValueError("Archived selected sample is outside its cohort")
        selected[cohort] = {
            "sample_id": sample,
            "complete_roster": roster.copy(),
            "selection": "ARCHIVED_OUTCOME_BLIND_FIXED_SELECTION",
        }
    return selected

def make_panel(records, graphs, selected_samples):
    selected = choose_samples(records, selected_samples)
    cases = []
    for (cohort, transition, direction) in (itertools.product(
        COHORTS, TRANSITIONS, DIRECTIONS
    )):
        matches = [
            graph
            for (graph) in (graphs)
            if (
                (
                    graph.get("cohort"),
                    graph.get("transition"),
                    graph.get("direction"),
                    graph.get("class_id"),
                    graph.get("kind"),
                )
                == (cohort, transition, direction, "ALL", "real")
            )
        ]
        if (len(matches) != 1):
            raise ValueError("Each complete ALL graph unit must match exactly once")
        graph = matches[0]
        expected = "XB__" + "__".join((transition, cohort, direction, "ALL"))
        if (
            graph["network_id"] != expected
            or graph["full_cohort_n"] != 36
            or graph["input_status"] != "READY"
        ):
            raise ValueError("Accepted graph identity or population differs")
        sample = selected[cohort]["sample_id"]
        cases.append(
            {
                "case_id": cohort + "|" + transition + "|" + direction + "|" + sample,
                "cohort": cohort,
                "sample_id": sample,
                "transition": transition,
                "direction": direction,
                "network_id": expected,
                "graph": binding(graph["canonical_edges"]),
                "graph_preparation": binding(graph["source_receipt"]),
                "complete_tumour_count": 36,
                "molecular_scope": "COMPLETE_GRAPH_MATCHED_COHORT_AND_NORMALS_REQUIRED",
                "graph_specific_run_spec": None,
            }
        )
    return {
        "status": "INPUT_ONLY_FIXED_REPRESENTATION_PANEL",
        "samples": selected,
        "cases": cases,
        "case_count": 8,
        "biological_samples": 2,
        "master_seeds": list(SEEDS),
        "representations": list(REPRESENTATIONS),
        "independent_process_repeats": 2,
        "logical_observations": 240,
        "native_observations_max": 240,
        "previous_outcomes_used_for_selection": False,
    }

def block_summary(comparisons):
    labels = [
        (representation, repeat)
        for (representation) in (REPRESENTATIONS)
        for (repeat) in (REPEATS)
    ]
    expected = set(itertools.combinations(labels, 2))
    by_pair = {}
    for (record) in (comparisons):
        pair = (tuple(record["left"]), tuple(record["right"]))
        if (pair not in expected or pair in by_pair):
            raise ValueError("Missing, duplicate or unexpected comparison pair")
        value = record["top10_fractional_L1"]
        if (value is not None and (
            not isinstance(value, (int, float, Fraction))
            or isinstance(value, bool)
            or (not math.isfinite(float(value)))
            or (value < 0)
            or (value > 20)
        )):
            raise ValueError("Invalid top10 distance")
        by_pair[pair] = value
    if (set(by_pair) != expected):
        raise ValueError("All fifteen comparison pairs required")

    def value(left, right):
        return by_pair[
            (left, right)
            if (labels.index(left) < labels.index(right))
            else (right, left)
        ]

    within = {
        representation: value((representation, 1), (representation, 2))
        for (representation) in (REPRESENTATIONS)
    }
    contrasts = {}
    for (representation) in (REPRESENTATIONS[1:]):
        cross = [
            value(("ORIGINAL", i), (representation, j))
            for ((i, j)) in (itertools.product(REPEATS, REPEATS))
        ]
        complete = all(
            (
                item is not None
                for (item) in (cross + [within["ORIGINAL"], within[representation]])
            )
        )
        cross_mean = sum(cross) / 4 if (complete) else None
        repeat_mean = (
            (within["ORIGINAL"] + within[representation]) / 2 if (complete) else None
        )
        contrasts[representation] = {
            "cross_distances": cross,
            "cross_mean": cross_mean,
            "mean_within_repeat": repeat_mean,
            "descriptive_excess": cross_mean - repeat_mean if (complete) else None,
            "status": "FINITE_DESCRIPTIVE_NOT_CAUSAL_OR_INDEPENDENT_INFERENCE"
            if (complete)
            else "UNAVAILABLE_PRESERVED",
        }
    return {"within_representation_repeat_distances": within, "contrasts": contrasts}

def case_summary(blocks):
    if ([block["master_seed"] for (block) in (blocks)] != list(SEEDS)):
        raise ValueError("All five ordered masters required")
    result = {
        "required_seeds": 5,
        "seeds_are_biological_replicates": False,
        "contrasts": {},
    }
    for (representation) in (REPRESENTATIONS[1:]):
        values = [
            block["summary"]["contrasts"][representation]["descriptive_excess"]
            for (block) in (blocks)
        ]
        available = [value for (value) in (values) if (value is not None)]
        result["contrasts"][representation] = {
            "complete_five_seed_mean": sum(available) / 5
            if (len(available) == 5)
            else None,
            "available_seed_count": len(available),
            "all_seed_values": values,
        }
    return result
