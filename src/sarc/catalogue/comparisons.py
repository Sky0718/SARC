from ..transfer import scoring as core
from .constants import (
    CONTEXT,
    ENDPOINTS,
    METHODS,
    NATIVE,
    SEEDS,
    TRANSITIONS,
    group_key,
)

CONTRAST_METHODS = ("PersonaDrive", "DawnRank", "PRODIGY")

def contrast_record(groups, name, terms, endpoint, bootstrap):
    selected = [(coefficient, groups[key]) for ((coefficient, key)) in (terms)]
    population = selected[0][1]["population"]
    if (any(
        group["population"] != population or key[0] != CONTEXT
        for (((coefficient, key), (unused, group))) in (zip(terms, selected))
    )):
        raise ValueError("Contrast populations differ")
    result = core.paired_contrast(
        population,
        [
            (coefficient, group["measurements"][endpoint])
            for ((coefficient, group)) in (selected)
        ],
    )
    return {
        "context": CONTEXT,
        "contrast": name,
        "support": endpoint[0],
        "reference": endpoint[1],
        "k": endpoint[2],
        "terms": [
            {"coefficient": coefficient, "group": list(key)}
            for ((coefficient, key)) in (terms)
        ],
        "estimate": result,
        "interval": bootstrap.contrast_interval(result),
    }

def comparisons(groups, bootstrap):
    records, reversals = [], []
    for (endpoint) in (ENDPOINTS):
        gaps = {}
        for (network) in (NATIVE):
            release = network.removeprefix("native_")
            for (index, first) in (enumerate(CONTRAST_METHODS)):
                for (second) in (CONTRAST_METHODS[index + 1 :]):
                    terms = [
                        (1, group_key(first, network)),
                        (-1, group_key(second, network)),
                    ]
                    record = contrast_record(
                        groups,
                        "method_gap::" + first + "::" + second + "::" + release,
                        terms,
                        endpoint,
                        bootstrap,
                    )
                    gaps[(first, second, release)] = record
                    records.append(record)
        for (transition) in (TRANSITIONS):
            old, new = transition.split("_to_")
            for (method) in (CONTRAST_METHODS):
                n0, n1 = (
                    group_key(method, "native_" + old),
                    group_key(method, "native_" + new),
                )
                p0, p1 = (
                    group_key(method, transition + "__persistent_old"),
                    group_key(method, transition + "__persistent_new"),
                )
                definitions = {
                    "native_release": [(1, n1), (-1, n0)],
                    "persistent_release": [(1, p1), (-1, p0)],
                    "old_support_restriction": [(1, p0), (-1, n0)],
                    "new_support_restriction": [(1, p1), (-1, n1)],
                    "native_minus_persistent_release": [
                        (1, n1),
                        (-1, n0),
                        (-1, p1),
                        (1, p0),
                    ],
                }
                if (method == "PRODIGY"):
                    confidence = group_key(method, transition + "__confidence_first")
                    topology = group_key(method, transition + "__topology_first")
                    definitions.update(
                        {
                            "confidence_then_topology_first_step": [
                                (1, confidence),
                                (-1, p0),
                            ],
                            "confidence_then_topology_second_step": [
                                (1, p1),
                                (-1, confidence),
                            ],
                            "topology_then_confidence_first_step": [
                                (1, topology),
                                (-1, p0),
                            ],
                            "topology_then_confidence_second_step": [
                                (1, p1),
                                (-1, topology),
                            ],
                            "confidence_topology_nonadditivity": [
                                (1, p1),
                                (-1, confidence),
                                (-1, topology),
                                (1, p0),
                            ],
                        }
                    )
                for (name, terms) in (definitions.items()):
                    records.append(
                        contrast_record(
                            groups,
                            name + "::" + method + "::" + transition,
                            terms,
                            endpoint,
                            bootstrap,
                        )
                    )
            for (index, first) in (enumerate(CONTRAST_METHODS)):
                for (second) in (CONTRAST_METHODS[index + 1 :]):
                    before, after = (
                        gaps[(first, second, old)],
                        gaps[(first, second, new)],
                    )
                    terms = [
                        (1, group_key(first, "native_" + new)),
                        (-1, group_key(second, "native_" + new)),
                        (-1, group_key(first, "native_" + old)),
                        (1, group_key(second, "native_" + old)),
                    ]
                    records.append(
                        contrast_record(
                            groups,
                            "interaction::" + first + "::" + second + "::" + transition,
                            terms,
                            endpoint,
                            bootstrap,
                        )
                    )
                    for (delta) in ((1, 0.5, 2)):
                        reversals.append(
                            {
                                "context": CONTEXT,
                                "first_method": first,
                                "second_method": second,
                                "transition": transition,
                                "support": endpoint[0],
                                "reference": endpoint[1],
                                "k": endpoint[2],
                                **core.reversals(
                                    before["estimate"]["value"],
                                    after["estimate"]["value"],
                                    before["interval"],
                                    after["interval"],
                                    delta,
                                ),
                            }
                        )
    if (len(records) != 1760 or len(reversals) != 576):
        raise ValueError("Complete fixed endpoint contrast family does not close")
    return (records, reversals)

def membership_records(completed, references):
    indexed = {
        (
            value["cell"]["method"],
            value["cell"]["network_id"],
            value["cell"]["master_seed"],
        ): value
        for (value) in (completed.values())
    }
    for (method) in (METHODS):
        for (transition) in (TRANSITIONS):
            old, new = transition.split("_to_")
            pairs = [
                ("native", "native_" + old, "native_" + new),
                (
                    "persistent",
                    transition + "__persistent_old",
                    transition + "__persistent_new",
                ),
            ]
            if (method == "PRODIGY"):
                pairs.extend(
                    (label, transition + "__" + first, transition + "__" + second)
                    for ((label, first, second)) in (
                        (
                            (
                                "confidence_first_step",
                                "persistent_old",
                                "confidence_first",
                            ),
                            (
                                "topology_after_confidence",
                                "confidence_first",
                                "persistent_new",
                            ),
                            ("topology_first_step", "persistent_old", "topology_first"),
                            (
                                "confidence_after_topology",
                                "topology_first",
                                "persistent_new",
                            ),
                        )
                    )
                )
            for (seed) in (SEEDS if (method == "PRODIGY") else (None,)):
                for (label, first, second) in (pairs):
                    before, after = (
                        indexed[(method, first, seed)],
                        indexed[(method, second, seed)],
                    )
                    for (sample, source) in (references.items()):
                        left, right = (
                            before["samples"][sample],
                            after["samples"][sample],
                        )
                        for (support) in (("common", "native")):
                            eligible = source[support + "_eligible_literal_gene_labels"]
                            for (k) in ((10, 1, 5, 20)):
                                old_membership = core.leading_membership(
                                    left["genes"], left["scores"], eligible, k
                                )
                                new_membership = core.leading_membership(
                                    right["genes"], right["scores"], eligible, k
                                )
                                yield (
                                    {
                                        "context": CONTEXT,
                                        "method": method,
                                        "sample_id": sample,
                                        "transition": transition,
                                        "network_contrast": label,
                                        "master_seed": seed,
                                        "support": support,
                                        "k": k,
                                        "old_member_count": len(old_membership)
                                        if (old_membership is not None)
                                        else None,
                                        "new_member_count": len(new_membership)
                                        if (new_membership is not None)
                                        else None,
                                        **core.leading_stability(
                                            old_membership, new_membership
                                        ),
                                    }
                                )

def variability_records(completed, groups, bootstrap):
    for (key, group) in (groups.items()):
        if (key[1] != "PRODIGY"):
            continue
        for (endpoint, rows) in (group["measurements"].items()):
            for (sample, result) in (rows.items()):
                yield (
                    {
                        "context": CONTEXT,
                        "method": "PRODIGY",
                        "network_id": key[2],
                        "parameter_id": key[3],
                        "sample_id": sample,
                        "support": endpoint[0],
                        "reference": endpoint[1],
                        "k": endpoint[2],
                        "seeds": SEEDS,
                        "completed_repeats": result["complete_count"],
                        "mean": result["hits"],
                        "minimum": result["repeat_minimum"],
                        "maximum": result["repeat_maximum"],
                        "range": result["repeat_maximum"] - result["repeat_minimum"]
                        if (result["repeat_maximum"] is not None)
                        else None,
                    }
                )

def seeded_contrasts(completed, bootstrap):
    indexed = {
        (value["cell"]["network_id"], value["cell"]["master_seed"]): value
        for (value) in (completed.values())
        if (value["cell"]["method"] == "PRODIGY")
    }
    for (transition) in (TRANSITIONS):
        old, new = ("native_" + part for (part) in (transition.split("_to_")))
        for (seed) in (SEEDS):
            before, after = indexed[(old, seed)], indexed[(new, seed)]
            population = before["cell"]["sample_order"]
            for (endpoint) in (ENDPOINTS):
                first = {
                    sample: before["samples"][sample]["measurements"][endpoint]
                    for (sample) in (population)
                }
                second = {
                    sample: after["samples"][sample]["measurements"][endpoint]
                    for (sample) in (population)
                }
                result = core.paired_contrast(population, [(1, second), (-1, first)])
                yield (
                    {
                        "context": CONTEXT,
                        "method": "PRODIGY",
                        "old_network": old,
                        "new_network": new,
                        "master_seed": seed,
                        "support": endpoint[0],
                        "reference": endpoint[1],
                        "k": endpoint[2],
                        "estimate": result,
                        "interval": bootstrap.contrast_interval(result),
                    }
                )
