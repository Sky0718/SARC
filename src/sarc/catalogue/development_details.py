from ..transfer.scoring import (
    bootstrap_indices,
    leading_membership,
    leading_stability,
    paired_contrast,
    paired_interval,
)
from .constants import FACTORS, SCALES, SEEDS

def membership_records(completed, references):
    records = []
    indexed = {
        (
            item["cell"]["context"],
            item["cell"]["method"],
            item["cell"]["network_id"],
            item["cell"]["parameter_id"],
            item["cell"]["master_seed"],
        ): item
        for (item) in (completed)
    }
    for (context) in (("COAD_CCLE", "LUAD_CCLE")):
        for (method) in (("PersonaDrive", "DawnRank", "PRODIGY")):
            parameter = {
                "PersonaDrive": "original",
                "DawnRank": "mu_3",
                "PRODIGY": "alpha_0.05",
            }[method]
            for (old, new) in ((("11_0", "11_5"), ("11_5", "12_0"))):
                transition = old + "_to_" + new
                network_pairs = [
                    ("native", "native_" + old, "native_" + new),
                    (
                        "persistent",
                        transition + "__persistent_old",
                        transition + "__persistent_new",
                    ),
                ]
                if (method == "PRODIGY"):
                    network_pairs.extend(
                        [
                            (
                                label,
                                transition + "__" + first,
                                transition + "__" + second,
                            )
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
                                    (
                                        "topology_first_step",
                                        "persistent_old",
                                        "topology_first",
                                    ),
                                    (
                                        "confidence_after_topology",
                                        "topology_first",
                                        "persistent_new",
                                    ),
                                )
                            )
                        ]
                    )
                for (seed) in (SEEDS if (method == "PRODIGY") else (None,)):
                    for (label, first, second) in (network_pairs):
                        before = indexed[(context, method, first, parameter, seed)]
                        after = indexed[(context, method, second, parameter, seed)]
                        for (sample, ref) in (references[context].items()):
                            left = before["samples"][sample]
                            right = after["samples"][sample]
                            for (support) in (("common", "native")):
                                eligible = ref[
                                    "common_eligible_literal_gene_labels"
                                    if (support == "common")
                                    else "native_eligible_literal_gene_labels"
                                ]
                                for (k) in (SCALES):
                                    old_membership = leading_membership(
                                        left["genes"], left["scores"], eligible, k
                                    )
                                    new_membership = leading_membership(
                                        right["genes"], right["scores"], eligible, k
                                    )
                                    records.append(
                                        {
                                            "context": context,
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
                                            **leading_stability(
                                                old_membership, new_membership
                                            ),
                                        }
                                    )
    return records

def variability_records(completed, groups):
    sample_records = []
    for (key, group) in (groups.items()):
        if (key[1] != "PRODIGY" or group["expected_seeds"] != SEEDS):
            continue
        for (projection, rows) in (group["measurements"].items()):
            for (sample, result) in (rows.items()):
                sample_records.append(
                    {
                        **dict(zip(FACTORS, key)),
                        "sample_id": sample,
                        "support": projection[0],
                        "reference": projection[1],
                        "k": projection[2],
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
    indexed = {
        (
            item["cell"]["context"],
            item["cell"]["network_id"],
            item["cell"]["master_seed"],
        ): item
        for (item) in (completed)
        if (
            item["cell"]["method"] == "PRODIGY"
            and item["cell"]["parameter_id"] == "alpha_0.05"
        )
    }
    paired_records = []
    for (context) in (("COAD_CCLE", "LUAD_CCLE")):
        for (old, new) in ((
            ("native_11_0", "native_11_5"),
            ("native_11_5", "native_12_0"),
        )):
            for (seed) in (SEEDS):
                before = indexed[(context, old, seed)]
                after = indexed[(context, new, seed)]
                population = before["cell"]["sample_order"]
                indices = bootstrap_indices(len(population))
                projections = before["samples"][population[0]]["measurements"]
                for (projection) in (projections):
                    first = {
                        sample: before["samples"][sample]["measurements"][projection]
                        for (sample) in (population)
                    }
                    second = {
                        sample: after["samples"][sample]["measurements"][projection]
                        for (sample) in (population)
                    }
                    result = paired_contrast(population, [(1, second), (-1, first)])
                    paired_records.append(
                        {
                            "context": context,
                            "method": "PRODIGY",
                            "old_network": old,
                            "new_network": new,
                            "master_seed": seed,
                            "support": projection[0],
                            "reference": projection[1],
                            "k": projection[2],
                            "estimate": result,
                            "interval": paired_interval(result, indices),
                        }
                    )
    return sample_records, paired_records
