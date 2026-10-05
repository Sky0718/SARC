from ..transfer.scoring import (
    bootstrap_indices,
    mean_records,
    paired_contrast,
    paired_interval,
    repeat_summary,
    reversals,
)
from .constants import REFERENCES, SCALES, SEEDS

def group_cells(completed):
    groups = {}
    for (value) in (completed):
        cell = value["cell"]
        key = tuple(
            cell[name]
            for (name) in (("context", "method", "network_id", "parameter_id"))
        )
        group = groups.setdefault(key, {})
        seed = cell["master_seed"]
        if (seed in group):
            raise ValueError("Duplicate logical method seed")
        group[seed] = value
    reduced = {}
    for (key, repeats) in (groups.items()):
        first = next(iter(repeats.values()))
        cell = first["cell"]
        expected = (
            SEEDS
            if (cell["method"] == "PRODIGY" and cell["parameter_id"] == "alpha_0.05")
            else (104729,)
            if (cell["method"] == "PRODIGY")
            else (None,)
        )
        if (set(repeats) != set(expected)):
            raise ValueError("Incomplete or unexpected stochastic seed identities")
        population = cell["sample_order"]
        if (any(
            value["cell"]["sample_order"] != population
            for (value) in (repeats.values())
        )):
            raise ValueError("Repeat population mismatch")
        measures = {}
        for (support) in (("common", "native")):
            for (reference) in (REFERENCES):
                for (k) in (SCALES):
                    projection = (support, reference, k)
                    records = {}
                    for (sample) in (population):
                        rows = [
                            repeats[seed]["samples"][sample]["measurements"][projection]
                            for (seed) in (expected)
                        ]
                        if (
                            len(
                                {
                                    (
                                        row["eligible_count"],
                                        row["positive_count"],
                                        row["k"],
                                    )
                                    for (row) in (rows)
                                }
                            )
                            != 1
                        ):
                            raise ValueError("Repeat evaluation support changed")
                        if (len(expected) == 5):
                            result = repeat_summary(rows)
                            positives = rows[0]["positive_count"]
                            records[sample] = result | {
                                "eligible_count": rows[0]["eligible_count"],
                                "positive_count": positives,
                                "k": k,
                                "recall": result["hits"] / positives
                                if (result["hits"] is not None and positives)
                                else None,
                            }
                        else:
                            records[sample] = rows[0]
                    measures[projection] = records
        reduced[key] = {
            "population": population,
            "expected_seeds": expected,
            "measurements": measures,
            "cohort_summaries": {
                projection: mean_records(list(rows.values()))
                for ((projection, rows)) in (measures.items())
            },
        }
    return reduced

def contrast_record(groups, context, name, terms, projection, indices):
    selected = [(coefficient, groups[key]) for ((coefficient, key)) in (terms)]
    population = selected[0][1]["population"]
    if (any(
        group["population"] != population or key[0] != context
        for (((coefficient, key), (unused, group))) in (zip(terms, selected))
    )):
        raise ValueError("Contrasts cannot mix cohorts or sample populations")
    values = [
        (coefficient, group["measurements"][projection])
        for ((coefficient, group)) in (selected)
    ]
    result = paired_contrast(population, values)
    return {
        "context": context,
        "contrast": name,
        "support": projection[0],
        "reference": projection[1],
        "k": projection[2],
        "terms": [
            {"coefficient": coefficient, "group": list(key)}
            for ((coefficient, key)) in (terms)
        ],
        "estimate": result,
        "interval": paired_interval(result, indices),
    }

def primary_key(context, method, network):
    parameter = {
        "PersonaDrive": "original",
        "DawnRank": "mu_3",
        "PRODIGY": "alpha_0.05",
    }[method]
    return (context, method, network, parameter)

def comparisons(groups):
    methods = ("PersonaDrive", "DawnRank", "PRODIGY")
    transitions = (("11_0", "11_5"), ("11_5", "12_0"))
    records = []
    reversal_records = []
    for (context) in (("COAD_CCLE", "LUAD_CCLE")):
        population = groups[primary_key(context, methods[0], "native_11_0")][
            "population"
        ]
        indices = bootstrap_indices(len(population))
        for (support) in (("common", "native")):
            for (reference) in (REFERENCES):
                for (k) in (SCALES):
                    projection = (support, reference, k)
                    gap_records = {}
                    for (release) in (("11_0", "11_5", "12_0")):
                        for (index, first) in (enumerate(methods)):
                            for (second) in (methods[index + 1 :]):
                                terms = [
                                    (
                                        1,
                                        primary_key(
                                            context, first, "native_" + release
                                        ),
                                    ),
                                    (
                                        -1,
                                        primary_key(
                                            context, second, "native_" + release
                                        ),
                                    ),
                                ]
                                record = contrast_record(
                                    groups,
                                    context,
                                    "method_gap::"
                                    + first
                                    + "::"
                                    + second
                                    + "::"
                                    + release,
                                    terms,
                                    projection,
                                    indices,
                                )
                                gap_records[(first, second, release)] = record
                                records.append(record)
                    for (old, new) in (transitions):
                        transition = old + "_to_" + new
                        for (method) in (methods):
                            n0 = primary_key(context, method, "native_" + old)
                            n1 = primary_key(context, method, "native_" + new)
                            p0 = primary_key(
                                context, method, transition + "__persistent_old"
                            )
                            p1 = primary_key(
                                context, method, transition + "__persistent_new"
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
                                confidence = primary_key(
                                    context, method, transition + "__confidence_first"
                                )
                                topology = primary_key(
                                    context, method, transition + "__topology_first"
                                )
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
                                        context,
                                        name + "::" + method + "::" + transition,
                                        terms,
                                        projection,
                                        indices,
                                    )
                                )
                        for (index, first) in (enumerate(methods)):
                            for (second) in (methods[index + 1 :]):
                                old_gap = gap_records[(first, second, old)]
                                new_gap = gap_records[(first, second, new)]
                                terms = [
                                    (1, primary_key(context, first, "native_" + new)),
                                    (-1, primary_key(context, second, "native_" + new)),
                                    (-1, primary_key(context, first, "native_" + old)),
                                    (1, primary_key(context, second, "native_" + old)),
                                ]
                                records.append(
                                    contrast_record(
                                        groups,
                                        context,
                                        "interaction::"
                                        + first
                                        + "::"
                                        + second
                                        + "::"
                                        + transition,
                                        terms,
                                        projection,
                                        indices,
                                    )
                                )
                                for (delta) in ((1, 0.5, 2)):
                                    reversal_records.append(
                                        {
                                            "context": context,
                                            "first_method": first,
                                            "second_method": second,
                                            "transition": transition,
                                            "support": support,
                                            "reference": reference,
                                            "k": k,
                                            **reversals(
                                                old_gap["estimate"]["value"],
                                                new_gap["estimate"]["value"],
                                                old_gap["interval"],
                                                new_gap["interval"],
                                                delta,
                                            ),
                                        }
                                    )
    return records, reversal_records

def recall_summary(rows):
    positive = [row for (row) in (rows.values()) if (row["positive_count"] > 0)]
    if (not positive):
        return {
            "full_population": len(rows),
            "positive_eligible_population": 0,
            "value": None,
            "status": "NO_REFERENCE_POSITIVES",
        }
    normalised = [
        {
            "hits": row["recall"],
            "lower": row["lower"] / row["positive_count"],
            "upper": row["upper"] / row["positive_count"],
        }
        for (row) in (positive)
    ]
    return {
        "full_population": len(rows),
        "positive_eligible_population": len(positive),
        "status": "FIXED_POSITIVE_ELIGIBILITY_SUBSET",
        **mean_records(normalised),
    }

def tuning_trials(completed):
    outputs = {item["cell"]["cell_id"]: item for (item) in (completed)}
    records = []
    for (method, parameters, default) in ((
        ("DawnRank", ("mu_1", "mu_3", "mu_10"), "mu_3"),
        ("PRODIGY", ("alpha_0.01", "alpha_0.05", "alpha_0.1"), "alpha_0.05"),
    )):
        for (release) in (("native_11_0", "native_11_5", "native_12_0")):
            trials = []
            for (parameter) in (parameters):
                cohort_values = []
                for (context) in (("COAD_CCLE", "LUAD_CCLE")):
                    selected = [
                        item
                        for (item) in (outputs.values())
                        if (
                            item["cell"]["context"] == context
                            and item["cell"]["method"] == method
                            and item["cell"]["network_id"] == release
                            and item["cell"]["parameter_id"] == parameter
                            and item["cell"]["master_seed"]
                            == (104729 if (method == "PRODIGY") else None)
                        )
                    ]
                    if (len(selected) != 1):
                        raise ValueError("Retuning trial is absent or duplicated")
                    item = selected[0]
                    rows = [
                        item["samples"][sample]["measurements"][
                            ("common", "NCG6_primary_all", 10)
                        ]
                        for (sample) in (item["cell"]["sample_order"])
                    ]
                    cohort_values.append(
                        {
                            "context": context,
                            "cell_id": item["cell"]["cell_id"],
                            **mean_records(rows),
                        }
                    )
                score = (
                    sum(item["value"] for (item) in (cohort_values)) / 2
                    if (all(item["value"] is not None for (item) in (cohort_values)))
                    else None
                )
                trials.append(
                    {
                        "parameter": parameter,
                        "equal_context_mean": score,
                        "cohorts": cohort_values,
                        "lower": sum(item["lower"] for (item) in (cohort_values)) / 2,
                        "upper": sum(item["upper"] for (item) in (cohort_values)) / 2,
                    }
                )
            complete = all(
                item["equal_context_mean"] is not None for (item) in (trials)
            )
            chosen = (
                sorted(
                    trials,
                    key = lambda item: (
                        -item["equal_context_mean"],
                        item["parameter"] != default,
                        item["parameter"],
                    ),
                )[0]["parameter"]
                if (complete)
                else None
            )
            records.append(
                {
                    "method": method,
                    "release": release,
                    "trial_count": 3,
                    "trials": trials,
                    "selected_parameter": chosen,
                    "status": "COMPLETE_DEVELOPMENT_SECONDARY_ADAPTATION"
                    if (complete)
                    else "UNAVAILABLE_COMPLETE_GRID_COMPARISON",
                    "protected_tuning": False,
                    "selection_criterion": "Equal-context mean common-support NCG6 hits@10; seed104729 only for PRODIGY; original parameter first on exact ties",
                }
            )
    return records
