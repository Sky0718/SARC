import itertools
from fractions import Fraction
import numpy as np
from ..transfer import scoring as CORE
from .original_prediction import MASTERS, METHODS, primary_event

REFERENCES = (
    "NCG6_primary_all",
    "NCG6_known_subgroup",
    "NCG6_candidate_subgroup",
    "CancerMine2019_secondary",
)
SCALES = (10, 1, 5, 20)
SUPPORTS = ("common", "native")
ENDPOINTS = tuple(
    (
        (support, reference, k)
        for (support) in (SUPPORTS)
        for (reference) in (REFERENCES)
        for (k) in (SCALES)
    )
)
CLASSES = ("PB", "PW", "DB", "DW", "PB_PW", "ALL")
CONTROLS = (271828, 314159, 173205)
COHORTS = {"COAD_CCLE": 36, "LUAD_CCLE": 36, "COAD_TCGA": 396}

def exact(value):
    if (value is None or isinstance(value, Fraction)):
        return value
    if (isinstance(value, dict) and set(value) >= {"numerator", "denominator"}):
        return Fraction(value["numerator"], value["denominator"])
    return Fraction(value)

def encode(value):
    if (isinstance(value, Fraction)):
        return {"numerator": value.numerator, "denominator": value.denominator}
    raise TypeError(type(value).__name__)

def seeds(method):
    if (method not in METHODS):
        raise ValueError("Unknown method")
    return MASTERS if (method == "PRODIGY") else (None,)

def endpoint_id(endpoint):
    return "::".join(map(str, endpoint))

def complete_mean(values):
    return (
        sum(values, Fraction(0)) / len(values)
        if (values and all((value is not None for (value) in (values))))
        else None
    )

def reference_person(person):
    for (support) in (SUPPORTS):
        labels = person[support + "_eligible_literal_gene_labels"]
        positives = person[support + "_reference_positive_labels"]
        if (len(labels) != len(set(labels)) or set(positives) != set(REFERENCES)):
            raise ValueError(
                "Full unique support and four inherited reference sets required"
            )
        if (any(
            (
                len(value) != len(set(value)) or not set(value).issubset(labels)
                for (value) in (positives.values())
            )
        )):
            raise ValueError("Reference positives escape support or contain duplicates")
    if (not set(person["common_eligible_literal_gene_labels"]).issubset(
        person["native_eligible_literal_gene_labels"]
    )):
        raise ValueError("Frozen common support must be a native subset")

def raw(record, person):
    state = record["status"]
    if (state not in ("SUCCESS", "SUCCESS_EMPTY", "FAILED", "ERROR", "UNAVAILABLE")):
        raise ValueError("Unknown preserved scientific outcome state")
    (genes, scores) = (record["genes"], record["scores"])
    if (state not in ("SUCCESS", "SUCCESS_EMPTY")):
        if (genes or scores):
            raise ValueError("Unavailable output contains ranked values")
        return (None, None)
    if ((state == "SUCCESS") != bool(genes)):
        raise ValueError("Successful output violates original empty status")
    CORE.score_list(
        genes, scores, person["native_eligible_literal_gene_labels"], [], 10
    )
    return (genes, scores)

def reduce_measurements(values):
    valid = all((row["hits"] is not None for (row) in (values)))
    return {
        "hits": complete_mean([row["hits"] for (row) in (values)]),
        "lower": complete_mean([row["lower"] for (row) in (values)]),
        "upper": complete_mean([row["upper"] for (row) in (values)]),
        "recall": complete_mean([row["recall"] for (row) in (values)]),
        "native_output_count": complete_mean(
            [row["native_output_count"] for (row) in (values)]
        ),
        "eligible_output_count": complete_mean(
            [row["common_output_count"] for (row) in (values)]
        ),
        "positive_output_count": complete_mean(
            [row["positive_output_count"] for (row) in (values)]
        ),
        "eligible_count": values[0]["eligible_count"],
        "positive_count": values[0]["positive_count"],
        "status": "AVAILABLE" if (valid) else "UNAVAILABLE",
        "seed_hits": [row["hits"] for (row) in (values)],
        "seed_statuses": [row["status"] for (row) in (values)],
    }

def compare_rankings(method, baseline, intervention, person, accepted_primary = None):
    reference_person(person)
    masters = seeds(method)
    if (set(baseline) != set(masters) or set(intervention) != set(masters)):
        raise ValueError("All exact paired algorithm seeds required")
    prepared = {"baseline": {}, "intervention": {}}
    for (name, ranks) in ((("baseline", baseline), ("intervention", intervention))):
        prepared[name] = {seed: raw(ranks[seed], person) for (seed) in (masters)}
    endpoints = {}
    for (endpoint) in (ENDPOINTS):
        (support, reference, k) = endpoint
        eligible = person[support + "_eligible_literal_gene_labels"]
        positives = person[support + "_reference_positive_labels"][reference]
        sides = {
            name: [
                CORE.score_list(*prepared[name][seed], eligible, positives, k)
                for (seed) in (masters)
            ]
            for (name) in (prepared)
        }
        (before, after) = (
            reduce_measurements(sides[name])
            for (name) in (("baseline", "intervention"))
        )
        if (
            endpoint == ("common", "NCG6_primary_all", 10)
            and accepted_primary is not None
        ):
            event = accepted_primary
        else:
            event = primary_event(
                method,
                dict(zip(masters, before["seed_hits"])),
                dict(zip(masters, after["seed_hits"])),
            )
        if (endpoint != ("common", "NCG6_primary_all", 10)):
            event.pop("event", None)
        endpoints[endpoint_id(endpoint)] = {
            "baseline": before,
            "intervention": after,
            "effect": event,
            "support_genes": sorted(eligible),
            "positive_genes": sorted(positives),
        }
    coverage = {}
    for (support) in (SUPPORTS):
        eligible = set(person[support + "_eligible_literal_gene_labels"])
        per_seed = []
        for (seed) in (masters):
            (old, new) = (
                prepared[name][seed] for (name) in (("baseline", "intervention"))
            )
            if (old[0] is None or new[0] is None):
                per_seed.append({"status": "UNAVAILABLE", "master_seed": seed})
                continue
            (old_set, new_set) = (set(pair[0]) & eligible for (pair) in ((old, new)))
            union = old_set | new_set
            top = {}
            for (k) in (SCALES):
                old_top = CORE.leading_membership(*old, eligible, k)
                new_top = CORE.leading_membership(*new, eligible, k)
                stability = CORE.leading_stability(old_top, new_top)
                keys = old_top.keys() | new_top.keys()
                top[str(k)] = stability | {
                    "old_membership": old_top,
                    "new_membership": new_top,
                    "fractional_membership_L1": sum(
                        (
                            abs(old_top.get(gene, 0) - new_top.get(gene, 0))
                            for (gene) in (keys)
                        ),
                        Fraction(0),
                    ),
                    "old_returned_mass": sum(old_top.values(), Fraction(0)),
                    "new_returned_mass": sum(new_top.values(), Fraction(0)),
                    "turnover_fraction": None
                    if (stability["fractional_membership_jaccard"] is None)
                    else 1 - stability["fractional_membership_jaccard"],
                }
            per_seed.append(
                {
                    "status": "AVAILABLE",
                    "master_seed": seed,
                    "eligible_count": len(eligible),
                    "old_candidate_count": len(old_set),
                    "new_candidate_count": len(new_set),
                    "old_coverage": Fraction(len(old_set), len(eligible))
                    if (eligible)
                    else None,
                    "new_coverage": Fraction(len(new_set), len(eligible))
                    if (eligible)
                    else None,
                    "gained_candidates": sorted(new_set - old_set),
                    "lost_candidates": sorted(old_set - new_set),
                    "candidate_jaccard": Fraction(len(old_set & new_set), len(union))
                    if (union)
                    else None,
                    "top": top,
                }
            )
        coverage[support] = {
            "seed_records": per_seed,
            "all_required_seeds_available": all(
                (row["status"] == "AVAILABLE" for (row) in (per_seed))
            ),
            "seeds_are_biological_replicates": False,
        }
    return {
        "method": method,
        "master_seeds": list(masters),
        "endpoints": endpoints,
        "coverage": coverage,
        "primary_event": endpoints["common::NCG6_primary_all::10"]["effect"],
    }

def validate_draws(indices, n, synthetic = False):
    if (
        len(indices) != (len(indices) if (synthetic) else 2000)
        or len(indices) < 2
        or any(
            (
                len(draw) != n
                or any(
                    (
                        type(index) is not int or index < 0 or index >= n
                        for (index) in (draw)
                    )
                )
                for (draw) in (indices)
            )
        )
    ):
        raise ValueError("Complete inherited paired biological-position bank required")
    return np.asarray(
        [np.bincount(draw, minlength = n) for (draw) in (indices)], dtype = np.float64
    )

def scalar_summary(population, values, counts, lower = None, upper = None, inference = True):
    if (
        len(population) != len(set(population))
        or set(values) != set(population)
        or counts.shape[1] != len(population)
    ):
        raise ValueError("Full paired biological denominator required")
    vector = [exact(values[sample]) for (sample) in (population)]
    available = [value for (value) in (vector) if (value is not None)]
    point = complete_mean(vector)
    interval = {
        "status": "UNAVAILABLE_FULL_POPULATION_CONTRAST",
        "lower_95": None,
        "upper_95": None,
        "draws": len(counts),
    }
    if (not inference):
        interval["status"] = "DESCRIPTIVE_SECONDARY_NO_INTERVAL_COMPUTED"
    if (point is not None and inference):
        draws = (
            counts
            @ np.asarray([float(value) for (value) in (vector)])
            / len(population)
        )
        bounds = np.quantile(draws, (0.025, 0.975), method = "linear")
        interval.update(
            {
                "status": "CONDITIONAL_FIXED_GRAPHS_AND_SAMPLE_COMPOSITION",
                "lower_95": float(bounds[0]),
                "upper_95": float(bounds[1]),
            }
        )
    return {
        "biological_total": len(population),
        "biological_available": len(available),
        "missing_sample_ids": [
            sample for ((sample, value)) in (zip(population, vector)) if (value is None)
        ],
        "value": point,
        "complete_subset_value_secondary": sum(available, Fraction(0)) / len(available)
        if (available)
        else None,
        "lower": complete_mean([exact(lower[sample]) for (sample) in (population)])
        if (lower is not None)
        else None,
        "upper": complete_mean([exact(upper[sample]) for (sample) in (population)])
        if (upper is not None)
        else None,
        "paired_interval": interval,
    }

def linear_contrast(population, terms, counts):
    rows = {}
    for (sample) in (population):
        measurements = [
            (coefficient, source[sample]) for ((coefficient, source)) in (terms)
        ]
        rows[sample] = {
            "value": complete_mean([Fraction(0)])
            if (all((row["hits"] is not None for ((coefficient, row)) in (measurements))))
            else None,
            "lower": sum(
                (
                    coefficient * exact(row["lower" if (coefficient > 0) else "upper"])
                    for ((coefficient, row)) in (measurements)
                ),
                Fraction(0),
            ),
            "upper": sum(
                (
                    coefficient * exact(row["upper" if (coefficient > 0) else "lower"])
                    for ((coefficient, row)) in (measurements)
                ),
                Fraction(0),
            ),
        }
        if (rows[sample]["value"] is not None):
            rows[sample]["value"] = sum(
                (
                    coefficient * exact(row["hits"])
                    for ((coefficient, row)) in (measurements)
                ),
                Fraction(0),
            )
    values = {sample: row["value"] for ((sample, row)) in (rows.items())}
    lower = {sample: row["lower"] for ((sample, row)) in (rows.items())}
    upper = {sample: row["upper"] for ((sample, row)) in (rows.items())}
    absolute = {
        sample: abs(value) if (value is not None) else None
        for ((sample, value)) in (values.items())
    }
    absolute_lower = {
        sample: Fraction(0)
        if (lower[sample] <= 0 <= upper[sample])
        else min(abs(lower[sample]), abs(upper[sample]))
        for (sample) in (population)
    }
    absolute_upper = {
        sample: max(abs(lower[sample]), abs(upper[sample])) for (sample) in (population)
    }
    return {
        "signed": scalar_summary(population, values, counts, lower, upper),
        "absolute": scalar_summary(
            population, absolute, counts, absolute_lower, absolute_upper
        ),
        "samples": rows,
    }

def panel_summaries(population, graphs, records, counts):
    by_slot = {
        (graph["kind"], graph["class_id"], graph["control_seed"]): graph
        for (graph) in (graphs)
    }
    output = []
    for (endpoint) in (ENDPOINTS):
        ep = endpoint_id(endpoint)
        for (method) in (METHODS):

            def side(slot, which):
                graph = by_slot[slot]
                return {
                    sample: records[graph["network_id"], method, sample]["endpoints"][
                        ep
                    ][which]
                    for (sample) in (population)
                }

            for (slot, graph) in (by_slot.items()):
                (before, after) = (side(slot, "baseline"), side(slot, "intervention"))
                output.append(
                    {
                        "analysis": "class_effect",
                        "endpoint": list(endpoint),
                        "method": method,
                        "network_id": graph["network_id"],
                        "kind": slot[0],
                        "class_id": slot[1],
                        "control_seed": slot[2],
                        "result": linear_contrast(
                            population, [(1, after), (-1, before)], counts
                        ),
                    }
                )
                levels = {
                    name: {
                        metric: scalar_summary(
                            population,
                            {sample: row[metric] for ((sample, row)) in (values.items())},
                            counts,
                            inference = False,
                        )
                        for (metric) in (
                            (
                                "hits",
                                "recall",
                                "native_output_count",
                                "eligible_output_count",
                                "positive_output_count",
                            )
                        )
                    }
                    for ((name, values)) in (
                        (
                            ("baseline", before),
                            ("intervention", after),
                        )
                    )
                }
                output.append(
                    {
                        "analysis": "inherited_endpoint_levels",
                        "endpoint": list(endpoint),
                        "method": method,
                        "network_id": graph["network_id"],
                        "levels": levels,
                    }
                )
            for (kind, seed) in ((
                ("real", None),
                *(("matched_control", value) for (value) in (CONTROLS)),
            )):
                terms = [
                    (1, side((kind, "PB_PW", seed), "intervention")),
                    (-1, side((kind, "PB", seed), "intervention")),
                    (-1, side((kind, "PW", seed), "intervention")),
                    (1, side((kind, "PB", seed), "baseline")),
                ]
                output.append(
                    {
                        "analysis": "PB_PW_nonadditivity",
                        "endpoint": list(endpoint),
                        "method": method,
                        "kind": kind,
                        "control_seed": seed,
                        "definition": "Q(PB_PW)-Q(PB)-Q(PW)+Q(baseline)",
                        "result": linear_contrast(population, terms, counts),
                    }
                )
            for (seed) in (CONTROLS):
                interaction_terms = []
                for (kind, graph_seed, factor) in ((
                    ("real", None, 1),
                    ("matched_control", seed, -1),
                )):
                    interaction_terms.extend(
                        [
                            (factor, side((kind, "PB_PW", graph_seed), "intervention")),
                            (-factor, side((kind, "PB", graph_seed), "intervention")),
                            (-factor, side((kind, "PW", graph_seed), "intervention")),
                            (factor, side((kind, "PB", graph_seed), "baseline")),
                        ]
                    )
                output.append(
                    {
                        "analysis": "real_minus_control_PB_PW_nonadditivity",
                        "endpoint": list(endpoint),
                        "method": method,
                        "control_seed": seed,
                        "result": linear_contrast(
                            population, interaction_terms, counts
                        ),
                        "finite_control_not_biological_replication": True,
                    }
                )
            for (name) in (CLASSES[:-1]):
                real_slot = ("real", name, None)
                for (seed) in (CONTROLS):
                    control_slot = ("matched_control", name, seed)
                    terms = [
                        (1, side(real_slot, "intervention")),
                        (-1, side(real_slot, "baseline")),
                        (-1, side(control_slot, "intervention")),
                        (1, side(control_slot, "baseline")),
                    ]
                    real_effect = linear_contrast(population, terms[:2], counts)
                    control_effect = linear_contrast(
                        population,
                        [
                            (-coefficient, values)
                            for ((coefficient, values)) in (terms[2:])
                        ],
                        counts,
                    )
                    absolute_difference = {
                        sample: abs(real_effect["samples"][sample]["value"])
                        - abs(control_effect["samples"][sample]["value"])
                        if (
                            real_effect["samples"][sample]["value"] is not None
                            and control_effect["samples"][sample]["value"] is not None
                        )
                        else None
                        for (sample) in (population)
                    }
                    output.append(
                        {
                            "analysis": "real_minus_matched_control",
                            "endpoint": list(endpoint),
                            "method": method,
                            "class_id": name,
                            "control_seed": seed,
                            "signed_difference": linear_contrast(
                                population, terms, counts
                            ),
                            "difference_of_absolute_effects": scalar_summary(
                                population, absolute_difference, counts
                            ),
                            "randomisation_p_value": None,
                            "finite_control_not_biological_replication": True,
                        }
                    )
        if (endpoint[0] == "common"):
            for (graph) in (graphs):
                for (first, second) in (itertools.combinations(METHODS, 2)):
                    cells = {
                        method: {
                            sample: records[graph["network_id"], method, sample][
                                "endpoints"
                            ][ep]
                            for (sample) in (population)
                        }
                        for (method) in ((first, second))
                    }
                    if (any(
                        (
                            cells[first][sample]["support_genes"]
                            != cells[second][sample]["support_genes"]
                            or cells[first][sample]["positive_genes"]
                            != cells[second][sample]["positive_genes"]
                            for (sample) in (population)
                        )
                    )):
                        raise ValueError(
                            "Method-gap endpoints do not share bound eligibility and positive labels"
                        )
                    terms = [
                        (
                            1,
                            {
                                sample: cells[first][sample]["intervention"]
                                for (sample) in (population)
                            },
                        ),
                        (
                            -1,
                            {
                                sample: cells[first][sample]["baseline"]
                                for (sample) in (population)
                            },
                        ),
                        (
                            -1,
                            {
                                sample: cells[second][sample]["intervention"]
                                for (sample) in (population)
                            },
                        ),
                        (
                            1,
                            {
                                sample: cells[second][sample]["baseline"]
                                for (sample) in (population)
                            },
                        ),
                    ]
                    output.append(
                        {
                            "analysis": "method_gap_change",
                            "endpoint": list(endpoint),
                            "methods": [first, second],
                            "network_id": graph["network_id"],
                            "result": linear_contrast(population, terms, counts),
                        }
                    )
    for (graph) in (graphs):
        for (method) in (METHODS):
            for (support) in (SUPPORTS):
                measurements = {}
                for (sample) in (population):
                    seed_rows = records[graph["network_id"], method, sample][
                        "coverage"
                    ][support]["seed_records"]
                    sample_values = {}
                    for (name) in ((
                        "old_candidate_count",
                        "new_candidate_count",
                        "old_coverage",
                        "new_coverage",
                        "candidate_jaccard",
                    )):
                        sample_values[name] = complete_mean(
                            [row.get(name) for (row) in (seed_rows)]
                        )
                    for (k) in (SCALES):
                        for (name) in ((
                            "fractional_membership_jaccard",
                            "fractional_membership_L1",
                            "old_returned_mass",
                            "new_returned_mass",
                            "turnover_fraction",
                        )):
                            sample_values[str(k) + "::" + name] = complete_mean(
                                [
                                    row.get("top", {}).get(str(k), {}).get(name)
                                    for (row) in (seed_rows)
                                ]
                            )
                    measurements[sample] = sample_values
                panels = {
                    name: scalar_summary(
                        population,
                        {
                            sample: values[name]
                            for ((sample, values)) in (measurements.items())
                        },
                        counts,
                        inference = False,
                    )
                    for (name) in (measurements[population[0]])
                }
                output.append(
                    {
                        "analysis": "candidate_coverage_and_membership",
                        "network_id": graph["network_id"],
                        "method": method,
                        "support": support,
                        "panels": panels,
                        "exact_seed_memberships_preserved_in_pieces": True,
                    }
                )
    return output
