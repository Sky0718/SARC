from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from ..precision.io import (
    allocate,
    binary_matrix,
    exact_score,
    read_axes,
    read_csv,
    read_json,
    require,
    write_csv,
    write_json,
)
from .coefficients import core_sets, state_name, state_names
from .functional import unpack_allocations
from .signed_mass import statistics, summarise

def fields():
    return (
        "R00",
        "R10",
        "R01",
        "R11",
        "Radd",
        "I_total",
        "I_score_reference",
        "I_decision_reference",
        "P_at_D0",
        "P_at_D1",
        "D_at_P0",
        "D_at_P1",
    )

def score_categories():
    return ("ALL", "BOTH_CORE", "NOT_BOTH_CORE")

def amounts(weights, shared, category):
    return sum(
        (
            weight
            for (gene, weight) in (weights.items())
            if (category == "ALL" or ((gene in shared) == (category == "BOTH_CORE")))
        ),
        Fraction(0),
    )

def decomposition(values):
    old, propagation, damping, both, additive = values
    return dict(
        zip(
            fields(),
            (
                old,
                propagation,
                damping,
                both,
                additive,
                both - propagation - damping + old,
                both - additive,
                additive - propagation - damping + old,
                propagation - old,
                both - damping,
                damping - old,
                both - propagation,
            ),
        )
    )

def intake(directory, category_document):
    directory = Path(directory)
    manifest = read_json(directory / "scores.json")
    require(
        manifest["mode"] == "components" and manifest["epsilon"] == 1e-12,
        "Strict component score source required",
    )
    require(set(manifest["states"]) == set(state_names()), "Incomplete score states")
    document = read_json(directory / "allocations.json")
    roster, allocations = unpack_allocations(document)
    genes, samples = read_axes(directory)
    require(
        len(genes) == 7399 and len(samples) == 85, "Complete first-transition domain required"
    )
    require(
        [row["sample_id"] for (row) in (roster)] == samples,
        "Score and allocation sample axes differ",
    )
    query = binary_matrix(directory / "mutation.bin", len(genes), len(samples))
    require(np.isin(query, [0, 1]).all(), "Nonbinary query")
    errors = {}
    for (row) in (read_csv(directory / "cutoff_certificates.csv")):
        if (int(row["budget"]) == 10):
            key = row["state"], row["sample_id"]
            require(
                key not in errors and row["certified"] == "TRUE",
                "Duplicate or uncertified original cutoff",
            )
            errors[key] = Fraction(row["error_bound"]) * Fraction(1000000000001, 1000000000000)
    require(set(errors) == set(allocations), "Certificate and allocation identities differ")
    scores = {
        state: binary_matrix(directory / (state + "_scores.bin"), len(genes), len(samples))
        for (state) in (state_names())
    }
    require(all((values >= 0).all() for (values) in (scores.values())), "Negative native score")
    organoid, cell_line = core_sets(category_document)
    shared = organoid & cell_line
    return genes, samples, roster, query, allocations, errors, scores, shared

def model_reference(index, identity, genes, query, accepted, errors, scores, shared, h):
    sample = identity["sample_id"]
    pairs = ((0, 0), (1, 0), (0, 1), (1, 1))
    names = [state_name((a, b, h)) for (a, b) in (pairs)]
    vectors = [[exact_score(value) for (value) in (scores[name][:, index])] for (name) in (names)]
    reference = [vectors[1][i] + vectors[2][i] - vectors[0][i] for (i) in (range(len(genes)))]
    interaction = [vectors[3][i] - reference[i] for (i) in (range(len(genes)))]
    positions = np.flatnonzero(query[:, index])
    state_errors = [errors[(name, sample)] for (name) in (names)]
    reference_error = sum(state_errors[:3], Fraction(0))
    interaction_error = sum(state_errors, Fraction(0))
    decisions, thresholds, cutoff_rows, allocation_rows = [], [], [], []
    for (state_index, vector) in (enumerate(vectors + [reference])):
        selected = {genes[position]: vector[position] for (position) in (positions)}
        weights, gap, threshold = allocate(selected, 10)
        decisions.append(weights)
        thresholds.append(threshold)
        name = names[state_index] if (state_index < 4) else f"ADDITIVE_H{h}"
        error = state_errors[state_index] if (state_index < 4) else reference_error
        certified = gap is None or gap > 2 * error
        if (state_index < 4):
            require(
                weights == accepted[(name, sample)] and certified,
                "Native allocation or cutoff mismatch",
            )
        cutoff_rows.append(
            {
                "sample_id": sample,
                "donor_id": identity["donor_id"],
                "h": h,
                "state": name,
                "query_n": len(selected),
                "cutoff_gap": gap,
                "error_bound": error,
                "gap_over_twice_error": float(gap / (2 * error))
                if (gap is not None and error > 0)
                else None,
                "complete_menu": gap is None,
                "certified": certified,
                "threshold": threshold,
            }
        )
        allocation_rows.append(
            {
                **identity,
                "h": h,
                "state": name,
                "weights": weights,
                "numerically_certified": certified,
            }
        )
    displacement = (
        sum(
            (
                abs(decisions[3].get(gene, Fraction(0)) - decisions[4].get(gene, Fraction(0)))
                for (gene) in (set(decisions[3]) | set(decisions[4]))
            ),
            Fraction(0),
        )
        / 2
    )
    summary = {
        **identity,
        "h": h,
        "query_n": len(positions),
        "full_n": len(genes),
        "score_interaction_error_bound": interaction_error,
        "full_score_L1": sum(map(abs, interaction), Fraction(0)),
        "full_score_Linf": max(map(abs, interaction)),
        "query_score_L1": sum((abs(interaction[i]) for (i) in (positions)), Fraction(0)),
        "query_score_Linf": max((abs(interaction[i]) for (i) in (positions)), default = Fraction(0)),
        "full_interaction_certified_nonzero_n": sum(
            abs(value) > interaction_error for (value) in (interaction)
        ),
        "query_interaction_certified_nonzero_n": sum(
            abs(interaction[i]) > interaction_error for (i) in (positions)
        ),
        "reference_negative_full_n": sum(value < 0 for (value) in (reference)),
        "reference_negative_query_n": sum(reference[i] < 0 for (i) in (positions)),
        "reference_score_sum": sum(reference, Fraction(0)),
        "actual_reference_half_L1": displacement,
        "actual_reference_changed": displacement > 0,
        "reference_cutoff_certified": cutoff_rows[-1]["certified"],
        "reference_cutoff_gap": cutoff_rows[-1]["cutoff_gap"],
        "reference_error_bound": reference_error,
    }
    query_rows, pattern_rows, model_rows = [], [], []
    patterns = ["".join(str(value) for (value) in (values)) for (values) in (product((0, 1), repeat = 4))]
    bins = {
        (category, pattern): {"gene_count": 0, "interaction_contribution": Fraction(0)}
        for (category) in (score_categories())
        for (pattern) in (patterns)
    }
    for (position) in (positions):
        gene = genes[position]
        weights = [decision.get(gene, Fraction(0)) for (decision) in (decisions)]
        pattern = (
            "".join(str(int(value)) for (value) in (weights[:4]))
            if (all(value in (0, 1) for (value) in (weights[:4])))
            else "(" + ",".join(str(value) for (value) in (weights[:4])) + ")"
        )
        item = {
            **identity,
            "h": h,
            "gene": gene,
            "shared_core": gene in shared,
            "pattern_00_10_01_11": pattern,
            "w00": weights[0],
            "w10": weights[1],
            "w01": weights[2],
            "w11": weights[3],
            "w_add": weights[4],
            "allocation_interaction": weights[3] - weights[1] - weights[2] + weights[0],
            "score_reference_term": weights[3] - weights[4],
            "decision_reference_term": weights[4] - weights[1] - weights[2] + weights[0],
            "score_interaction": interaction[position],
            "score_interaction_error_bound": interaction_error,
            "score_sign_certified": abs(interaction[position]) > interaction_error,
        }
        for (j, label) in (enumerate(("00", "10", "01", "11", "add"))):
            vector = vectors[j] if (j < 4) else reference
            item["score_" + label] = vector[position]
            item["cutoff_margin_" + label] = (
                vector[position] - thresholds[j] if (thresholds[j] is not None) else None
            )
        query_rows.append(item)
        for (category) in (score_categories()):
            if (category != "ALL" and (gene in shared) != (category == "BOTH_CORE")):
                continue
            entry = bins.setdefault(
                (category, pattern), {"gene_count": 0, "interaction_contribution": Fraction(0)}
            )
            entry["gene_count"] += 1
            entry["interaction_contribution"] += item["allocation_interaction"]
    for (category) in (score_categories()):
        terms = decomposition([amounts(weights, shared, category) for (weights) in (decisions)])
        require(
            terms["I_total"] == terms["I_score_reference"] + terms["I_decision_reference"],
            "Signed decomposition failed",
        )
        require(
            terms["I_total"]
            == sum(
                (
                    entry["interaction_contribution"]
                    for ((cat, pattern), entry) in (bins.items())
                    if (cat == category)
                ),
                Fraction(0),
            ),
            "Pattern decomposition failed",
        )
        model_rows.append({**identity, "h": h, "category": category, **terms})
    for ((category, pattern), value) in (sorted(bins.items())):
        pattern_rows.append(
            {**identity, "h": h, "category": category, "pattern": pattern, **value}
        )
    return model_rows, summary, query_rows, cutoff_rows, pattern_rows, allocation_rows

def donor_accounting(model_rows, roster):
    counts = Counter(row["donor_id"] for (row) in (roster))
    lookup = {(row["sample_id"], row["h"], row["category"]): row for (row) in (model_rows)}
    donor_rows, population, distributions, crosses = [], [], [], []
    for (h) in ((0, 1)):
        for (category) in (score_categories()):
            donors = []
            for (donor) in (sorted(counts)):
                samples = [row["sample_id"] for (row) in (roster) if (row["donor_id"] == donor)]
                values = {
                    field: sum(
                        (lookup[(sample, h, category)][field] for (sample) in (samples)),
                        Fraction(0),
                    )
                    / len(samples)
                    for (field) in (fields())
                }
                row = {
                    "donor_id": donor,
                    "models": len(samples),
                    "h": h,
                    "category": category,
                    **values,
                }
                require(
                    row["I_total"] == row["I_score_reference"] + row["I_decision_reference"],
                    "Donor decomposition failed",
                )
                donors.append(row)
                donor_rows.append(row)
            population.append(
                {
                    "h": h,
                    "category": category,
                    **{
                        field: sum((row[field] for (row) in (donors)), Fraction(0)) / len(donors)
                        for (field) in (fields())
                    },
                }
            )
            for (field) in (fields()[5:]):
                distributions.append(
                    {
                        "h": h,
                        "category": category,
                        "quantity": field,
                        **statistics([row[field] for (row) in (donors)]),
                    }
                )
            for (p_sign, d_sign) in (product(("NEGATIVE", "ZERO", "POSITIVE"), repeat = 2)):
                selected = sum(
                    direction(row["P_at_D0"]) == p_sign and direction(row["D_at_P0"]) == d_sign
                    for (row) in (donors)
                )
                crosses.append(
                    {
                        "h": h,
                        "category": category,
                        "P_at_D0_sign": p_sign,
                        "D_at_P0_sign": d_sign,
                        "donors": selected,
                    }
                )
    return donor_rows, population, distributions, crosses

def direction(value):
    return "POSITIVE" if (value > 0) else "NEGATIVE" if (value < 0) else "ZERO"

def population_patterns(rows, roster):
    counts = Counter(row["donor_id"] for (row) in (roster))
    result = []
    for (h) in ((0, 1)):
        for (category) in (score_categories()):
            group = [row for (row) in (rows) if (row["h"] == h and row["category"] == category)]
            for (pattern) in (sorted({row["pattern"] for (row) in (group)})):
                selected = [row for (row) in (group) if (row["pattern"] == pattern)]
                result.append(
                    {
                        "h": h,
                        "category": category,
                        "pattern": pattern,
                        "donor_weighted_query_gene_count": sum(
                            (
                                Fraction(row["gene_count"], 83 * counts[row["donor_id"]])
                                for (row) in (selected)
                            ),
                            Fraction(0),
                        ),
                        "interaction_contribution": sum(
                            (
                                row["interaction_contribution"] / (83 * counts[row["donor_id"]])
                                for (row) in (selected)
                            ),
                            Fraction(0),
                        ),
                        "models_present": sum(row["gene_count"] > 0 for (row) in (selected)),
                    }
                )
    return result

def score_summaries(rows, models, model_rows):
    counts = Counter(row["donor_id"] for (row) in (models))
    lookup = {(row["sample_id"], row["h"], row["category"]): row for (row) in (model_rows)}
    output = []
    for (h) in ((0, 1)):
        selected = [row for (row) in (rows) if (row["h"] == h)]
        weighted = lambda field: sum(
            (row[field] / (83 * counts[row["donor_id"]]) for (row) in (selected)), Fraction(0)
        )
        ratios = [
            row["reference_cutoff_gap"] / (2 * row["reference_error_bound"])
            for (row) in (selected)
            if (row["reference_cutoff_gap"] is not None and row["reference_error_bound"] > 0)
        ]
        output.append(
            {
                "h": h,
                "models": len(selected),
                "top10_changed_models": sum(
                    row["actual_reference_changed"] for (row) in (selected)
                ),
                "changed_model_fraction": Fraction(
                    sum(row["actual_reference_changed"] for (row) in (selected)), 85
                ),
                "reference_certified_models": sum(
                    row["reference_cutoff_certified"] for (row) in (selected)
                ),
                "reference_uncertified_models": sum(
                    not row["reference_cutoff_certified"] for (row) in (selected)
                ),
                "donor_weighted_half_L1": weighted("actual_reference_half_L1"),
                "query_score_interaction_nonzero_models": sum(
                    row["query_interaction_certified_nonzero_n"] > 0 for (row) in (selected)
                ),
                "negative_reference_models": sum(
                    row["reference_negative_full_n"] > 0 for (row) in (selected)
                ),
                "focal_R_changed_models": sum(
                    lookup[(row["sample_id"], h, "BOTH_CORE")]["I_score_reference"] != 0
                    for (row) in (selected)
                ),
                "full_score_L1_mean": weighted("full_score_L1"),
                "query_score_L1_mean": weighted("query_score_L1"),
                "full_score_Linf_maximum": max(row["full_score_Linf"] for (row) in (selected)),
                "reference_gap_over_twice_error_minimum": min(ratios, default = None),
            }
        )
    return output

def analyse(directory, category_document, output):
    genes, samples, roster, query, accepted, errors, scores, shared = intake(
        directory, category_document
    )
    model_rows, score_rows, query_rows, cutoff_rows, pattern_rows, allocations = (
        [],
        [],
        [],
        [],
        [],
        [],
    )
    for (h) in ((0, 1)):
        for (index, identity) in (enumerate(roster)):
            values = model_reference(
                index, identity, genes, query, accepted, errors, scores, shared, h
            )
            model_rows.extend(values[0])
            score_rows.append(values[1])
            query_rows.extend(values[2])
            cutoff_rows.extend(values[3])
            pattern_rows.extend(values[4])
            allocations.extend(values[5])
    donors, population, distributions, crosses = donor_accounting(model_rows, roster)
    patterns = population_patterns(pattern_rows, roster)
    signed, vectors = summarise(donors)
    tables = {
        "model_allocation": model_rows,
        "model_score_diagnostics": score_rows,
        "query_gene_diagnostics": query_rows,
        "cutoff_certificates": cutoff_rows,
        "model_patterns": pattern_rows,
        "donor_allocation": donors,
        "population_allocation": population,
        "donor_distributions": distributions,
        "donor_simple_effect_cross": crosses,
        "population_patterns": patterns,
        "signed_mass_summary": signed,
        "donor_vectors": vectors,
    }
    for (name, rows) in (tables.items()):
        write_csv(output / (name + ".csv"), rows)
    write_json(output / "allocation_reference.json", {"states": allocations})
    result = {
        "models": 85,
        "donors": 83,
        "backgrounds": 2,
        "population": population,
        "donor_distributions": distributions,
        "score_summary": score_summaries(score_rows, roster, model_rows),
        "pattern_contributions": patterns,
        "signed_mass_summary": signed,
        "interpretation": "Exact native-score reference accounting, not unique causal percentages or a new graph",
    }
    write_json(output / "results.json", result)
    return result
