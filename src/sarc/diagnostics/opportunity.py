from collections import Counter, defaultdict
from fractions import Fraction

SEEDS = [104729, 130363, 155921, 196613, 228017]

CELLS = ["G0_H0", "G0_H1", "G1_H0", "G1_H1"]

MASS_FIELDS = [
    "pro_mass",
    "dawn_mass",
    "intersection",
    "union",
    "menu_multiset_upper",
    "baseline_capacity_upper",
    "combined_upper",
    "uniform_expected_intersection",
    "uniform_expected_union",
    "common_query_expected_intersection",
    "common_query_expected_union",
]

def need(value, message):
    if (not value):
        raise ValueError(message)

def exact(value):
    need(
        isinstance(value, (str, int, Fraction)) and not isinstance(value, bool),
        "Exact fraction required",
    )
    return Fraction(value)

def weights(rows):
    need(
        rows is not None and len({row[0] for (row) in (rows)}) == len(rows),
        "Complete unique allocation required",
    )
    result = {gene: exact(value) for (gene, value) in (rows)}
    need(all(0 <= value <= 1 for (value) in (result.values())), "Allocation outside unit capacity")
    return result

def ratio(a, b):
    return str(a / b) if (b) else None

def action_record(p, d, pmenu, dmenu, query, pold, dold, action):
    pmenu, dmenu, query = set(pmenu), set(dmenu), set(query)
    need(pmenu <= query and dmenu <= query, "Emitted menu outside query")
    need(all(value >= 0 for (value) in ([*p.values(), *d.values()])), "Negative action mass")
    need(
        {gene for (gene, value) in (p.items()) if (value)} <= pmenu, "PRO action outside emitted menu"
    )
    need(
        {gene for (gene, value) in (d.items()) if (value)} <= dmenu,
        "Dawn action outside emitted menu",
    )
    pvalues = sorted((value for (value) in (p.values()) if (value)), reverse = True)
    dvalues = sorted((value for (value) in (d.values()) if (value)), reverse = True)
    ptotal, dtotal = sum(pvalues, Fraction()), sum(dvalues, Fraction())
    intersection = sum((min(p.get(gene, 0), d.get(gene, 0)) for (gene) in (query)), Fraction())
    union = ptotal + dtotal - intersection
    shared_count = len(pmenu & dmenu)
    relaxed = sum(
        (min(a, b) for (a, b) in (list(zip(pvalues, dvalues))[:shared_count])), Fraction()
    )
    pcaps = {
        gene: 1 - pold.get(gene, 0) if (action == "addition") else pold.get(gene, 0)
        for (gene) in (pmenu)
    }
    dcaps = {
        gene: 1 - dold.get(gene, 0) if (action == "addition") else dold.get(gene, 0)
        for (gene) in (dmenu)
    }
    need(
        all(value <= pcaps.get(gene, 0) for (gene, value) in (p.items())),
        "PRO baseline capacity violated",
    )
    need(
        all(value <= dcaps.get(gene, 0) for (gene, value) in (d.items())),
        "Dawn baseline capacity violated",
    )
    capacity = min(
        ptotal,
        dtotal,
        sum((min(pcaps[gene], dcaps[gene]) for (gene) in (pmenu & dmenu)), Fraction()),
    )
    combined = min(relaxed, capacity)
    need(intersection <= combined, "Observed overlap exceeds a relaxation")
    pair_min = sum((min(a, b) for (a) in (pvalues) for (b) in (dvalues)), Fraction())
    expected = (
        pair_min * Fraction(shared_count, len(pmenu) * len(dmenu))
        if (pmenu and dmenu)
        else Fraction()
    )
    qexpected = pair_min / len(query) if (query) else Fraction()
    need(
        expected <= relaxed and qexpected <= min(ptotal, dtotal),
        "Uniform expectation exceeds maximum",
    )
    values = dict(
        zip(
            MASS_FIELDS,
            [
                ptotal,
                dtotal,
                intersection,
                union,
                relaxed,
                capacity,
                combined,
                expected,
                ptotal + dtotal - expected,
                qexpected,
                ptotal + dtotal - qexpected,
            ],
        )
    )
    result = {key: str(value) for (key, value) in (values.items())}
    result.update(
        query_count = len(query),
        pro_menu_count = len(pmenu),
        dawn_menu_count = len(dmenu),
        common_menu_count = shared_count,
        pro_positive_capacity_count = sum(value > 0 for (value) in (pcaps.values())),
        dawn_positive_capacity_count = sum(value > 0 for (value) in (dcaps.values())),
        observed_jaccard = ratio(intersection, union),
        observed_fraction_of_combined_upper = ratio(intersection, combined),
        uniform_ratio_of_expected_masses = ratio(expected, values["uniform_expected_union"]),
        common_query_ratio_of_expected_masses = ratio(
            qexpected, values["common_query_expected_union"]
        ),
        gene_ledger = [
            {
                "gene": gene,
                "pro_action": str(p.get(gene, 0)),
                "dawn_action": str(d.get(gene, 0)),
                "pro_menu": gene in pmenu,
                "dawn_menu": gene in dmenu,
                "pro_baseline_capacity": str(pcaps.get(gene, 0)),
                "dawn_baseline_capacity": str(dcaps.get(gene, 0)),
            }
            for (gene) in (sorted(query))
        ],
    )
    return result

def aggregate(rows, row_weights):
    need(
        len(rows) == len(row_weights) and sum(row_weights, Fraction()) == 1,
        "Complete normalised aggregation required",
    )
    values = {
        key: sum(
            (exact(row[key]) * weight for (row, weight) in (zip(rows, row_weights))), Fraction()
        )
        for (key) in (MASS_FIELDS)
    }
    output = {key: str(value) for (key, value) in (values.items())}
    for (key) in ([
        "query_count",
        "pro_menu_count",
        "dawn_menu_count",
        "common_menu_count",
        "pro_positive_capacity_count",
        "dawn_positive_capacity_count",
    ]):
        output["weighted_" + key] = str(
            sum(
                (Fraction(row[key]) * weight for (row, weight) in (zip(rows, row_weights))),
                Fraction(),
            )
        )
    output.update(
        observed_jaccard = ratio(values["intersection"], values["union"]),
        observed_fraction_of_combined_upper = ratio(
            values["intersection"], values["combined_upper"]
        ),
        uniform_ratio_of_expected_masses = ratio(
            values["uniform_expected_intersection"], values["uniform_expected_union"]
        ),
        common_query_ratio_of_expected_masses = ratio(
            values["common_query_expected_intersection"], values["common_query_expected_union"]
        ),
    )
    return output

def analyse(prepared):
    donor_counts = Counter(row["model"]["donor_id"] for (row) in (prepared))
    output = {
        "schema": "joint_menu_opportunity_diagnostics_v1",
        "counts": {
            "models": len(prepared),
            "donors": len(donor_counts),
            "cells_per_method": len(prepared) * 4,
        },
        "methods": ["PRODIGY", "DawnRank"],
        "contrast": "joint",
        "models": [],
        "donors": [],
        "population": {},
        "interpretation": "Exploratory opportunity accounting. Menu-conditioned relaxations and mechanical uniform references are not pipeline-reachable effects or biological randomisation tests.",
    }
    for (item) in (prepared):
        model = item["model"]
        query = set(model["eligible_genes"])
        record = {key: model[key] for (key) in (["model_id", "sample_id", "donor_id"])}
        record["model_population_weight"] = str(
            Fraction(1, len(donor_counts) * donor_counts[model["donor_id"]])
        )
        record["query_count"] = len(query)
        record["cells"] = {
            cell: {
                "pro_emitted_menu_count": len(item["pro"][cell]["menu"]),
                "pro_seed_menu_counts": item["pro"][cell]["seed_menu_counts"],
                "pro_seed_menus_identical": item["pro"][cell]["seed_menus_identical"],
                "pro_statuses": item["pro"][cell]["statuses"],
                "dawn_emitted_menu_count": len(item["dawn"][cell]["menu"]),
                "pro_returned_mass": str(
                    sum(item["pro"][cell]["weights"].values(), Fraction())
                ),
                "dawn_returned_mass": str(
                    sum(item["dawn"][cell]["weights"].values(), Fraction())
                ),
            }
            for (cell) in (CELLS)
        }
        record["actions"] = {}
        for (action, sign, menu_cell) in ([("addition", 1, "G1_H1"), ("removal", -1, "G0_H0")]):
            vectors = {
                name: {
                    gene: max(
                        sign
                        * (
                            item[name]["G1_H1"]["weights"].get(gene, 0)
                            - item[name]["G0_H0"]["weights"].get(gene, 0)
                        ),
                        0,
                    )
                    for (gene) in (query)
                }
                for (name) in (["pro", "dawn"])
            }
            record["actions"][action] = action_record(
                vectors["pro"],
                vectors["dawn"],
                item["pro"][menu_cell]["menu"],
                item["dawn"][menu_cell]["menu"],
                query,
                item["pro"]["G0_H0"]["weights"],
                item["dawn"]["G0_H0"]["weights"],
                action,
            )
        record["dawn_cutoff_diagnostics"] = item["cutoff"]
        output["models"].append(record)
    for (donor) in (sorted(donor_counts)):
        rows = [model for (model) in (output["models"]) if (model["donor_id"] == donor)]
        output["donors"].append(
            {
                "donor_id": donor,
                "model_ids": [row["model_id"] for (row) in (rows)],
                "actions": {
                    action: aggregate(
                        [row["actions"][action] for (row) in (rows)],
                        [Fraction(1, len(rows))] * len(rows),
                    )
                    for (action) in (["addition", "removal"])
                },
            }
        )
    row_weights = [exact(row["model_population_weight"]) for (row) in (output["models"])]
    output["population"] = {
        action: aggregate([row["actions"][action] for (row) in (output["models"])], row_weights)
        for (action) in (["addition", "removal"])
    }
    output["cutoff_summary"] = {}
    for (cell) in (CELLS):
        rows = [
            row
            for (model) in (output["models"])
            for (row) in (model["dawn_cutoff_diagnostics"])
            if (row["cell"] == cell)
        ]
        gaps = sorted(exact(row["cutoff_gap"]) for (row) in (rows) if (row["cutoff_gap"] is not None))
        need(len(rows) == len(prepared), "One cutoff record per model and cell required")
        output["cutoff_summary"][cell] = {
            "model_count": len(rows),
            "gap_evaluable_models": len(gaps),
            "fewer_than_eleven_candidates": len(rows) - len(gaps),
            "boundary_ties": sum(row["boundary_tie"] is True for (row) in (rows)),
            "positive_gaps": sum(value > 0 for (value) in (gaps)),
            "minimum_recorded_gap": str(gaps[0]) if (gaps) else None,
            "median_recorded_gap": str((gaps[(len(gaps) - 1) // 2] + gaps[len(gaps) // 2]) / 2)
            if (gaps)
            else None,
            "maximum_recorded_gap": str(gaps[-1]) if (gaps) else None,
            "iterations_min": min(row["iterations"] for (row) in (rows)),
            "iterations_max": max(row["iterations"] for (row) in (rows)),
            "fixed_point_or_cutoff_certificates": sum(
                row.get("cutoff_certified") is True for (row) in (rows)
            ),
        }
    return output
