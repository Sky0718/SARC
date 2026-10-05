import struct
from collections import Counter, defaultdict
from fractions import Fraction
from .hypothesis_codec import require, prizes, mask_operator

RELEASES = ("native_11_0", "native_11_5", "native_12_0")
PAIRS = (
    (RELEASES[0], RELEASES[1]),
    (RELEASES[1], RELEASES[2]),
    (RELEASES[0], RELEASES[2]),
)
BUDGETS = (1, 5, 10, 20)
QUERIES = ("RAW_LINK", "AGGREGATE_INCLUDED_LINK")
STRATA = ("positive_retained", "equal_positive_weight", "whole_W_equal")
THRESHOLD = struct.unpack("<d", bytes.fromhex("9a9999999999e93f"))[0]

def token_equal(a, b):
    return a == b

def candidate_union(links, dictionary, unresolved):
    unknown = sorted(
        (p for (p) in (links) if (p not in dictionary or dictionary[p] is None))
    )
    members = (
        set().union(
            *(
                dictionary[p]
                for (p) in (links)
                if (p in dictionary and dictionary[p] is not None)
            )
        )
        if (links)
        else set()
    )
    return {
        "status": "COMPLETE_OBSERVED_UNIVERSE"
        if (not unresolved and (not unknown))
        else "PARTIAL_KNOWN_UNIVERSE",
        "members": members,
        "unknown_link_rows": unresolved,
        "unavailable_dictionary_rows": unknown,
    }

def universe_difference(old, new):
    complete = old["status"] == new["status"] == "COMPLETE_OBSERVED_UNIVERSE"
    return {
        "status": "COMPLETE" if (complete) else "PARTIAL_KNOWN_ONLY",
        "equal": old["members"] == new["members"] if (complete) else None,
        "old_known_members": old["members"],
        "new_known_members": new["members"],
        "known_member_entry": new["members"] - old["members"],
        "known_member_exit": old["members"] - new["members"],
        "old_state": old["status"],
        "new_state": new["status"],
        "old_unknown_link_rows": old["unknown_link_rows"],
        "new_unknown_link_rows": new["unknown_link_rows"],
        "old_unavailable_dictionary_rows": old["unavailable_dictionary_rows"],
        "new_unavailable_dictionary_rows": new["unavailable_dictionary_rows"],
    }

def gene_contrast(old, new, gene, dictionary):
    require(
        gene in old["columns"] and gene in new["columns"],
        "Retained nominated gene lacks native influence column",
    )
    old_index = {p: i for ((i, p)) in (enumerate(old["rows"]))}
    new_index = {p: i for ((i, p)) in (enumerate(new["rows"]))}
    common = sorted(old_index.keys() & new_index.keys())
    old_column, new_column = (old["columns"][gene], new["columns"][gene])
    same_support = set(old["rows"]) == set(new["rows"])
    equal_common = all(
        (
            token_equal(
                old_column["tokens"][old_index[p]], new_column["tokens"][new_index[p]]
            )
            for (p) in (common)
        )
    )
    mask_changes = [p for (p) in (common) if ((p in old["mask"]) != (p in new["mask"]))]
    cell_term = mask_term = support_term = Fraction()
    maximum = None
    query = {
        q: {
            "gain": [],
            "loss": [],
            "both_link": 0,
            "both_absent": 0,
            "unresolved": [],
            "crossings": [],
        }
        for (q) in (QUERIES)
    }
    prize_equal = []
    for (p) in (common):
        i, j = (old_index[p], new_index[p])
        v0, v1 = (old_column["values"][i], new_column["values"][j])
        x0, x1 = (Fraction(old_column["X"][i]), Fraction(new_column["X"][j]))
        m0, m1 = (int(p in old["mask"]), int(p in new["mask"]))
        cell_term += Fraction(m0 + m1, 2) * (x1 - x0)
        mask_term += (x0 + x1) * Fraction(m1 - m0, 2)
        finite = old_column["states"][i] == new_column["states"][j] == "F"
        if (finite):
            delta = Fraction(v1) - Fraction(v0)
            maximum = max(maximum or Fraction(), abs(delta))
        for (q) in (QUERIES):
            target = query[q]
            if (not finite):
                target["unresolved"].append(p)
                continue
            a = v0 > THRESHOLD and (q == "RAW_LINK" or m0)
            b = v1 > THRESHOLD and (q == "RAW_LINK" or m1)
            if (a == b):
                target["both_link" if (a) else "both_absent"] += 1
            else:
                target["gain" if (b) else "loss"].append(p)
                target["crossings"].append(
                    {
                        "pathway": p,
                        "old": old_column["tokens"][i],
                        "new": new_column["tokens"][j],
                        "old_minus_threshold": Fraction(v0) - Fraction(THRESHOLD),
                        "new_minus_threshold": Fraction(v1) - Fraction(THRESHOLD),
                        "old_mask": bool(m0),
                        "new_mask": bool(m1),
                    }
                )
        a, b = (prizes(old, p, dictionary.get(p)), prizes(new, p, dictionary.get(p)))
        available = (
            a["status"] == b["status"] == "RECOVERED_STORED_VECTOR_NO_ENRICHMENT"
        )
        prize_equal.append(
            {
                "pathway": p,
                "status": "AVAILABLE" if (available) else "UNAVAILABLE",
                "ordered_vector_equal": a["ordered"] == b["ordered"]
                if (available)
                else None,
                "name_keyed_vector_equal": dict(a["ordered"]) == dict(b["ordered"])
                if (available)
                else None,
                "normalising_exact_sum_equal": a["exact_sum"] == b["exact_sum"]
                if (available)
                else None,
                "old": {k: v for ((k, v)) in (a.items()) if (k != "ordered")},
                "new": {k: v for ((k, v)) in (b.items()) if (k != "ordered")},
                "old_prize_count": len(a["ordered"]) if (available) else None,
                "new_prize_count": len(b["ordered"]) if (available) else None,
            }
        )
    for (p) in (set(new["rows"]) - set(old["rows"])):
        if (p in new["mask"]):
            support_term += Fraction(new_column["X"][new_index[p]])
    for (p) in (set(old["rows"]) - set(new["rows"])):
        if (p in old["mask"]):
            support_term -= Fraction(old_column["X"][old_index[p]])
    exact_delta = new["totals"][gene][0] - old["totals"][gene][0]
    require(
        cell_term + mask_term + support_term == exact_delta,
        "Exact symmetric decomposition does not close",
    )
    monotone_eligible = (
        same_support
        and bool(common)
        and all(
            (
                old_column["states"][old_index[p]]
                == new_column["states"][new_index[p]]
                == "F"
                for (p) in (common)
            )
        )
    )
    nondecrease = (
        all(
            (
                Fraction(new_column["X"][new_index[p]])
                >= Fraction(old_column["X"][old_index[p]])
                for (p) in (common)
            )
        )
        if (monotone_eligible)
        else None
    )
    monotone_decline = nondecrease and exact_delta < 0 if (monotone_eligible) else None
    if (monotone_decline):
        require(bool(mask_changes), "Monotone-column decline must involve changed mask")
    for (q) in (QUERIES):
        universes = []
        for (cp, column) in (((old, old_column), (new, new_column))):
            links = {
                p
                for ((i, p)) in (enumerate(cp["rows"]))
                if (
                    column["states"][i] == "F"
                    and column["values"][i] > THRESHOLD
                    and (q == "RAW_LINK" or p in cp["mask"])
                )
            }
            unresolved = [
                p for ((i, p)) in (enumerate(cp["rows"])) if (column["states"][i] != "F")
            ]
            universes.append(candidate_union(links, dictionary, unresolved))
        query[q]["observed_change"] = bool(query[q]["gain"] or query[q]["loss"])
        query[q]["common_rows_fully_comparable"] = bool(common) and (
            not query[q]["unresolved"]
        )
        query[q]["universe"] = universe_difference(*universes)
    return {
        "gene": gene,
        "common_pathway_count": len(common),
        "old_pathway_exit": sorted(old_index.keys() - new_index.keys()),
        "new_pathway_entry": sorted(new_index.keys() - old_index.keys()),
        "native_mutation_axis_equal": old["genes"] == new["genes"],
        "old_mutation_axis": old["genes"],
        "new_mutation_axis": new["genes"],
        "max_abs_common_finite_delta": maximum,
        "same_native_pathway_support": same_support,
        "common_column_bit_state_equal": equal_common if (common) else None,
        "complete_column_bit_state_equal": same_support and equal_common,
        "changed_effective_mask_rows": mask_changes,
        "strict_nonlocal_mask_witness": same_support
        and equal_common
        and bool(mask_changes),
        "weak_common_column_mask_observation": bool(common)
        and equal_common
        and bool(mask_changes),
        "monotone_classification_eligible": monotone_eligible,
        "MONOTONE_OWN_COLUMN_AGGREGATE_DECLINE": monotone_decline,
        "old_native_branch": old["branch"],
        "new_native_branch": new["branch"],
        "own_cell_term": cell_term,
        "mask_term": mask_term,
        "support_term": support_term,
        "exact_aggregate_delta": exact_delta,
        "prize_comparisons": prize_equal,
        "queries": query,
    }

def common_rectangle(old, new):
    rows = sorted(set(old["rows"]) & set(new["rows"]))
    genes = sorted(set(old["genes"]) & set(new["genes"]))
    results = []
    for (cp) in ((old, new)):
        idx = {p: i for ((i, p)) in (enumerate(cp["rows"]))}
        columns = [[cp["columns"][g]["X"][idx[p]] for (p) in (rows)] for (g) in (genes)]
        try:
            mask, branch, totals = mask_operator(rows, genes, columns, False)
            results.append(
                {
                    "status": "AVAILABLE",
                    "branch": branch,
                    "analytical_rectangle_mask": mask,
                    "native_restricted_mask": cp["mask"] & set(rows),
                    "changed_from_native": mask ^ cp["mask"] & set(rows),
                }
            )
        except ValueError as error:
            results.append(
                {
                    "status": "UNAVAILABLE_LITERAL_SOURCE_BRANCH",
                    "reason": str(error),
                    "native_restricted_mask": cp["mask"] & set(rows),
                }
            )
    return {
        "status": "RECTANGULAR_ALGEBRA_NOT_NATIVE_RUN",
        "pathways": rows,
        "genes": genes,
        "old": results[0],
        "new": results[1],
    }

def model_weights(models):
    counts = Counter((r["study_donor_id"] for (r) in (models)))
    result = {
        r["sample_ID"]: Fraction(1, len(counts) * counts[r["study_donor_id"]])
        for (r) in (models)
    }
    require(sum(result.values()) == 1, "Donor-derived weights do not total1")
    return result

def pair_record(a, b, k, metadata, get_contrast):
    w0, w1 = (a["weights"][k], b["weights"][k])
    base = dict(
        metadata,
        budget = k,
        offered_budget = k,
        status = "AVAILABLE",
        old_returned_mass = None if (w0 is None) else sum(w0.values(), Fraction()),
        new_returned_mass = None if (w1 is None) else sum(w1.values(), Fraction()),
        old_menu_count = None if (w0 is None) else len(w0),
        new_menu_count = None if (w1 is None) else len(w1),
    )
    if (w0 is None or w1 is None):
        return dict(
            base,
            status = "UNAVAILABLE_RANK",
            unavailable_capacity_bound = Fraction(k),
            strata = {},
        )
    whole = w0 == w1
    rho = {
        g: min(w0.get(g, Fraction()), w1.get(g, Fraction()))
        for (g) in (w0.keys() | w1.keys())
    }
    rho = {g: x for ((g, x)) in (rho.items()) if (x > 0)}
    base.update(
        {
            "whole_W_equal": whole,
            "whole_W_equal_semantics": "SAME_DOMAIN_W_EQUAL_SUFFICIENT_NOT_NECESSARY_ZERO_EXTENDED_ACTION_EQUALITY",
            "rho": rho,
            "retained_mass": sum(rho.values(), Fraction()),
            "old_positive_exit": sorted(
                (g for (g) in (w0) if (w0[g] > 0 and w1.get(g, 0) == 0))
            ),
            "new_positive_entry": sorted(
                (g for (g) in (w1) if (w1[g] > 0 and w0.get(g, 0) == 0))
            ),
            "strata": {},
            "retained_gene_contrast_ids": [],
        }
    )
    available = a["status"] == b["status"] == "AVAILABLE"
    contrasts = {g: get_contrast(g) for (g) in (rho)} if (available) else {}
    base["retained_gene_contrast_ids"] = [
        c["contrast_id"] for (c) in (contrasts.values())
    ]
    if (not available):
        base["status"] = "UNAVAILABLE_HYPOTHESIS_DERIVATION"
    for (stratum) in (STRATA):
        selected = {
            g: v
            for ((g, v)) in (rho.items())
            if (
                stratum == "positive_retained"
                or (stratum == "equal_positive_weight" and w0[g] == w1[g])
                or (stratum == "whole_W_equal" and whole)
            )
        }
        total = sum(selected.values(), Fraction())
        base["strata"][stratum] = {}
        for (q) in (QUERIES):
            item = {
                "retained_mass": total,
                "observed_changed_mass": Fraction(),
                "comparable_retained_mass": Fraction(),
                "unavailable_retained_mass": Fraction(),
                "partially_observed_mass": Fraction(),
                "no_common_pathway_mass": Fraction(),
                "shared_support_universe": None,
            }
            for (g, mass) in (selected.items()):
                if (not available):
                    item["unavailable_retained_mass"] += mass
                    continue
                contrast = contrasts[g]
                value = contrast["queries"][q]
                if (value["observed_change"]):
                    item["observed_changed_mass"] += mass
                if (value["common_rows_fully_comparable"]):
                    item["comparable_retained_mass"] += mass
                else:
                    item["unavailable_retained_mass"] += mass
                    if (contrast["common_pathway_count"] == 0):
                        item["no_common_pathway_mass"] += mass
                    else:
                        item["partially_observed_mass"] += mass
            if (available):
                universes = [
                    contrasts[g]["queries"][q]["universe"] for (g) in (selected)
                ]
                old = (
                    set().union(*(u["old_known_members"] for (u) in (universes)))
                    if (universes)
                    else set()
                )
                new = (
                    set().union(*(u["new_known_members"] for (u) in (universes)))
                    if (universes)
                    else set()
                )
                complete = all((u["status"] == "COMPLETE" for (u) in (universes)))
                item["shared_support_universe"] = {
                    "status": "COMPLETE" if (complete) else "PARTIAL_KNOWN_ONLY",
                    "old_members": old,
                    "new_members": new,
                    "entry": new - old,
                    "exit": old - new,
                    "equal": old == new if (complete) else None,
                    "semantics": "TIE_INCLUSIVE_SUPPORT_UNIVERSE_NOT_REALIZED_OR_EXPECTED_PANEL",
                }
            base["strata"][stratum][q] = item
    return base

def summarize(records, weights):
    groups = {}
    for (pair) in (PAIRS):
        for (k) in (BUDGETS):
            rows = [
                r
                for (r) in (records)
                if (tuple(r["release_pair"]) == pair and r["budget"] == k)
            ]
            require(len(rows) == len(weights) * 5, "Incomplete pair/budget rows")
            result = {
                "release_pair": pair,
                "budget": k,
                "models": len(weights),
                "states": len(rows),
                "rank_unavailable_population_mass": Fraction(),
                "rank_unavailable_capacity_bound": Fraction(),
                "no_retained_population_mass": Fraction(),
                "strata": {
                    s: {q: defaultdict(Fraction) for (q) in (QUERIES)}
                    for (s) in (STRATA)
                },
            }
            for (row) in (rows):
                weight = weights[row["sample"]] / 5
                if (row["status"] == "UNAVAILABLE_RANK"):
                    result["rank_unavailable_population_mass"] += weight
                    result["rank_unavailable_capacity_bound"] += weight * k
                    continue
                if (row["retained_mass"] == 0):
                    result["no_retained_population_mass"] += weight
                for (s) in (STRATA):
                    for (q) in (QUERIES):
                        for (key, value) in (row["strata"][s][q].items()):
                            if (isinstance(value, Fraction)):
                                result["strata"][s][q][key] += weight * value
            for (s) in (STRATA):
                for (q) in (QUERIES):
                    entry = result["strata"][s][q]
                    total = entry["retained_mass"]
                    entry["conditional_observed_changed_fraction"] = (
                        entry["observed_changed_mass"] / total if (total) else None
                    )
                    entry["secondary_per_offered_budget_scale"] = (
                        entry["observed_changed_mass"] / k
                    )
            groups["|".join((*pair, str(k)))] = result
    return groups
