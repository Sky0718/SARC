from collections import Counter, defaultdict
from fractions import Fraction
from .retained_hypotheses import PAIRS, BUDGETS, STRATA

QUERIES = ("RAW_LINK", "AGGREGATE_INCLUDED_LINK")

def decode(value):
    return (
        Fraction(value["numerator"], value["denominator"])
        if (set(value) == {"numerator", "denominator"})
        else value
    )

def tag(value):
    return (
        "NULL"
        if (value is None)
        else str(value).upper()
        if (isinstance(value, bool))
        else str(value)
    )

def joint(*values):
    return "|".join((tag(v) for (v) in (values)))

def empirical_quartiles(counts):
    size = sum(counts.values())
    result = {"Q25": None, "Q50": None, "Q75": None}
    if (size == 0):
        return result
    ranks = {"Q25": (size + 3) // 4, "Q50": (size + 1) // 2, "Q75": (3 * size + 3) // 4}
    cumulative = 0
    for (value, count) in (sorted(counts.items())):
        cumulative += count
        for (key, rank) in (ranks.items()):
            if (result[key] is None and cumulative >= rank):
                result[key] = value
    return result

class Census:
    def __init__(self):
        self.records = 0
        self.categories = defaultdict(Counter)
        self.numeric = {}
        self.frequencies = defaultdict(Counter)

    def add(self, key, value):
        if (isinstance(value, (int, Fraction)) and (not isinstance(value, bool))):
            if (key not in self.numeric):
                self.numeric[key] = {
                    "n": 0,
                    "negative": 0,
                    "zero": 0,
                    "positive": 0,
                    "min": value,
                    "max": value,
                }
            item = self.numeric[key]
            item["n"] += 1
            item[
                "negative" if (value < 0) else "positive" if (value > 0) else "zero"
            ] += 1
            item["min"], item["max"] = (
                min(item["min"], value),
                max(item["max"], value),
            )
            self.frequencies[key][value] += 1
        else:
            self.categories[key][tag(value)] += 1

    def feature(self, features):
        self.records += 1
        for (key, value) in (features.items()):
            self.add(key, value)

    def result(self):
        numeric = {
            key: dict(value, **empirical_quartiles(self.frequencies[key]))
            for ((key, value)) in (self.numeric.items())
        }
        return {
            "records": self.records,
            "categorical_counts": dict(self.categories),
            "exact_numeric_ranges_and_sign_counts": numeric,
            "quartile_definition": "Inverse unweighted empirical CDF: sorted observation at one-based ceil(p*n), p=1/4,1/2,3/4; no interpolation; null values excluded from n and separately counted",
        }

def profile(items, key):
    return (
        ",".join(sorted({tag(item[key]) for (item) in (items)})) if (items) else "EMPTY"
    )

def gene_features(value):
    scalar = (
        "common_pathway_count",
        "native_mutation_axis_equal",
        "max_abs_common_finite_delta",
        "same_native_pathway_support",
        "common_column_bit_state_equal",
        "complete_column_bit_state_equal",
        "strict_nonlocal_mask_witness",
        "weak_common_column_mask_observation",
        "monotone_classification_eligible",
        "MONOTONE_OWN_COLUMN_AGGREGATE_DECLINE",
        "old_native_branch",
        "new_native_branch",
        "own_cell_term",
        "mask_term",
        "support_term",
        "exact_aggregate_delta",
    )
    result = {key: value[key] for (key) in (scalar)}
    for (key) in ((
        "old_pathway_exit",
        "new_pathway_entry",
        "old_mutation_axis",
        "new_mutation_axis",
        "changed_effective_mask_rows",
    )):
        result[key + "_count"] = len(value[key])
    prizes = value["prize_comparisons"]
    result["prize_comparison_count"] = len(prizes)
    for (key) in ((
        "status",
        "ordered_vector_equal",
        "name_keyed_vector_equal",
        "normalising_exact_sum_equal",
    )):
        result["prize_" + key + "_observed_states"] = profile(prizes, key)
    for (query) in (QUERIES):
        q = value["queries"][query]
        for (key) in ((
            "observed_change",
            "common_rows_fully_comparable",
            "both_link",
            "both_absent",
        )):
            result[query + "." + key] = q[key]
        for (key) in (("gain", "loss", "unresolved", "crossings")):
            result[query + "." + key + "_count"] = len(q[key])
        universe = q["universe"]
        for (key) in (("status", "equal", "old_state", "new_state")):
            result[query + ".universe." + key] = universe[key]
        for (key) in ((
            "old_known_members",
            "new_known_members",
            "known_member_entry",
            "known_member_exit",
            "old_unknown_link_rows",
            "new_unknown_link_rows",
            "old_unavailable_dictionary_rows",
            "new_unavailable_dictionary_rows",
        )):
            result[query + ".universe." + key + "_count"] = len(universe[key])
        result[query + ".joint_common_link_universe"] = joint(
            q["observed_change"],
            q["common_rows_fully_comparable"],
            universe["status"],
            universe["equal"],
        )
    result["joint_support_axis_prize_raw_included_strict_monotone"] = joint(
        value["same_native_pathway_support"],
        value["native_mutation_axis_equal"],
        profile(prizes, "name_keyed_vector_equal"),
        value["queries"]["RAW_LINK"]["observed_change"],
        value["queries"]["AGGREGATE_INCLUDED_LINK"]["observed_change"],
        value["strict_nonlocal_mask_witness"],
        value["MONOTONE_OWN_COLUMN_AGGREGATE_DECLINE"],
    )
    return result

def gene_events(census, value):
    for (prize) in (value["prize_comparisons"]):
        for (key) in ((
            "status",
            "ordered_vector_equal",
            "name_keyed_vector_equal",
            "normalising_exact_sum_equal",
            "old_prize_count",
            "new_prize_count",
        )):
            census.add("prize_pathway_occurrence." + key, prize.get(key))
        census.add(
            "prize_pathway_occurrence.joint",
            joint(
                prize["status"],
                prize["ordered_vector_equal"],
                prize["name_keyed_vector_equal"],
                prize["normalising_exact_sum_equal"],
            ),
        )
    for (query) in (QUERIES):
        q = value["queries"][query]
        gains, losses = (set(q["gain"]), set(q["loss"]))
        for (crossing) in (q["crossings"]):
            direction = (
                "GAIN"
                if (crossing["pathway"] in gains)
                else "LOSS"
                if (crossing["pathway"] in losses)
                else "UNASSIGNED"
            )
            prefix = query + ".crossing." + direction
            census.add(
                prefix + ".mask_pair", joint(crossing["old_mask"], crossing["new_mask"])
            )
            for (key) in (("old_minus_threshold", "new_minus_threshold")):
                census.add(prefix + "." + key, crossing[key])

def pair_features(row, stratum, query):
    result = {"record_status": row["status"], "whole_W_equal": row.get("whole_W_equal")}
    item = row["strata"].get(stratum, {}).get(query)
    if (item is None):
        return dict(result, stratum_state = "NOT_AVAILABLE")
    result["stratum_state"] = (
        "ZERO_RETAINED" if (item["retained_mass"] == 0) else "POSITIVE_RETAINED"
    )
    universe = item["shared_support_universe"]
    result["shared_universe_status"] = (
        "NOT_AVAILABLE" if (universe is None) else universe["status"]
    )
    result["shared_universe_equal"] = None if (universe is None) else universe["equal"]
    result["joint_retention_shared_universe"] = joint(
        result["stratum_state"],
        result["shared_universe_status"],
        result["shared_universe_equal"],
    )
    for (key) in ((
        "retained_mass",
        "observed_changed_mass",
        "comparable_retained_mass",
        "unavailable_retained_mass",
        "partially_observed_mass",
        "no_common_pathway_mass",
    )):
        result[key] = item[key]
    for (key) in (("old_members", "new_members", "entry", "exit")):
        result["shared_universe_" + key + "_count"] = (
            None if (universe is None) else len(universe[key])
        )
    return result

def rectangle_features(row):
    result = {
        "status": row["status"],
        "common_pathway_count": len(row["pathways"]),
        "common_gene_count": len(row["genes"]),
    }
    for (side) in (("old", "new")):
        value = row[side]
        for (key) in (("status", "branch", "reason")):
            result[side + "." + key] = value.get(key)
        for (key) in ((
            "analytical_rectangle_mask",
            "native_restricted_mask",
            "changed_from_native",
        )):
            result[side + "." + key + "_count"] = (
                None if (key not in value) else len(value[key])
            )
    result["joint_side_status"] = joint(row["old"]["status"], row["new"]["status"])
    return result

def restore(value):
    if (isinstance(value, dict)):
        if (set(value) == {"numerator", "denominator"}):
            return decode(value)
        return {key: restore(item) for ((key, item)) in (value.items())}
    if (isinstance(value, list)):
        return [restore(item) for (item) in (value)]
    return value

def analyse(derivative):
    data = restore(derivative)
    gene_groups = {
        "|".join((*pair, str(budget))): Census()
        for (pair) in (PAIRS)
        for (budget) in ((*BUDGETS, "ALL_BUDGET_UNION"))
    }
    pair_groups = {
        "|".join((*pair, str(budget), stratum, query)): Census()
        for (pair) in (PAIRS)
        for (budget) in (BUDGETS)
        for (stratum) in (STRATA)
        for (query) in (QUERIES)
    }
    rectangle_groups = {"|".join(pair): Census() for (pair) in (PAIRS)}
    references, reference_occurrences, totals = defaultdict(set), Counter(), Counter()
    output = []
    for (row) in (data["pair_seed_budget"]):
        totals["pair_seed_budget_records"] += 1
        group = "|".join((*row["release_pair"], str(row["budget"])))
        if (group not in gene_groups):
            raise ValueError("Unexpected release pair or budget")
        for (identity) in (row.get("retained_gene_contrast_ids", [])):
            references[identity].add(group)
            reference_occurrences[group] += 1
        for (stratum) in (STRATA):
            for (query) in (QUERIES):
                pair_groups[group + "|" + stratum + "|" + query].feature(
                    pair_features(row, stratum, query)
                )
    for (row) in (data["gene_contrasts"]):
        totals["distinct_lineage_gene_records"] += 1
        identity = row["contrast_id"]
        parts = identity.split("|")
        groups = references.pop(identity)
        if (not groups):
            raise ValueError("Unreferenced gene contrast")
        features = gene_features(row)
        output.append(
            {
                "contrast_id": identity,
                "referenced_pair_budget_groups": sorted(groups),
                "features": features,
            }
        )
        groups.add("|".join((*parts[1:3], "ALL_BUDGET_UNION")))
        for (group) in (groups):
            gene_groups[group].feature(features)
            gene_events(gene_groups[group], row)
    if (references):
        raise ValueError("Referenced contrast absent from complete gene carrier")
    for (row) in (data["common_rectangles"]):
        totals["common_rectangle_records"] += 1
        pair = "|".join(row["identity"].split("|")[1:3])
        rectangle_groups[pair].feature(rectangle_features(row))
    return {
        "semantics": "Unweighted descriptive record census; checkpoint contrasts deduplicated within pair and budget; not biological replications, effect estimates or p-values",
        "input_record_counts": dict(totals),
        "gene_reference_occurrences_before_within_group_deduplication": dict(
            reference_occurrences
        ),
        "gene_contrast_groups": {
            key: value.result() for ((key, value)) in (gene_groups.items())
        },
        "pair_seed_budget_stratum_query_groups": {
            key: value.result() for ((key, value)) in (pair_groups.items())
        },
        "common_rectangle_groups": {
            key: value.result() for ((key, value)) in (rectangle_groups.items())
        },
        "gene_field_census": output,
    }
