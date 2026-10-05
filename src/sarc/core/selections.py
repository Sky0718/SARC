import csv
import gzip
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

from ..transfer.policies import consensus
from ..transfer.scoring import leading_membership
from .io import SEEDS, load_json, resolve_path

METHODS = ("DawnRank", "PRODIGY", "PersonaDrive")
POOLS = ("all34", "distal16", "proximal18")
RELEASES = ("11_0", "11_5", "12_0")
POLICIES = ("FIXED_V0", "FIXED_V1", "FIXED_V2", "RC", "STABLE_ONLY_3", "V0_DEGREE")
BUDGETS = (1, 5, 10, 20)
EIGHT = tuple(
    "ACH-" + value
    for (value) in (
        (
            "003332",
            "003334",
            "003337",
            "003346",
            "003387",
            "003400",
            "003401",
            "003402",
        )
    )
)
THREE = ("HCM-SANG-0265-C18", "HCM-SANG-0270-C20", "HCM-SANG-0286-C20")
SPLITS = ("development", "held_out", "held_out")

def require(value, message):
    if (not value):
        raise ValueError(message)

def read_rows(path):
    path = Path(path)
    opener = gzip.open if (path.suffix == ".gz") else open
    suffixes = path.suffixes
    with opener(path, "rt", encoding = "utf-8-sig") as stream:
        result = (
            [json.loads(line) for (line) in (stream) if (line.strip())]
            if (".jsonl" in suffixes)
            else json.load(stream)
        )
    require(isinstance(result, list), "Expected an explicit list of rank records")
    return result

def source_ids(config, fields):
    identifiers = config["source_ids"]
    require(set(identifiers) == set(fields), "All scientific source IDs are required")
    require(
        all(
            isinstance(value, str) and value.strip() == value and value
            for (value) in (identifiers.values())
        ),
        "Source IDs must be non-empty literal identifiers",
    )
    return dict(identifiers)

def mutation_matrix(path):
    opener = gzip.open if (str(path).endswith(".gz")) else open
    with opener(path, "rt", encoding = "utf-8-sig", newline = "") as stream:
        reader = csv.reader(stream)
        require(
            next(reader)[1:] == list(EIGHT),
            "The exact eight-model mutation order is required",
        )
        result = {}
        for (row) in (reader):
            require(
                len(row) == 9
                and row[0]
                and row[0].strip() == row[0]
                and row[0] not in result,
                "Invalid mutation gene axis",
            )
            require(
                all(value in ("0", "1", "0.0", "1.0") for (value) in (row[1:])),
                "Mutation calls must be binary",
            )
            result[row[0]] = tuple(int(float(value)) for (value) in (row[1:]))
    return result

def build_eligibility(counts, persona, source_ledger):
    shared = set(counts) & set(persona)
    require(
        all(counts[gene] == persona[gene] for (gene) in (shared)),
        "Shared mutation calls differ",
    )
    merged = dict(counts, **persona)
    samples = {
        sample: {
            "native_eligible_literal_gene_labels": sorted(
                gene for ((gene, calls)) in (merged.items()) if (calls[index])
            ),
            "common_eligible_literal_gene_labels": sorted(
                gene for (gene) in (shared) if (counts[gene][index])
            ),
        }
        for ((index, sample)) in (enumerate(EIGHT))
    }
    require(
        set(source_ledger["complete_qualifying_mapped_genes"]) == set(EIGHT),
        "Complete source mutation populations are required",
    )
    exclusions = []
    for (index, sample) in (enumerate(EIGHT)):
        positives = set(samples[sample]["native_eligible_literal_gene_labels"])
        source_positive = set(source_ledger["complete_qualifying_mapped_genes"][sample])
        require(
            positives <= source_positive,
            "A native positive is absent from the source mutation ledger",
        )
        for (gene) in (sorted(source_positive)):
            in_counts, in_persona = gene in counts, gene in persona
            require(
                not in_counts or counts[gene][index] == 1,
                "Count-axis mutation differs from source call",
            )
            require(
                not in_persona or persona[gene][index] == 1,
                "TPM-axis mutation differs from source call",
            )
            exclusions.append(
                {
                    "sample_id": sample,
                    "gene": gene,
                    "in_count_input_axis": in_counts,
                    "in_tpm_input_axis": in_persona,
                    "native_eligible": gene in positives,
                    "disposition": "ADMITTED_NATIVE_INPUT"
                    if (gene in positives)
                    else "QUALIFYING_SOURCE_CALL_OUTSIDE_BOTH_NATIVE_INPUT_AXES",
                }
            )
    return {
        "schema": "NEXTGEN8_INPUT_ONLY_NONRELEASE_ELIGIBILITY_V1",
        "samples": list(EIGHT),
        "pools": {pool: samples for (pool) in (POOLS)},
        "native_rule": "Union of count-source and TPM-source mutation-positive axes before method/network restriction; degree zero retained",
        "common_rule": "Intersection of count-source and TPM-source axes before graph restriction, used for inherited RC",
        "count_axis_genes": len(counts),
        "persona_axis_genes": len(persona),
        "common_axis_genes": len(shared),
        "complete_source_positive_input_ledger": exclusions,
        "functional_values_opened": False,
        "ranking_inputs_used": False,
        "graph_support_used": False,
    }

def eight_model_eligibility(config, base):
    identifiers = source_ids(
        config, ("count_mutation", "persona_mutation", "mutation_gene_ledger")
    )
    counts = mutation_matrix(resolve_path(config["count_mutation"], base))
    persona = mutation_matrix(resolve_path(config["persona_mutation"], base))
    result = build_eligibility(
        counts, persona, load_json(resolve_path(config["mutation_gene_ledger"], base))
    )
    result["source_ids"] = identifiers
    return result

def degree_from_rows(rows):
    degree, seen = Counter(), set()
    for (row) in (rows):
        first, second = row["gene1"], row["gene2"]
        require(first and second and first != second, "Invalid baseline graph edge")
        edge = tuple(sorted((first, second)))
        require(edge not in seen, "Repeated undirected baseline edge")
        seen.add(edge)
        degree[first] += 1
        degree[second] += 1
    return dict(degree)

def rank_index(rows, samples, eight):
    expected = {
        (method, pool, "native_" + release, seed, sample)
        for (method) in (METHODS)
        for (pool) in (POOLS)
        for (release) in (RELEASES)
        for (seed) in (SEEDS if (method == "PRODIGY") else (None,))
        for (sample) in (samples)
    }
    output, states, reasons = {}, {}, {}
    for (row) in (rows):
        key = tuple(row["key"])
        if (not eight and key not in expected):
            continue
        require(
            key in expected and key not in output,
            "Unexpected or duplicate native rank group",
        )
        genes, scores = row["genes"], row["scores"]
        status = row.get("native_status", row.get("status"))
        if (not eight and status is None):
            status = (
                "UNAVAILABLE"
                if (genes is None)
                else "SUCCESS"
                if (genes)
                else "SUCCESS_EMPTY"
            )
        if (genes is None):
            require(
                scores is None and status in ("FAILED", "UNAVAILABLE", "MISSING_GROUP"),
                "Unavailable native ranks require an explicit unavailable state",
            )
        else:
            require(
                isinstance(genes, list)
                and isinstance(scores, list)
                and len(genes) == len(scores)
                and len(set(genes)) == len(genes),
                "Malformed native rank vectors",
            )
            require(
                all(
                    isinstance(gene, str) and gene and gene == gene.strip()
                    for (gene) in (genes)
                ),
                "Ranked gene labels must be literal unique identifiers",
            )
            require(
                all(
                    type(value) in (int, float) and math.isfinite(value)
                    for (value) in (scores)
                ),
                "Native rank scores must be finite",
            )
            require(
                all(left >= right for ((left, right)) in (zip(scores, scores[1:]))),
                "Native ranks must retain descending score order",
            )
            require(
                status == ("SUCCESS" if (genes) else "SUCCESS_EMPTY"),
                "Native rank status differs from its full vectors",
            )
        output[key] = (genes, scores)
        states[key] = status
        reasons[key] = row.get("reason")
    require(set(output) == expected, "The complete native rank grid is required")
    return output, states, reasons

def stable_weights(releases):
    require(len(releases) == 3, "Exactly three release selections are required")
    if (any(value is None for (value) in (releases))):
        return None
    genes = set().union(*(set(value) for (value) in (releases)))
    values = {
        gene: min(Fraction(value.get(gene, 0)) for (value) in (releases))
        for (gene) in (sorted(genes))
    }
    return {gene: value for ((gene, value)) in (values.items()) if (value > 0)}

def serialized_weights(weights):
    if (weights is None):
        return [], None
    entries = []
    for (gene) in (sorted(weights)):
        value = Fraction(weights[gene])
        require(
            0 < value <= 1, "Fractional membership is outside its permitted interval"
        )
        entries.append(
            {
                "gene": gene,
                "weight": float(value),
                "weight_numerator": value.numerator,
                "weight_denominator": value.denominator,
            }
        )
    return entries, sum(
        (
            Fraction(entry["weight_numerator"], entry["weight_denominator"])
            for (entry) in (entries)
        ),
        Fraction(0),
    )

def normalise_three(rows):
    output = []
    status_map = {
        "available_nonempty": "AVAILABLE",
        "available_empty": "AVAILABLE",
        "unavailable": "UNAVAILABLE",
    }
    for (raw) in (rows):
        row = dict(
            raw,
            raw_policy = raw["policy"],
            raw_status = raw["status"],
            raw_native_status = raw["native_status"],
        )
        row["policy"] = "V0DEGREE" if (raw["policy"] == "V0_DEGREE") else raw["policy"]
        row["status"] = status_map[raw["status"]]
        for (gene) in (row["genes"]):
            require(
                float(Fraction(gene["weight_numerator"], gene["weight_denominator"]))
                == gene["weight"],
                "Fractional membership encoding differs",
            )
        output.append(row)
    return output

def assemble_groups(
    ranks, states, reasons, eligibility, degree, metadata, samples, eight, identifiers
):
    result = []
    for (sample) in (samples):
        for (method) in (METHODS):
            for (pool) in (POOLS):
                eligible = eligibility["pools"][pool][sample][
                    "native_eligible_literal_gene_labels"
                ]
                common = eligibility["pools"][pool][sample][
                    "common_eligible_literal_gene_labels"
                ]
                require(
                    len(set(eligible)) == len(eligible)
                    and len(set(common)) == len(common)
                    and set(common) <= set(eligible),
                    "Invalid native/common molecular eligibility",
                )
                ordered = sorted(
                    eligible, key = lambda gene: (-degree.get(gene, 0), gene)
                )
                for (seed) in (SEEDS if (method == "PRODIGY") else (None,)):
                    keys = [
                        (method, pool, "native_" + release, seed, sample)
                        for (release) in (RELEASES)
                    ]
                    raw = [ranks[key] for (key) in (keys)]
                    rc = consensus(raw, common)
                    for (budget) in (BUDGETS):
                        fixed = [
                            leading_membership(genes, scores, eligible, budget)
                            for ((genes, scores)) in (raw)
                        ]
                        selections = dict(zip(POLICIES[:3], fixed))
                        selections["RC"] = leading_membership(
                            rc["genes"], rc["scores"], common, budget
                        )
                        selections["STABLE_ONLY_3"] = stable_weights(fixed)
                        selections["V0_DEGREE"] = leading_membership(
                            ordered,
                            [degree.get(gene, 0) for (gene) in (ordered)],
                            eligible,
                            budget,
                        )
                        for (policy) in (POLICIES):
                            weights = selections[policy]
                            entries, mass = serialized_weights(weights)
                            require(
                                mass is None or mass <= budget,
                                "Returned mass exceeds offered budget",
                            )
                            row = {
                                "sample_id": sample,
                                "method": method,
                                "normal_pool": pool,
                                "master_seed": seed,
                                "policy": policy,
                                "k": budget,
                                "support": "native",
                                "genes": entries,
                                "status": "unavailable"
                                if (weights is None)
                                else "available_nonempty"
                                if (entries)
                                else "available_empty",
                                "native_status": "SOURCE_BOUND_UNTRAINED_V0_DEGREE"
                                if (policy == "V0_DEGREE")
                                else "UNAVAILABLE"
                                if (weights is None)
                                else "SUCCESS"
                                if (entries)
                                else "SUCCESS_EMPTY",
                                "source_release_statuses": {
                                    release: states[key]
                                    for ((release, key)) in (zip(RELEASES, keys))
                                },
                                "source_release_reasons": {
                                    release: reasons[key]
                                    for ((release, key)) in (zip(RELEASES, keys))
                                },
                                "reason": "Native failure/unavailability is not empty or zero"
                                if (weights is None)
                                else None,
                                "returned_mass": None
                                if (mass is None)
                                else float(mass),
                                "returned_numerator": None
                                if (mass is None)
                                else mass.numerator,
                                "returned_denominator": None
                                if (mass is None)
                                else mass.denominator,
                                "degree_alias_not_independent_seed_or_method_evidence": policy
                                == "V0_DEGREE",
                                "source_ids": identifiers,
                            }
                            if (eight):
                                row.update(
                                    model_id = sample,
                                    author_alias = metadata["author_aliases"][sample],
                                    returned_mass_numerator = row["returned_numerator"],
                                    returned_mass_denominator = row[
                                        "returned_denominator"
                                    ],
                                    persona_pool_alias_not_independent = method
                                    == "PersonaDrive",
                                    native_cohort_size = 8,
                                )
                            else:
                                row.update(
                                    historical_split = metadata["historical_splits"][
                                        sample
                                    ],
                                    historical_persona_cohort_size = 85
                                    if (method == "PersonaDrive")
                                    else None,
                                )
                            result.append(row)
    require(
        len(result) == (4032 if (eight) else 1512), "Complete policy grid is required"
    )
    order = {value: index for ((index, value)) in (enumerate(POLICIES))}
    result.sort(
        key = lambda row: (
            samples.index(row["sample_id"]),
            METHODS.index(row["method"]),
            POOLS.index(row["normal_pool"]),
            -1 if (row["master_seed"] is None) else row["master_seed"],
            order[row["policy"]],
            BUDGETS.index(row["k"]),
        )
    )
    return result if (eight) else normalise_three(result)

def continuous_selections(config, base):
    require(
        config["cohort"] in ("eight_model", "three_model"),
        "Unknown continuous-validation cohort",
    )
    identifiers = source_ids(
        config, ("ranks", "eligibility", "baseline_graph", "metadata")
    )
    eight = config["cohort"] == "eight_model"
    samples = EIGHT if (eight) else THREE
    metadata = load_json(resolve_path(config["metadata"], base))
    require(
        metadata["native_cohort_size"] == (8 if (eight) else 85),
        "Original complete native cohort size is required",
    )
    if (eight):
        require(
            set(metadata["author_aliases"]) == set(EIGHT)
            and all(
                isinstance(value, str) and value
                for (value) in (metadata["author_aliases"].values())
            ),
            "All eight official source aliases are required",
        )
    else:
        require(
            set(metadata["historical_splits"]) == set(THREE)
            and tuple(metadata["historical_splits"][sample] for (sample) in (THREE))
            == SPLITS,
            "Original three-model exposure splits are required",
        )
    ranks, states, reasons = rank_index(
        read_rows(resolve_path(config["ranks"], base)), samples, eight
    )
    eligibility = load_json(resolve_path(config["eligibility"], base))
    require(
        set(eligibility["pools"]) == set(POOLS),
        "All three normal-pool eligibility records are required",
    )
    with Path(resolve_path(config["baseline_graph"], base)).open(
        encoding = "utf-8-sig", newline = ""
    ) as stream:
        degree = degree_from_rows(csv.DictReader(stream, delimiter = "\t"))
    return assemble_groups(
        ranks,
        states,
        reasons,
        eligibility,
        degree,
        metadata,
        samples,
        eight,
        identifiers,
    )
