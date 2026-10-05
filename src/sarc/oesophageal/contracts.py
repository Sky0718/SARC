import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

STATES = ("G0H0", "G1H0", "G0H1", "G1H1", "degree_G0", "degree_G1")
SEEDS = [104729, 130363, 155921, 196613, 228017]
CONTRASTS = {
    "Delta": {"G1H1": 1, "G0H0": -1},
    "Theta": {"G1H1": 1, "G0H0": -1, "degree_G1": -1, "degree_G0": 1},
    "graph_H0": {"G1H0": 1, "G0H0": -1},
    "graph_H1": {"G1H1": 1, "G0H1": -1},
    "evidence_G0": {"G0H1": 1, "G0H0": -1},
    "evidence_G1": {"G1H1": 1, "G1H0": -1},
    "interaction": {"G1H1": 1, "G1H0": -1, "G0H1": -1, "G0H0": 1},
}

def require(condition, message):
    if (not condition):
        raise ValueError(message)

def unique_pairs(pairs):
    result = {}
    for (key, value) in (pairs):
        require(key not in result, "Duplicate JSON key: " + str(key))
        result[key] = value
    return result

def read(path):
    with Path(path).open(encoding = "utf-8-sig") as stream:
        return json.load(stream, object_pairs_hook = unique_pairs)

def exact(value):
    require(type(value) in (str, int), "Exact integer or rational string required")
    return Fraction(value)

def resolve(base, path):
    return (Path(base) / path).resolve()

def new_directory(path):
    path = Path(path).resolve()
    require(not path.exists(), "A new output directory is required")
    path.mkdir(parents = True)
    return path

def save(path, value):
    with Path(path).open("x", encoding = "utf-8", newline = "\n") as stream:
        json.dump(value, stream, indent = 2, ensure_ascii = False, allow_nan = False)
        stream.write("\n")

def validate_models(models):
    require(len(models) == 59, "Complete 59-model population required")
    for (field) in (("sample_id", "model_id")):
        require(len({row[field] for (row) in (models)}) == 59, "Unique model identities required")
    require(
        all(isinstance(row["donor_id"], str) and row["donor_id"] for (row) in (models)),
        "Literal donor identities required",
    )
    counts = Counter(row["donor_id"] for (row) in (models))
    require(len(counts) == 58, "Complete 58-donor population required")
    return counts

def validate_scope(scope, allocation):
    require(
        allocation["budget"] == 10 and allocation["seeds"] == SEEDS,
        "Ten slots and five deterministic references required",
    )
    models = allocation["models"]
    donors = validate_models(models)
    keys = ("sample_id", "donor_id", "model_id")
    require(
        [tuple(row[key] for (key) in (keys)) for (row) in (models)]
        == [tuple(row[key] for (key) in (keys)) for (row) in (scope["models"])],
        "Scope and allocation rosters differ",
    )
    mapping = scope["gene_map"]
    require(
        mapping and set(scope["categories"]) == {"core", "role"},
        "Complete gene map and both partitions required",
    )
    for (field) in (("ensembl_id", "hgnc_id")):
        require(all(row.get(field) for (row) in (mapping.values())), "Missing gene identity")
        require(
            len({str(row[field]) for (row) in (mapping.values())}) == len(mapping),
            "Gene mapping is not one-to-one",
        )
    allowed = {
        "core": {"BOTH_CORE", "ORGANOID_ONLY", "CELL_LINE_ONLY", "NOT_LISTED_IN_EITHER"},
        "role": {"ONCOGENE", "TSG", "DUAL_ONCOGENE_TSG", "OTHER_OR_UNSPECIFIED", "UNLISTED"},
    }
    for (partition, categories) in (scope["categories"].items()):
        require(
            set(categories) == set(mapping) and set(categories.values()) <= allowed[partition],
            "Invalid complete category partition",
        )
    for (model) in (models):
        query = model["query"]
        require(
            len(query) == len(set(query)) and set(query) <= set(mapping),
            "Invalid complete query",
        )
        require(
            set(model["states"]) == set(STATES), "Four method and two degree states required"
        )
        for (state) in (STATES):
            record = model["states"][state]
            require(
                record["status"] == "AVAILABLE", "Unavailable predictions cannot be imputed"
            )
            weights = {gene: exact(weight) for (gene, weight) in (record["weights"].items())}
            require(
                all(gene in query and 0 < weight <= 1 for (gene, weight) in (weights.items())),
                "Invalid query allocation",
            )
            require(
                sum(weights.values(), Fraction()) == min(10, len(query)),
                "Incomplete offered-budget allocation",
            )
    table = scope["table5_source"]
    require(
        table["format"] == "zip_space_csv" and table["field"] == "is_depleted",
        "Original supplied binary source required",
    )
    require(
        table["member"] == "supplementary_table_5_revision.csv"
        and table["member_bytes"] == 171694934,
        "Unexpected Table 5 member",
    )
    require(
        scope["library_policy"] == "one_retained_library_per_sample",
        "Original source library rule required",
    )
    return donors

def freeze(allocation, scope):
    donors = validate_scope(scope, allocation)
    requests = []
    coefficients = []
    for (model) in (allocation["models"]):
        sample = model["sample_id"]
        union = sorted({gene for (row) in (model["states"].values()) for (gene) in (row["weights"])})
        requests.extend(
            {"sample_id": sample, "gene": gene, **scope["gene_map"][gene]} for (gene) in (union)
        )
        for (contrast, signs) in (CONTRASTS.items()):
            terms = {}
            for (state, sign) in (signs.items()):
                for (gene, weight) in (model["states"][state]["weights"].items()):
                    terms[gene] = terms.get(gene, Fraction()) + sign * exact(weight)
            for (gene, coefficient) in (sorted(terms.items())):
                if (coefficient):
                    coefficients.append(
                        {
                            "sample_id": sample,
                            "donor_id": model["donor_id"],
                            "contrast": contrast,
                            "gene": gene,
                            "model_coefficient": str(coefficient),
                            "population_coefficient": str(
                                coefficient / (len(donors) * donors[model["donor_id"]])
                            ),
                        }
                    )
    return {
        "donor_ids": sorted(donors),
        "donor_model_counts": dict(donors),
        "contrasts": CONTRASTS,
        "coefficients": coefficients,
        "requests": sorted(requests, key = lambda row: (row["sample_id"], row["gene"])),
        "allocation": allocation,
        "scope": scope,
    }
