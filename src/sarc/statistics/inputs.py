import csv
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

from .endpoint import as_fraction, require
from .paired import CELLS, METHODS, SEEDS, combine

POPULATIONS = {"CRC": (85, 83), "ESCA": (59, 58)}
CORE = ("BOTH_CORE", "ORGANOID_ONLY", "CELL_LINE_ONLY", "NOT_LISTED_IN_EITHER")
ROLES = ("ONCOGENE", "TSG", "DUAL_ONCOGENE_TSG", "OTHER_OR_UNSPECIFIED", "UNLISTED")
CATEGORIES = {"all": ("ALL",), "core": CORE, "role": ROLES}

def read_json(path):
    return json.loads(Path(path).read_text(encoding = "utf-8-sig"))

def serial(value):
    if (isinstance(value, Fraction)):
        return str(value)
    if (isinstance(value, dict)):
        return {str(key): serial(item) for (key, item) in (value.items())}
    if (isinstance(value, (list, tuple))):
        return [serial(item) for (item) in (value)]
    return value

def save_json(path, value):
    with Path(path).open("x", encoding = "utf-8") as stream:
        json.dump(serial(value), stream, indent = 2, ensure_ascii = False, allow_nan = False)
        stream.write("\n")

def read_csv(path):
    with Path(path).open(encoding = "utf-8-sig", newline = "") as stream:
        return list(csv.DictReader(stream))

def save_csv(path, rows, fields = None):
    if (fields is None):
        require(bool(rows), "CSV rows are empty")
        fields = list(dict.fromkeys(key for (row) in (rows) for (key) in (row)))
    with Path(path).open("x", encoding = "utf-8", newline = "") as stream:
        writer = csv.DictWriter(stream, fieldnames = fields)
        writer.writeheader()
        for (row) in (rows):
            writer.writerow(
                {
                    key: json.dumps(serial(value))
                    if (isinstance(value, (list, tuple, dict)))
                    else serial(value)
                    for (key, value) in (row.items())
                }
            )

def resolve(base, value):
    path = Path(value)
    return path if (path.is_absolute()) else base / path

def validate_models(models, context):
    expected_models, expected_donors = POPULATIONS[context]
    require(
        (len(models) == len({row["model_id"] for (row) in (models)}) == expected_models),
        "Full distinct model population required",
    )
    require(
        (len({row["sample_id"] for (row) in (models)}) == expected_models),
        "Distinct sample identities required",
    )
    require(
        (len({row["donor_id"] for (row) in (models)}) == expected_donors),
        "Full donor population required",
    )
    require(
        all(row["model_id"] and row["sample_id"] and row["donor_id"] for (row) in (models)),
        "Empty model identity",
    )
    return Counter(row["donor_id"] for (row) in (models))

def weights(source, query, capacity = None):
    pairs = list(source.items()) if (isinstance(source, dict)) else source
    require((len(pairs) == len({gene for (gene, value) in (pairs)})), "Duplicate gene allocation")
    output = {gene: as_fraction(value) for (gene, value) in (pairs)}
    require((set(output) <= set(query)), "Allocation outside full model query")
    require(
        all(0 <= value <= 1 for (value) in (output.values())), "Allocation outside unit capacity"
    )
    expected = Fraction(min(10, len(query))) if (capacity is None) else as_fraction(capacity)
    require((0 <= expected <= min(10, len(query))), "Invalid allocation capacity")
    require(
        (sum(output.values(), Fraction()) == expected),
        "Allocation mass differs from admitted capacity",
    )
    return {gene: value for (gene, value) in (output.items()) if (value)}

def method_weights(row, query, method):
    status = row.get("status", "AVAILABLE")
    require((status in ("AVAILABLE", "SUCCESS", "SUCCESS_EMPTY")), "Unavailable allocation")
    native_status = row.get("native_status", status)
    require(
        (native_status in ("AVAILABLE", "SUCCESS", "SUCCESS_EMPTY")), "Unavailable native state"
    )
    default_menu = (
        [] if (native_status == "SUCCESS_EMPTY") else query if (method == "DawnRank") else None
    )
    menu = row.get("menu", default_menu)
    require((menu is not None), "Complete emitted menu required for PRODIGY")
    require((len(menu) == len(set(menu))), "Duplicate emitted menu gene")
    eligible = sorted(set(menu) & set(query))
    capacity = min(10, len(eligible))
    if ("capacity" in row):
        require(
            (as_fraction(row["capacity"]) == capacity), "Recorded emitted-menu capacity differs"
        )
    if (native_status == "SUCCESS_EMPTY"):
        require((not menu and capacity == 0), "Native empty state has an emitted menu")
    return weights(row["weights"], eligible, capacity)

def load_calls(paths):
    calls = {}
    states = {
        "MEASURED",
        "MEASURED_POSITIVE",
        "MEASURED_NEGATIVE",
        "ABSENT",
        "EXPLICIT_NULL",
        "IDENTITY_CONFLICT",
        "CONFLICTING_CALLS",
    }
    for (path) in (paths):
        for (row) in (read_json(path)["calls"]):
            key = (row["sample_id"], row["gene"])
            state, value = row["state"], row["value"]
            require((state in states), "Unsupported label state")
            measured = state in {"MEASURED", "MEASURED_POSITIVE", "MEASURED_NEGATIVE"}
            require(
                (type(value) is int and value in (0, 1)) if (measured) else (value is None),
                "Label value inconsistent with state",
            )
            if (state == "MEASURED_POSITIVE"):
                require((value == 1), "Positive call mismatch")
            if (state == "MEASURED_NEGATIVE"):
                require((value == 0), "Zero call mismatch")
            normal = ("MEASURED" if (measured) else state, value)
            require((key not in calls or calls[key] == normal), "Conflicting selected labels")
            calls[key] = normal
    return calls

def require_support(models, policies, calls):
    for (model) in (models):
        sample, mid = model["sample_id"], model["model_id"]
        for (states) in (policies[mid].values()):
            for (vector) in (states.values()):
                require(
                    all((sample, gene) in calls for (gene) in (vector)),
                    "Missing selected label record",
                )

def load_crc(base, config):
    degree = read_json(resolve(base, config["degree_allocations"]))
    models = degree["models"]
    counts = validate_models(models, "CRC")
    by_id = {row["model_id"]: row for (row) in (models)}
    degrees, allocations = {}, {}
    degree_rows = {}
    for (row) in (degree.get("allocations", [])):
        if (row["budget"] != 10):
            continue
        key = (row["model_id"], row.get("state", row.get("cell")))
        require((key[0] in by_id and key[1] in ("G0", "G1")), "Unknown degree model or state")
        require((key not in degree_rows), "Duplicate degree allocation")
        degree_rows[key] = row["weights"]
    for (model) in (models):
        mid = model["model_id"]
        query = model.get("query", model.get("eligible_genes"))
        require(
            (query is not None and len(query) == len(set(query))),
            "Unique complete query required",
        )
        model["query"] = query
        if ("alpha" in model):
            require(
                (as_fraction(model["alpha"]) == Fraction(1, 83 * counts[model["donor_id"]])),
                "Equal-donor weighting mismatch",
            )
        source = (
            model["full_degree"]
            if ("full_degree" in model)
            else {graph: degree_rows[(mid, graph)] for (graph) in (("G0", "G1"))}
        )
        degrees[mid] = {graph: weights(source[graph], query) for (graph) in (("G0", "G1"))}
    for (source) in (config["allocations"]):
        document = read_json(resolve(base, source["path"]))
        method = source["method"]
        require((method in METHODS), "Unsupported colorectal method")
        source_models = {row["model_id"]: row for (row) in (document["models"])}
        require((len(source_models) == len(document["models"])), "Duplicate source model")
        require((set(source_models) == set(by_id)), "Allocation model membership mismatch")
        for (mid, model) in (source_models.items()):
            require(
                all(model[key] == by_id[mid][key] for (key) in (("sample_id", "donor_id"))),
                "Allocation model identity mismatch",
            )
        records = document["allocations"]
        for (row) in (records):
            if (row.get("method", method) != method or row["budget"] != 10):
                continue
            mid = row["model_id"]
            cell = row.get("cell", row.get("state")).replace("_", "")
            require((mid in by_id and cell in CELLS), "Unsupported allocation cell")
            require(
                (row.get("numerically_certified", True) is True),
                "Uncertified allocation cutoff",
            )
            seed = str(row["seed"]) if (method == "PRODIGY") else "deterministic"
            require(
                (method != "PRODIGY" or seed in {str(value) for (value) in (SEEDS)}),
                "Unexpected PRODIGY seed",
            )
            vector = method_weights(row, by_id[mid]["query"], method)
            key = (method, mid, cell, seed)
            require(
                (key not in allocations or method == "DawnRank" and allocations[key] == vector),
                "Duplicate or inconsistent allocation",
            )
            allocations[key] = vector
    means, policies = {}, {}
    for (model) in (models):
        mid = model["model_id"]
        policies[mid] = {"degree": degrees[mid]}
        for (method) in (METHODS):
            seeds = SEEDS if (method == "PRODIGY") else ("deterministic",)
            policies[mid][method] = {}
            for (cell) in (CELLS):
                keys = [(method, mid, cell, str(seed)) for (seed) in (seeds)]
                require(
                    all(key in allocations for (key) in (keys)),
                    "Incomplete method-cell-seed population",
                )
                vector = combine([(Fraction(1, len(seeds)), allocations[key]) for (key) in (keys)])
                means[(method, mid, cell)] = vector
                policies[mid][method][cell] = vector
    calls = load_calls([resolve(base, path) for (path) in (config["labels"])])
    require_support(models, policies, calls)
    return models, allocations, means, degrees, calls, policies

def load_esca(base, config):
    source = read_json(resolve(base, config["allocations"]))
    models = source["models"]
    validate_models(models, "ESCA")
    policies = {}
    for (model) in (models):
        query = model["query"]
        require((len(query) == len(set(query))), "Duplicate query gene")
        states = model["states"]
        require(
            (set(states) == set(CELLS) | {"degree_G0", "degree_G1"}),
            "Complete six-state allocation required",
        )
        require(
            all(
                state.get("status", "AVAILABLE") in ("AVAILABLE", "SUCCESS", "SUCCESS_EMPTY")
                for (state) in (states.values())
            ),
            "Unavailable oesophageal allocation",
        )
        policies[model["model_id"]] = {
            "DawnRank": {
                cell: method_weights(states[cell], query, "DawnRank") for (cell) in (CELLS)
            },
            "degree": {
                graph: weights(states["degree_" + graph]["weights"], query)
                for (graph) in (("G0", "G1"))
            },
        }
    calls = load_calls([resolve(base, path) for (path) in (config["labels"])])
    require_support(models, policies, calls)
    return models, policies, calls

def load_categories(path, models, context):
    source = read_json(path)
    genes = {gene for (model) in (models) for (gene) in (model["query"])}
    if (context == "ESCA"):
        mappings = source["categories"]
    else:
        organoid, cell_line = set(source["sets"]["organoid"]), set(source["sets"]["cell_line"])
        mappings = {"core": {}, "role": {}}
        for (gene) in (genes):
            org, cell = gene in organoid, gene in cell_line
            mappings["core"][gene] = (
                "BOTH_CORE"
                if (org and cell)
                else "ORGANOID_ONLY"
                if (org)
                else "CELL_LINE_ONLY"
                if (cell)
                else "NOT_LISTED_IN_EITHER"
            )
            record = source["roles"].get(gene)
            tokens = set(record["cosmic_tokens"]) if (record is not None) else set()
            onc, tsg = "oncogene" in tokens, "TSG" in tokens
            mappings["role"][gene] = (
                "UNLISTED"
                if (record is None)
                else "DUAL_ONCOGENE_TSG"
                if (onc and tsg)
                else "ONCOGENE"
                if (onc)
                else "TSG"
                if (tsg)
                else "OTHER_OR_UNSPECIFIED"
            )
    for (partition) in (("core", "role")):
        require((genes <= set(mappings[partition])), "Incomplete category membership")
        require(
            all(mappings[partition][gene] in CATEGORIES[partition] for (gene) in (genes)),
            "Unknown category",
        )
    return mappings

def load_draws(path, donors):
    document = (
        np.load(path, allow_pickle = False) if (Path(path).suffix == ".npy") else read_json(path)
    )
    require((len(document) == 10000), "All 10000 saved paired donor draws required")
    array = np.asarray(document)
    require(
        (array.shape == (10000, len(donors)) and np.issubdtype(array.dtype, np.integer)),
        "Incorrect donor draw dimensions or type",
    )
    require(bool(np.all((array >= 0) & (array < len(donors)))), "Donor draw index out of range")
    return array.tolist()
