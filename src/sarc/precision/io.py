import csv
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

def require(condition, message):
    if (not condition):
        raise ValueError(message)

def unique_pairs(pairs):
    result = {}
    for (key, value) in (pairs):
        require(key not in result, "Duplicate JSON key: " + key)
        result[key] = value
    return result

def reject_nonfinite(value):
    raise ValueError("Nonfinite JSON value: " + value)

def read_json(path):
    with Path(path).open(encoding = "utf-8-sig") as stream:
        return json.load(
            stream, object_pairs_hook = unique_pairs, parse_constant = reject_nonfinite
        )

def encode(value):
    if (isinstance(value, Fraction)):
        return str(value)
    if (isinstance(value, np.generic)):
        return value.item()
    raise TypeError(type(value).__name__)

def write_json(path, value):
    with Path(path).open("x", encoding = "utf-8", newline = "\n") as stream:
        json.dump(value, stream, indent = 2, allow_nan = False, default = encode)
        stream.write("\n")

def read_csv(path):
    with Path(path).open(encoding = "utf-8-sig", newline = "") as stream:
        return list(csv.DictReader(stream))

def write_csv(path, records):
    require(bool(records), "Empty output table")
    with Path(path).open("x", encoding = "utf-8", newline = "") as stream:
        writer = csv.DictWriter(stream, fieldnames = list(records[0]))
        writer.writeheader()
        writer.writerows(records)

def resolve(base, value):
    return (Path(base) / value).resolve()

def binary_matrix(path, rows, columns):
    require(Path(path).stat().st_size == rows * columns * 8, "Binary64 matrix size differs")
    values = np.fromfile(path, dtype = "<f8").reshape((rows, columns), order = "F")
    require(np.isfinite(values).all(), "Nonfinite binary64 values")
    return values

def read_axes(directory):
    directory = Path(directory)
    genes = (directory / "genes.txt").read_text(encoding = "utf-8").splitlines()
    samples = (directory / "models.txt").read_text(encoding = "utf-8").splitlines()
    require(bool(genes) and len(genes) == len(set(genes)), "Invalid gene axis")
    require(bool(samples) and len(samples) == len(set(samples)), "Invalid sample axis")
    return genes, samples

def allocate(values, budget):
    require(type(budget) is int and budget >= 0, "Invalid budget")
    require(
        all(isinstance(value, Fraction) for (value) in (values.values())),
        "Exact rational scores required",
    )
    ordered = sorted(values, key = lambda gene: (-values[gene], gene))
    if (not ordered or budget == 0):
        return {}, None, None
    if (len(ordered) <= budget):
        return {gene: Fraction(1) for (gene) in (ordered)}, None, values[ordered[-1]]
    threshold = values[ordered[budget - 1]]
    above = [gene for (gene) in (ordered) if (values[gene] > threshold)]
    tied = [gene for (gene) in (ordered) if (values[gene] == threshold)]
    weights = {gene: Fraction(1) for (gene) in (above)}
    weights.update({gene: Fraction(budget - len(above), len(tied)) for (gene) in (tied)})
    return weights, threshold - values[ordered[budget]], threshold

def exact_score(value):
    require(math.isfinite(float(value)), "Nonfinite score")
    return Fraction.from_float(float(value))

def roster_rows(path, samples, expected_models, expected_donors):
    document = read_json(path)
    rows = document["models"] if (isinstance(document, dict)) else document
    required = {"sample_id", "model_id", "donor_id"}
    require(all(required <= set(row) for (row) in (rows)), "Incomplete model roster")
    require(len(rows) == expected_models, "Model population differs")
    require(len({row["sample_id"] for (row) in (rows)}) == len(rows), "Duplicate samples")
    require(len({row["model_id"] for (row) in (rows)}) == len(rows), "Duplicate model identities")
    require(
        len({row["donor_id"] for (row) in (rows)}) == expected_donors, "Donor population differs"
    )
    require(
        set(samples) == {row["sample_id"] for (row) in (rows)},
        "Prepared sample axis and roster differ",
    )
    by_sample = {row["sample_id"]: row for (row) in (rows)}
    return [
        {key: by_sample[sample][key] for (key) in (("sample_id", "model_id", "donor_id"))}
        for (sample) in (samples)
    ]
