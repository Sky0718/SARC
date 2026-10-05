import csv
import json
import math
from pathlib import Path

SEEDS = (104729, 130363, 155921, 196613, 228017)

def load_json(path):
    return json.loads(Path(path).read_text(encoding = "utf-8-sig"))

def save_json(path, value):
    with Path(path).open("x", encoding = "utf-8", newline = "\n") as stream:
        json.dump(value, stream, ensure_ascii = False, indent = 2, allow_nan = False)
        stream.write("\n")

def resolve_path(value, base):
    path = Path(value).expanduser()
    return str((base / path).resolve() if (not path.is_absolute()) else path.resolve())

def unique_axis(values, label):
    if (
        (not values)
        or len(values) != len(set(values))
        or any(not value for (value) in (values))
    ):
        raise ValueError("Expected a non-empty unique " + label + " axis")
    return values

def read_matrix(path):
    with Path(path).open(newline = "", encoding = "utf-8-sig") as stream:
        reader = csv.reader(stream)
        header = next(reader)
        samples = unique_axis(header[1:], "sample")
        rows = {}
        for (row) in (reader):
            if (len(row) != len(header) or not row[0] or row[0] in rows):
                raise ValueError("Invalid or duplicated matrix row")
            numbers = [float(value) for (value) in (row[1:])]
            if (any(not math.isfinite(value) for (value) in (numbers))):
                raise ValueError("Non-finite matrix input")
            rows[row[0]] = row[1:]
    unique_axis(list(rows), "gene")
    return header[0], samples, rows

def write_matrix(path, first_header, samples, genes, rows):
    with Path(path).open("x", newline = "", encoding = "utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow([first_header] + samples)
        writer.writerows([gene] + rows[gene] for (gene) in (genes))

def validate_seeds(config):
    scope = config.get("parameter_scope", "primary")
    if (scope not in ("primary", "catalogue_secondary")):
        raise ValueError("Unknown parameter scope")
    method = config["method"]
    expected = SEEDS
    if (scope == "catalogue_secondary"):
        if (
            config["mode"] != "native"
            or config["context"] not in ("COAD_CCLE", "LUAD_CCLE")
            or config["expected_models"] != 36
            or config["network_id"] not in ("native_11_0", "native_11_5", "native_12_0")
        ):
            raise ValueError(
                "Secondary parameters require the original development native cells"
            )
        if (method == "PRODIGY"):
            if (config["alpha"] not in (0.01, 0.1)):
                raise ValueError("Secondary PRODIGY alpha must be 0.01 or 0.1")
            expected = (104729,)
        elif (method == "DawnRank"):
            if (config.get("mu") not in (1, 10)):
                raise ValueError("Secondary DawnRank mu must be 1 or 10")
        else:
            raise ValueError("No secondary grid exists for this method")
    elif (method == "PRODIGY" and config["alpha"] != 0.05):
        raise ValueError("Primary PRODIGY alpha must be 0.05")
    elif (method == "DawnRank" and config.get("mu", 3) != 3):
        raise ValueError("Primary DawnRank mu must be 3")
    if (tuple(config["seeds"]) != expected):
        raise ValueError("The exact parameter-scope seed order is required")
    return expected
