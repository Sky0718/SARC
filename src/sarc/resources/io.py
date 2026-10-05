import csv
import gzip
import json
from fractions import Fraction
from pathlib import Path

def check(condition, message):
    if (not condition):
        raise ValueError(message)

def quoted(value):
    return "'" + str(value).replace("'", "''") + "'"

def read_json(path):
    return json.loads(Path(path).read_text(encoding = "utf-8-sig"))

def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents = True, exist_ok = True)
    path.write_text(
        json.dumps(value, indent = 2, ensure_ascii = False, allow_nan = False, default = encode)
        + "\n",
        encoding = "utf-8",
    )

def encode(value):
    if (isinstance(value, Fraction)):
        return {"numerator": value.numerator, "denominator": value.denominator}
    if (isinstance(value, (set, frozenset))):
        return sorted(value)
    if (isinstance(value, Path)):
        return str(value)
    raise TypeError(type(value).__name__)

def read_rows(path):
    path = Path(path)
    if (path.suffix == ".parquet"):
        import pyarrow.parquet as pq

        yield from pq.read_table(path).to_pylist()
    elif (path.suffix == ".json"):
        yield from read_json(path)
    else:
        opener = gzip.open if (path.suffix == ".gz") else open
        with opener(path, "rt", encoding = "utf-8-sig", newline = "") as handle:
            if (".jsonl" in path.name):
                for (line) in (handle):
                    if (line.strip()):
                        yield json.loads(line)
            else:
                yield from csv.DictReader(
                    handle, delimiter = "\t" if (".tsv" in path.name) else ","
                )

def write_parquet(path, rows, schema = None):
    import pyarrow as pa
    import pyarrow.parquet as pq

    path = Path(path)
    path.parent.mkdir(parents = True, exist_ok = True)
    pq.write_table(pa.Table.from_pylist(rows, schema = schema), path, compression = "zstd")

def resolve_paths(value, base):
    if (isinstance(value, list)):
        return [resolve_paths(item, base) for (item) in (value)]
    path = Path(value)
    return path if (path.is_absolute()) else (Path(base) / path).resolve()
