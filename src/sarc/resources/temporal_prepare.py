import pyarrow.parquet as pq

def literal(value):
    return "'" + str(value).replace("\\", "/").replace("'", "''") + "'"

def parquet_source(paths):
    return "read_parquet([" + ",".join((literal(path) for (path) in (paths))) + "])"

def read_entities(paths, columns):
    records = {}
    for (path) in (paths):
        for (row) in (pq.read_table(
            path, columns = list(columns), use_threads = False
        ).to_pylist()):
            if (row["id"] in records):
                raise ValueError("Duplicate entity identifier")
            records[row["id"]] = row
    return records
