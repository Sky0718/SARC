import csv
import gzip

def quoted(value):
    return "'" + str(value).replace("'", "''") + "'"

def check(condition, message):
    if (not condition):
        raise ValueError(message)

def single_row(connection, sql):
    cursor = connection.execute(sql)
    return dict(zip((item[0] for (item) in (cursor.description)), cursor.fetchone()))

def info_rows(path):
    with gzip.open(path, "rt", encoding = "utf-8", newline = "") as handle:
        records = csv.DictReader(handle, delimiter = "\t")
        rows = [
            {"anchor_id": row["#string_protein_id"], "label": row["preferred_name"]}
            for (row) in (records)
        ]
    check(
        (len(rows) == len({row["anchor_id"] for (row) in (rows)})),
        "Duplicate protein metadata identity",
    )
    check(
        all(row["anchor_id"].startswith("9606.") for (row) in (rows)),
        "Human protein identities are required",
    )
    return rows

def copy_query(connection, query, path):
    connection.execute(
        f"COPY ({query}) TO {quoted(path.as_posix())} (FORMAT PARQUET, COMPRESSION ZSTD)"
    )
