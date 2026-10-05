import csv
import io

CLASS_SETS = {
    "PB": {"PB"},
    "PW": {"PW"},
    "DB": {"DB"},
    "DW": {"DW"},
    "PB_PW": {"PB", "PW"},
    "ALL": {"PB", "PW", "DB", "DW"},
}
SEEDS = [104729, 130363, 155921, 196613, 228017]

def edges(path):
    with path.open(encoding = "utf-8", newline = "") as stream:
        reader = csv.reader(stream, delimiter = "\t")
        if (next(reader) != ["gene1", "gene2", "combined_score"]):
            raise ValueError("Unexpected edge header")
        result = [(row[0], row[1], int(row[2])) for (row) in (reader)]
    if (any(
        (not a or a >= b or (not 700 < score <= 1000) for ((a, b, score)) in (result))
    )):
        raise ValueError("Invalid canonical edge")
    if (len({row[:2] for (row) in (result)}) != len(result)):
        raise ValueError("Duplicate canonical edge")
    return result

def mutation_summary(path, samples):
    counts = {}
    totals = [0] * len(samples)
    with path.open(encoding = "utf-8-sig", newline = "") as stream:
        reader = csv.reader(stream)
        header = next(reader)
        if (header[1:] != samples):
            raise ValueError(
                "Mutation sample order differs from accepted complete cohort"
            )
        for (row) in (reader):
            if (len(row) != len(samples) + 1 or not row[0] or row[0] in counts):
                raise ValueError("Invalid mutation matrix identity")
            values = [float(value) for (value) in (row[1:])]
            if (any((value not in (0.0, 1.0) for (value) in (values)))):
                raise ValueError("Mutation input is not binary")
            counts[row[0]] = sum((value == 1.0 for (value) in (values)))
            totals = [total + int(value) for ((total, value)) in (zip(totals, values))]
    return (counts, totals)

def classify(edge, anchors, modules):
    (a, b) = edge
    return ("P" if (a in anchors or b in anchors) else "D") + (
        "W" if (modules[a] == modules[b]) else "B"
    )

def intervene(baseline, other, class_map, selected):
    baseline_scores = {row[:2]: row[2] for (row) in (baseline)}
    other_scores = {row[:2]: row[2] for (row) in (other)}
    removed = {
        edge
        for (edge) in (baseline_scores.keys() - other_scores.keys())
        if (class_map[edge] in selected)
    }
    added = {
        edge
        for (edge) in (other_scores.keys() - baseline_scores.keys())
        if (class_map[edge] in selected)
    }
    result = [row for (row) in (baseline) if (row[:2] not in removed)]
    result.extend(((*edge, other_scores[edge]) for (edge) in (sorted(added))))
    if (len(result) != len(baseline) - len(removed) + len(added)):
        raise ValueError("Wrong intervention size")
    if (len({row[:2] for (row) in (result)}) != len(result)):
        raise ValueError("Duplicate intervention edge")
    return (result, {"removed": len(removed), "added": len(added)})

def graph_bytes(rows):
    stream = io.StringIO(newline = "")
    writer = csv.writer(stream, delimiter = "\t", lineterminator = "\n")
    writer.writerow(["gene1", "gene2", "combined_score"])
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")
