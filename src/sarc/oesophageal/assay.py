import csv
import io
import math
import zipfile
from collections import Counter
from pathlib import Path
from .contracts import require

HEADERS = ["sample_ID", "gene", "ensembl_id", "HGNC_ID", "library", "is_depleted"]
UNKNOWN = {"ABSENT", "EXPLICIT_NULL", "IDENTITY_CONFLICT", "CONFLICTING_CALLS"}

def hgnc(value):
    if (value is None):
        return ""
    text = str(value).strip()
    if (text.startswith("HGNC:")):
        text = text[5:]
    if (text.endswith(".0") and text[:-2].isdigit()):
        text = text[:-2]
    return text if (text.isdigit()) else ""

def ensembl(value):
    if (value is None):
        return ""
    text = str(value).strip()
    base, separator, version = text.partition(".")
    if (separator and not version.isdigit()):
        return ""
    return base if (base.startswith("ENSG") and base[4:].isdigit()) else ""

def supplied_call(value):
    if (value is None or (
        isinstance(value, str)
        and value.strip() in ("", "NA", "N/A", "NaN", "nan", "NULL", "null", "None")
    )):
        return "EXPLICIT_NULL", None
    require(type(value) is not bool, "Boolean is not a supplied numerical classification")
    if (type(value) is float and not math.isfinite(value)):
        return "EXPLICIT_NULL", None
    if (type(value) in (int, float) and value in (0, 1)):
        return ("MEASURED_POSITIVE", 1) if (value == 1) else ("MEASURED_NEGATIVE", 0)
    if (isinstance(value, str) and value.strip() in ("0", "0.0", "1", "1.0")):
        return (
            ("MEASURED_POSITIVE", 1)
            if (value.strip().startswith("1"))
            else ("MEASURED_NEGATIVE", 0)
        )
    if (isinstance(value, str)):
        try:
            if (not math.isfinite(float(value))):
                return "EXPLICIT_NULL", None
        except ValueError:
            pass
    raise ValueError("Selected is_depleted is neither 0, 1 nor an explicit missing value")

def source_rows(spec):
    path = Path(spec["path"])
    if (spec["format"] == "zip_space_csv"):
        with zipfile.ZipFile(path) as archive:
            members = [item for (item) in (archive.infolist()) if (item.filename == spec["member"])]
            require(
                len(members) == 1 and members[0].file_size == spec["member_bytes"],
                "Exact unique archived Table5 member required",
            )
            with archive.open(members[0]) as raw:
                with io.TextIOWrapper(raw, encoding = "utf-8-sig", newline = "") as stream:
                    reader = csv.DictReader(stream, delimiter = " ", skipinitialspace = True)
                    require(
                        reader.fieldnames == HEADERS, "Exact admitted Table5 fields required"
                    )
                    for (index, row) in (enumerate(reader, start = 2)):
                        require(
                            None not in row and all(name in row for (name) in (HEADERS)),
                            "Malformed archived Table5 row",
                        )
                        yield index, row
    elif (spec["format"] == "whitespace"):
        with path.open(encoding = "utf-8-sig") as stream:
            header = next(stream).split()
            require(header == HEADERS, "Exact whitespace Table5 header required")
            for (index, line) in (enumerate(stream, start = 2)):
                if (not line.strip()):
                    continue
                fields = line.split()
                require(len(fields) == len(header), "Table5 row width differs")
                yield index, dict(zip(header, fields))
    else:
        raise ValueError("Unsupported admitted Table5 format")

def select_rows(rows, requests, scope):
    wanted = {(row["sample_id"], row["gene"]): row for (row) in (requests)}
    require(len(wanted) == len(requests), "Duplicate frozen assay request")
    samples = {row["sample_id"] for (row) in (scope["models"])}
    require(
        all(sample in samples for (sample, gene) in (wanted)),
        "Request outside complete frozen cohort",
    )
    exclusions = {tuple(row) for (row) in (scope["excluded_sample_library_pairs"])}
    observed = {}
    libraries = {sample: set() for (sample) in (samples)}
    by_ensembl = {ensembl(row["ensembl_id"]): gene for (gene, row) in (scope["gene_map"].items())}
    by_hgnc = {hgnc(row["hgnc_id"]): gene for (gene, row) in (scope["gene_map"].items())}
    source_maps = scope.get("source_identity_maps")
    if (source_maps is None):
        source_maps = {
            "symbol": {gene: [gene] for (gene) in (scope["gene_map"])},
            "ensembl": {code: [gene] for (code, gene) in (by_ensembl.items())},
            "hgnc": {code: [gene] for (code, gene) in (by_hgnc.items())},
        }
    require(
        set(source_maps) == {"symbol", "ensembl", "hgnc"},
        "Complete source identity-map types required",
    )
    ambiguous = scope.get("ambiguous_source_identifiers", {})
    ambiguous_ensembl = set(ambiguous.get("ensembl", []))
    ambiguous_hgnc = {hgnc(value) for (value) in (ambiguous.get("hgnc", []))}
    ambiguous_symbol = set(ambiguous.get("symbol", []))
    counts = Counter()
    issues = []
    for (source_row, row) in (rows):
        counts["physically_parsed_rows"] += 1
        sample = row.get("sample_ID")
        if (sample not in samples):
            continue
        counts["frozen_cohort_source_rows"] += 1
        library = row.get("library")
        require(
            isinstance(library, str) and library.strip() == library and library,
            "Missing or nonliteral library identity",
        )
        if ((sample, library) in exclusions):
            counts["source_excluded_library_rows"] += 1
            continue
        libraries[sample].add(library)
        require(
            len(libraries[sample]) <= 1,
            "Multiple retained source libraries for sample " + sample,
        )
        candidates = set()
        candidates.update(source_maps["symbol"].get(row.get("gene"), []))
        candidates.update(source_maps["ensembl"].get(ensembl(row.get("ensembl_id")), []))
        candidates.update(source_maps["hgnc"].get(hgnc(row.get("HGNC_ID")), []))
        relevant = [gene for (gene) in (candidates) if ((sample, gene) in wanted)]
        if (not relevant):
            continue
        counts["selected_source_rows"] += 1
        identity_valid = (
            len(candidates) == 1
            and ensembl(row.get("ensembl_id")) not in ambiguous_ensembl
            and hgnc(row.get("HGNC_ID")) not in ambiguous_hgnc
            and row.get("gene") not in ambiguous_symbol
        )
        if (not identity_valid):
            state, value = "IDENTITY_CONFLICT", None
            issues.append(
                {
                    "source_row": source_row,
                    "sample_id": sample,
                    "genes": sorted(candidates),
                    "reason": "SOURCE_IDENTIFIER_DISAGREEMENT",
                    "ensembl_id": row.get("ensembl_id"),
                    "hgnc_id": row.get("HGNC_ID"),
                    "library": library,
                }
            )
        else:
            state, value = supplied_call(row.get("is_depleted"))
        representation = [row.get(name) for (name) in (HEADERS)]
        for (gene) in (relevant):
            key = sample, gene
            current = {
                "state": state,
                "value": value,
                "library": library,
                "source_rows": [source_row],
                "multiplicity": 1,
                "raw_representations": [representation],
            }
            if (key not in observed):
                observed[key] = current
            else:
                previous = observed[key]
                previous["multiplicity"] += 1
                previous["source_rows"].append(source_row)
                if (representation not in previous["raw_representations"]):
                    previous["raw_representations"].append(representation)
                    previous["state"], previous["value"] = "CONFLICTING_CALLS", None
                    issues.append(
                        {
                            "source_row": source_row,
                            "sample_id": sample,
                            "gene": gene,
                            "reason": "NONIDENTICAL_MAPPED_DUPLICATE",
                            "library": library,
                        }
                    )
                else:
                    counts["identical_selected_duplicates_collapsed"] += 1
    require(
        all(len(values) == 1 for (values) in (libraries.values())),
        "Every frozen sample must have exactly one retained source library",
    )
    selected = []
    for (key, request) in (sorted(wanted.items())):
        call = observed.get(
            key,
            {
                "state": "ABSENT",
                "value": None,
                "library": next(iter(libraries[key[0]]), None),
                "source_rows": [],
                "multiplicity": 0,
                "raw_representations": [],
            },
        )
        selected.append({"sample_id": key[0], "gene": key[1], **call})
    return {
        "calls": selected,
        "issues": issues,
        "counts": dict(counts),
        "retained_libraries": {key: sorted(value) for (key, value) in (libraries.items())},
        "state_counts": dict(Counter(row["state"] for (row) in (selected))),
        "exposure": "Source rows are physically parsed to locate frozen identities and enforce source library rules. Only frozen requested is_depleted cells are admitted to endpoint evaluation; no unrequested assay endpoints are exported or used.",
    }
