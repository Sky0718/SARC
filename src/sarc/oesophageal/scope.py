import csv
from collections import defaultdict
from pathlib import Path
from .assay import ensembl, hgnc
from .contracts import read, require

def make_scope(allocation, crosswalk_path, metadata_path, category_path, table, exclusions):
    molecular = {
        (row["stable_ensembl"], row["hgnc_symbol"])
        for (row) in (read(crosswalk_path))
        if (row["status"] == "ONE_TO_ONE_SHARED")
    }
    symbols = defaultdict(set)
    ensembl_codes = defaultdict(set)
    hgnc_codes = defaultdict(set)
    with Path(metadata_path).open(encoding = "utf-8-sig", newline = "") as stream:
        for (row) in (csv.DictReader(stream)):
            if (row["hgnc_status"] != "current"):
                continue
            code = row["ensembl_gene_id"].split(".")[0]
            symbol = row["hgnc_symbol"]
            identifier = row["hgnc_id"]
            if (code and symbol and identifier):
                symbols[symbol].add((code, identifier))
                ensembl_codes[code].add(symbol)
                hgnc_codes[identifier].add(symbol)
    gene_map = {}
    for (code, symbol) in (sorted(molecular)):
        matches = symbols[symbol]
        if (len(matches) == 1):
            matched_code, identifier = next(iter(matches))
            if (
                matched_code == code
                and len(ensembl_codes[code]) == 1
                and len(hgnc_codes[identifier]) == 1
            ):
                gene_map[symbol] = {"ensembl_id": code, "hgnc_id": identifier}
    needed = {gene for (model) in (allocation["models"]) for (gene) in (model["query"])}
    require(needed <= set(gene_map), "Query lacks a one-to-one source gene identity")
    require(
        table["format"] == "zip_space_csv" and table["field"] == "is_depleted",
        "Original supplied binary Table 5 required",
    )
    categories = category_partitions(read(category_path), gene_map)
    return {
        "schema": "ESCA_EVALUATION_SCOPE_V1",
        "models": [
            {key: model[key] for (key) in (("model_id", "sample_id", "donor_id"))}
            for (model) in (allocation["models"])
        ],
        "gene_map": gene_map,
        "table5_source": table,
        "excluded_sample_library_pairs": exclusions,
        "library_policy": "one_retained_library_per_sample",
        "categories": categories,
        "source_identity_maps": {
            "symbol": {symbol: [symbol] for (symbol) in (sorted(symbols))},
            "ensembl": {
                ensembl(code): sorted(values) for (code, values) in (sorted(ensembl_codes.items()))
            },
            "hgnc": {hgnc(code): sorted(values) for (code, values) in (sorted(hgnc_codes.items()))},
        },
        "ambiguous_source_identifiers": {
            "ensembl": sorted(
                code for (code, values) in (ensembl_codes.items()) if (len(values) != 1)
            ),
            "hgnc": sorted(code for (code, values) in (hgnc_codes.items()) if (len(values) != 1)),
            "symbol": sorted(code for (code, values) in (symbols.items()) if (len(values) != 1)),
        },
    }

def category_partitions(source, genes):
    organoid = set(source["sets"]["organoid"])
    cell_line = set(source["sets"]["cell_line"])
    require(
        len(organoid) == 751 and len(cell_line) == 1121,
        "Original complete core gene sets required",
    )
    categories = {"core": {}, "role": {}}
    for (gene) in (genes):
        if (gene in organoid and gene in cell_line):
            core = "BOTH_CORE"
        elif (gene in organoid):
            core = "ORGANOID_ONLY"
        elif (gene in cell_line):
            core = "CELL_LINE_ONLY"
        else:
            core = "NOT_LISTED_IN_EITHER"
        record = source["roles"].get(gene)
        if (record is None):
            role = "UNLISTED"
        else:
            oncogene = "oncogene" in record["cosmic_tokens"]
            suppressor = "TSG" in record["cosmic_tokens"]
            if (oncogene and suppressor):
                role = "DUAL_ONCOGENE_TSG"
            elif (oncogene):
                role = "ONCOGENE"
            elif (suppressor):
                role = "TSG"
            else:
                role = "OTHER_OR_UNSPECIFIED"
        categories["core"][gene] = core
        categories["role"][gene] = role
    return categories
