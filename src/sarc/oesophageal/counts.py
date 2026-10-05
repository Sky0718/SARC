import csv
import gzip
import io
import zipfile
from collections import defaultdict
from pathlib import Path
from .contracts import require, resolve, save
from .identity import crosswalk, stable

def integer(value):
    number = int(value)
    require(
        number >= 0 and str(number) == str(value).strip(), "Nonnegative integer counts required"
    )
    return number

def tumour_rows(specification, base):
    path = resolve(base, specification["path"])
    if (specification["format"] == "zip_csv"):
        with zipfile.ZipFile(path) as archive:
            matches = [
                item for (item) in (archive.infolist()) if (item.filename == specification["member"])
            ]
            require(len(matches) == 1, "Unique RNA source member required")
            with (
                archive.open(matches[0]) as raw,
                io.TextIOWrapper(raw, encoding = "utf-8-sig", newline = "") as stream,
            ):
                yield from csv.DictReader(stream)
    elif (specification["format"] == "csv_gz"):
        with gzip.open(path, "rt", encoding = "utf-8-sig", newline = "") as stream:
            yield from csv.DictReader(stream)
    else:
        raise ValueError("Unsupported RNA source format")

def read_tumour(specification, models, base):
    ids = [row["model_id"] for (row) in (models)]
    columns = {model: index for (index, model) in (enumerate(ids))}
    axes = {model: set() for (model) in (ids)}
    identities = defaultdict(set)
    values = {}
    for (row) in (tumour_rows(specification, base)):
        model = row["model_id"]
        if (model not in columns):
            continue
        key = (row["gene_id"], row["ensembl_gene_id"], row["gene_symbol"])
        require(key not in axes[model], "Duplicate tumour gene row")
        require(row["data_source"] == "Sanger", "Unexpected tumour count source")
        axes[model].add(key)
        gene = stable(row["ensembl_gene_id"])
        identities[gene].add(key)
        vector = values.setdefault(key, [None] * len(models))
        vector[columns[model]] = integer(row["htseq_read_count"])
    first = axes[ids[0]]
    require(
        len(first) == 41122 and all(axis == first for (axis) in (axes.values())),
        "Complete common tumour source axes required",
    )
    require(
        all(all(value is not None for (value) in (row)) for (row) in (values.values())),
        "Missing tumour count cell",
    )
    return identities, values

def read_normal(path, attributes_path):
    identities = defaultdict(list)
    values = {}
    with gzip.open(path, "rt", encoding = "utf-8-sig", newline = "") as stream:
        require(next(stream).strip() == "#1.2", "Unexpected normal GCT format")
        dimensions = [int(value) for (value) in (next(stream).split())]
        reader = csv.reader(stream, delimiter = "\t")
        header = next(reader)
        samples = header[2:]
        require(len(samples) == len(set(samples)) == 403, "Complete 403-normal roster required")
        for (row) in (reader):
            require(
                len(row) == len(header) and row[0] not in values,
                "Ragged or repeated normal source row",
            )
            identities[stable(row[0])].append((row[0], row[1]))
            values[row[0]] = [integer(value) for (value) in (row[2:])]
    require(
        dimensions == [len(values), len(samples)] and len(values) == 74628,
        "Normal source dimensions differ",
    )
    attributes = {}
    with Path(attributes_path).open(encoding = "utf-8-sig", newline = "") as stream:
        for (row) in (csv.DictReader(stream, delimiter = "\t")):
            if (row["SAMPID"] in samples):
                require(row["SAMPID"] not in attributes, "Duplicate normal metadata")
                attributes[row["SAMPID"]] = row
    require(set(attributes) == set(samples), "Missing normal sample metadata")
    require(
        all(
            row["SMTSD"] == "Esophagus - Gastroesophageal Junction"
            for (row) in (attributes.values())
        ),
        "Unexpected normal tissue",
    )
    donors = ["-".join(sample.split("-")[:2]) for (sample) in (samples)]
    require(len(set(donors)) == 403, "Normal donor population differs")
    return identities, values, samples, donors

def write_matrix(path, genes, samples, values):
    require(
        all(sum(row[index] for (row) in (values)) > 0 for (index) in (range(len(samples)))),
        "Empty count library",
    )
    with gzip.open(path, "xt", encoding = "utf-8", newline = "") as stream:
        writer = csv.writer(stream, lineterminator = "\n")
        writer.writerow(["gene"] + samples)
        writer.writerows([gene] + row for (gene, row) in (zip(genes, values, strict = True)))

def build_counts(specification, models, base, output):
    tumour_ids, tumour_values = read_tumour(specification["tumour_rna"], models, base)
    normal_ids, normal_values, samples, donors = read_normal(
        resolve(base, specification["normal_gct"]),
        resolve(base, specification["normal_attributes"]),
    )
    require(
        samples == specification["normal_samples"],
        "Source normal sample order differs from the complete roster",
    )
    original = (
        resolve(base, specification["common_genes"]).read_text(encoding = "utf-8").splitlines()
    )
    require(
        len(original) == len(set(original)) == 7399, "Original complete common domain required"
    )
    records, chosen = crosswalk(
        tumour_ids, normal_ids, resolve(base, specification["gene_metadata"]), original
    )
    genes = [row["hgnc_symbol"] for (row) in (chosen)]
    tumour = [tumour_values[tuple(row["cmp_identities"][0])] for (row) in (chosen)]
    normal = [normal_values[row["gtex_identities"][0][0]] for (row) in (chosen)]
    write_matrix(
        output / "tumour_htseq_counts.csv.gz",
        genes,
        [row["sample_id"] for (row) in (models)],
        tumour,
    )
    write_matrix(output / "normal_rnaseqc_counts.csv.gz", genes, samples, normal)
    save(output / "gene_crosswalk.json", records)
    save(
        output / "normal_roster.json",
        [
            {"sample_id": sample, "donor_id": donor}
            for (sample, donor) in (zip(samples, donors, strict = True))
        ],
    )
    return {
        "normal_samples": samples,
        "normal_donors": donors,
        "count_gene_symbols": genes,
        "shared_original_query": [gene for (gene) in (original) if (gene in set(genes))],
        "missing_original_query": sorted(set(original) - set(genes)),
    }
