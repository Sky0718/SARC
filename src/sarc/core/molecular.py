import csv
import io
import re
import zipfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

DATE_LABEL = re.compile(
    r"(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\.\d{2}\Z", re.IGNORECASE
)

@dataclass
class Matrix:
    genes: tuple
    samples: tuple
    values: np.ndarray

    def select(self, genes, samples):
        row_map = {name: index for ((index, name)) in (enumerate(self.genes))}
        column_map = {name: index for ((index, name)) in (enumerate(self.samples))}
        if (len(row_map) != len(self.genes) or len(column_map) != len(self.samples)):
            raise ValueError("Selection requires unique identities")
        if (len(set(genes)) != len(genes) or len(set(samples)) != len(samples)):
            raise ValueError("Repeated requested identity")
        rows = [row_map[name] for (name) in (genes)]
        columns = [column_map[name] for (name) in (samples)]
        return Matrix(
            tuple(genes), tuple(samples), self.values[np.ix_(rows, columns)].copy()
        )

def valid_gene_labels(genes):
    counts = Counter(genes)
    return {
        gene
        for ((gene, count)) in (counts.items())
        if (
            isinstance(gene, str)
            and gene
            and gene == gene.strip()
            and not any(letter.isspace() for (letter) in (gene))
            and count == 1
            and not DATE_LABEL.fullmatch(gene)
        )
    }

def load_matrix(path, sample_map):
    with Path(path).open("r", encoding = "utf-8-sig", newline = "") as stream:
        reader = csv.reader(stream)
        header, first = next(reader), next(reader)
    implicit = len(first) == len(header) + 1
    original_samples = header if (implicit) else header[1:]
    if (
        not original_samples
        or len(original_samples) != len(set(original_samples))
        or set(sample_map) != set(original_samples)
    ):
        raise ValueError("Incomplete or duplicated source sample mapping")
    samples = tuple(sample_map[name] for (name) in (original_samples))
    if (len(samples) != len(set(samples))):
        raise ValueError("Canonical sample collision")
    with Path(path).open("r", encoding = "utf-8-sig", newline = "") as stream:
        reader = csv.reader(stream)
        next(reader)
        if (any(len(row) != len(samples) + 1 for (row) in (reader))):
            raise ValueError(
                "Source matrix row width differs from complete sample roster"
            )
    frame = pd.read_csv(
        path,
        header = None,
        skiprows = 1,
        names = list(range(len(samples) + 1)),
        dtype = {0: str},
        index_col = False,
        keep_default_na = False,
        float_precision = "round_trip",
        encoding = "utf-8-sig",
    )
    values = frame.iloc[:, 1:].to_numpy(dtype = float)
    if (
        values.shape != (len(frame.index), len(samples))
        or not np.isfinite(values).all()
        or (values < 0).any()
    ):
        raise ValueError("Invalid quantitative input domain")
    return Matrix(tuple(frame.iloc[:, 0]), samples, values)

def remove_unresolved(matrix, binary = False):
    valid = valid_gene_labels(matrix.genes)
    rows = [index for ((index, gene)) in (enumerate(matrix.genes)) if (gene in valid)]
    values = matrix.values[rows, :].copy()
    if (binary and not np.isin(values, [0, 1]).all()):
        raise ValueError("Expected unchanged binary mutation or DEG input")
    excluded = [
        {"source_row": index + 2, "literal_label": gene}
        for ((index, gene)) in (enumerate(matrix.genes))
        if (gene not in valid)
    ]
    return (
        Matrix(
            tuple(matrix.genes[index] for (index) in (rows)), matrix.samples, values
        ),
        excluded,
    )

def expression_scale(matrix, transformation):
    if (transformation not in ("identity", "log2p1")):
        raise ValueError("Unrecognised expression transformation")
    values = (
        matrix.values.copy()
        if (transformation == "identity")
        else np.log2(matrix.values + 1)
    )
    if (not np.isfinite(values).all()):
        raise ValueError("Non-finite transformed expression")
    return Matrix(matrix.genes, matrix.samples, values)

def measured_genes(archive_path, member):
    with zipfile.ZipFile(archive_path) as archive:
        with (
            archive.open(member) as raw,
            io.TextIOWrapper(raw, encoding = "utf-8-sig", newline = "") as stream,
        ):
            reader = csv.reader(stream)
            next(reader)
            return valid_gene_labels([row[0] for (row) in (reader)])

def sample_mapping(columns):
    result = {
        column["original_column"]: column["canonical_sample"] for (column) in (columns)
    }
    if (len(result) != len(columns) or len(set(result.values())) != len(columns)):
        raise ValueError("Sample mapping is not one-to-one")
    return result

def quantitative_inputs(root, contract, observed_genes):
    dawn = contract["DawnRank"]
    prodigy = contract["PRODIGY"]
    raw_tumour = load_matrix(
        Path(root) / dawn["tumour_file"]["path"], sample_mapping(dawn["tumour_columns"])
    )
    raw_normal = load_matrix(
        Path(root) / dawn["normal_file"]["path"], sample_mapping(dawn["normal_columns"])
    )
    raw_counts = load_matrix(
        Path(root) / prodigy["read_counts_file"]["path"],
        sample_mapping(prodigy["columns_in_original_order"]),
    )
    tumour, tumour_excluded = remove_unresolved(raw_tumour)
    normal, normal_excluded = remove_unresolved(raw_normal)
    counts, counts_excluded = remove_unresolved(raw_counts)
    samples = tuple(contract["original_tumour_order"])
    if (set(tumour.samples) != set(samples) or set(tumour.samples) & set(normal.samples)):
        raise ValueError("Incorrect tumour roster or canonical specimen collision")
    tumour = tumour.select(
        [gene for (gene) in (tumour.genes) if (gene in observed_genes)], samples
    )
    normal = normal.select(
        [gene for (gene) in (normal.genes) if (gene in observed_genes)], normal.samples
    )
    tumour = expression_scale(
        tumour, "log2p1" if (contract["cohort"].endswith("_TCGA")) else "identity"
    )
    normal = expression_scale(normal, "log2p1")
    origins = tuple(
        column["sample_origin_for_original_PRODIGY"]
        for (column) in (prodigy["columns_in_original_order"])
    )
    if (len(origins) != len(counts.samples) or set(origins) != {"tumor", "normal"}):
        raise ValueError("Invalid explicit count-column origins")
    count_tumours = {
        sample
        for ((sample, origin)) in (zip(counts.samples, origins))
        if (origin == "tumor")
    }
    if (count_tumours != set(samples)):
        raise ValueError("Count matrix does not retain the entire tumour population")
    original_tumours = {
        item["canonical_sample"]: item["original_column"]
        for (item) in (dawn["tumour_columns"])
    }
    original_normals = {
        item["canonical_sample"]: item["original_column"]
        for (item) in (dawn["normal_columns"])
    }
    normalisation_names = {
        "tumour": tuple(original_tumours[sample] for (sample) in (samples)),
        "normal": tuple(original_normals[sample] for (sample) in (normal.samples)),
    }
    return {
        "tumour": tumour,
        "normal": normal,
        "counts": counts,
        "count_origins": origins,
        "normalisation_names": normalisation_names,
        "excluded": {
            "tumour": tumour_excluded,
            "normal": normal_excluded,
            "counts": counts_excluded,
        },
    }

def dawn_projection(quantitative, mutation, network_nodes):
    tumour, normal = quantitative["tumour"], quantitative["normal"]
    genes = sorted(
        set(tumour.genes) & set(normal.genes) & set(mutation.genes) & set(network_nodes)
    )
    if (not genes):
        raise ValueError("Empty DawnRank input support")
    return {
        "tumour": tumour.select(genes, tumour.samples),
        "normal": normal.select(genes, normal.samples),
        "mutation": mutation.select(genes, tumour.samples),
        "genes": genes,
        "normalisation_names": quantitative["normalisation_names"],
    }

def prodigy_projection(quantitative, mutation, network_nodes):
    counts = quantitative["counts"]
    genes = sorted(set(counts.genes) & set(mutation.genes) & set(network_nodes))
    if (not genes):
        raise ValueError("Empty PRODIGY input support")
    return {
        "counts": counts.select(genes, counts.samples),
        "mutation": mutation.select(genes, mutation.samples),
        "sample_origins": quantitative["count_origins"],
        "genes": genes,
    }

def common_nonrelease_genes(quantitative, mutation):
    return (
        set(mutation.genes)
        & set(quantitative["tumour"].genes)
        & set(quantitative["normal"].genes)
        & set(quantitative["counts"].genes)
    )
