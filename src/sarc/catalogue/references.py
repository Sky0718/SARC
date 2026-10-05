import csv
import re
from collections import Counter
from pathlib import Path

SPECS = {
    "COAD_CCLE": ("COAD", "colorectal", "colorectal cancer", "DOID:9256"),
    "LUAD_CCLE": ("LUAD", "lung", "lung cancer", "DOID:1324"),
    "COAD_TCGA": ("COAD", "colorectal", "colorectal cancer", "DOID:9256"),
}

def table(path):
    with Path(path).open(encoding = "utf-8-sig", newline = "") as stream:
        return list(csv.DictReader(stream, delimiter = "\t"))

def literal_states(labels):
    counts = Counter(labels)
    result = {}
    for (label) in (labels):
        reasons = []
        if (not label):
            reasons.append("EMPTY_LITERAL_LABEL")
        if (label != label.strip() or any(
            character.isspace() for (character) in (label)
        )):
            reasons.append("WHITESPACE_LITERAL_LABEL")
        if (counts[label] != 1):
            reasons.append("DUPLICATE_LITERAL_LABEL")
        if (re.fullmatch(
            r"(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\.\d{2}",
            label,
            flags = re.IGNORECASE,
        )):
            reasons.append("UNRESOLVED_DATE_ARTIFACT")
        result[label] = reasons
    return result

def build_catalogues(primary_rows, secondary_rows, known):
    result = {}
    for (context, (cancer, substring, exact_name, disease_id)) in (SPECS.items()):
        selected_primary = {
            row["symbol"]
            for (row) in (primary_rows)
            if (substring in row["primary_site"])
        }
        secondary_by_entrez = {}
        secondary_rows_selected = []
        for (row) in (secondary_rows):
            if (row["cancer_normalized"] != exact_name):
                continue
            if (row["cancer_id"] != disease_id):
                raise ValueError("Frozen CancerMine exact disease identity conflict")
            if (row["role"] not in {"Driver", "Oncogene", "Tumor_Suppressor"}):
                continue
            if (not row["citation_count"].isdigit() or int(row["citation_count"]) <= 0):
                continue
            if (not row["gene_entrez_id"].isdigit() or int(row["gene_entrez_id"]) <= 0):
                raise ValueError("Invalid CancerMine Entrez identity")
            secondary_by_entrez.setdefault(row["gene_entrez_id"], set()).add(
                row["gene_normalized"]
            )
            secondary_rows_selected.append(row)
        if (any(len(names) != 1 for (names) in (secondary_by_entrez.values()))):
            raise ValueError("Multiple CancerMine names for one Entrez identity")
        secondary_names = {
            next(iter(names)) for (names) in (secondary_by_entrez.values())
        }
        if (len(secondary_names) != len(secondary_by_entrez) or any(
            literal_states(sorted(secondary_names)).values()
        )):
            raise ValueError("Unresolved CancerMine literal identity")
        result[context] = {
            "NCG6_primary_all": selected_primary,
            "NCG6_known_subgroup": selected_primary & known,
            "NCG6_candidate_subgroup": selected_primary - known,
            "CancerMine2019_secondary": secondary_names,
        }
    return result

def mutation_sets(path, samples):
    with Path(path).open(encoding = "utf-8-sig", newline = "") as stream:
        reader = csv.reader(stream)
        header = next(reader)
        first = next(reader)
        columns = header if (len(first) == len(header) + 1) else header[1:]
        canonical = [name.replace(".", "-") for (name) in (columns)]
        if (canonical != samples or len(set(canonical)) != len(samples)):
            raise ValueError("Original complete mutation sample order changed")
        rows = [first] + list(reader)
    if (any(len(row) != len(samples) + 1 for (row) in (rows))):
        raise ValueError("Ragged mutation matrix")
    labels = [row[0] for (row) in (rows)]
    states = literal_states(labels)
    mutated = {sample: set() for (sample) in (samples)}
    for (row) in (rows):
        for (sample, text) in (zip(samples, row[1:])):
            value = float(text)
            if (value not in {0.0, 1.0}):
                raise ValueError("Original MUT value is not binary")
            if (value != 0):
                mutated[sample].add(row[0])
    return mutated, states, columns, len(rows)
