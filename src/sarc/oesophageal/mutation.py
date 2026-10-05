import csv
import gzip
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from .contracts import require, resolve, save

def mapping(path):
    by_gene, by_symbol = defaultdict(set), defaultdict(set)
    with Path(path).open(encoding = "utf-8-sig", newline = "") as stream:
        for (row) in (csv.DictReader(stream)):
            gene = re.sub(r"\.\d+$", "", row["ensembl_gene_id"])
            symbol = row["hgnc_symbol"]
            if (gene and symbol and (row["hgnc_status"] == "current")):
                by_gene[gene].add(symbol)
                by_symbol[symbol].add(gene)
    unique = {
        gene: next(iter(symbols))
        for (gene, symbols) in (by_gene.items())
        if ((len(symbols) == 1) and (len(by_symbol[next(iter(symbols))]) == 1))
    }
    return unique, by_symbol

def parse(line, samples, normal):
    row = line.rstrip("\r\n").split("\t")
    if ((len(row) != (9 + len(samples))) or (len(samples) != 2)):
        raise ValueError("Unexpected source column count")
    if (
        (not re.fullmatch("[ACGTN]+", row[3]))
        or (not re.fullmatch("[ACGTN]+", row[4]))
        or (row[3] == row[4])
        or (not row[1].isdigit())
        or (int(row[1]) < 1)
    ):
        raise ValueError(
            "Unanticipated allele or coordinate representation: " + "|".join(row[:5])
        )
    info = {}
    for (part) in (row[7].split(";")):
        if (part == "."):
            continue
        key, value = part.split("=", 1) if ("=" in part) else (part, True)
        if (key in info):
            raise ValueError("Duplicate INFO key")
        info[key] = value
    fmt = row[8].split(":")
    if (len(fmt) != len(set(fmt))):
        raise ValueError("Duplicate FORMAT key")
    tumour = next(sample for (sample) in (samples) if (sample != normal))
    result = {"row": row, "info": info, "tumour": tumour, "normal": normal}
    for (role, name) in ((("organoid", tumour), ("normal", normal))):
        raw = row[9 + samples.index(name)]
        fields = raw.split(":")
        if ((raw != ".") and (len(fields) != len(fmt))):
            raise ValueError("Unanticipated incomplete FORMAT sample representation")
        result[role + "_raw"] = raw
        result[role + "_format"] = dict(zip(fmt, fields)) if (raw != ".") else {}
    return result

def classify(record, rule, reverse, symbol_candidates):
    row, info = record["row"], record["info"]
    data = record["organoid_format"]
    reasons = []
    if (row[6] not in rule["global_FILTER_exact_allowlist"]):
        reasons.append("GLOBAL_FILTER_NOT_ADMITTED")
    if (data.get("OFS") != "PASS"):
        reasons.append(
            "ORGANOID_OFS_MISSING"
            if (data.get("OFS") in (None, "", ".", "NA"))
            else "ORGANOID_OFS_FAILED"
        )
    if (("gDRV" in info) or ("GRM" in row[6].split(";"))):
        reasons.append("SOURCE_GERMLINE")
    if ("NPGL" in info):
        reasons.append("POPULATION_PANEL_NPGL")
    try:
        vaf = float(data.get("VAF", "nan"))
    except (TypeError, ValueError):
        vaf = float("nan")
    if (not math.isfinite(vaf)):
        reasons.append("ORGANOID_VAF_UNKNOWN")
    elif ((vaf < 0) or (vaf > 1)):
        raise ValueError("Unanticipated VAF outside probability domain")
    elif (vaf <= 0.05):
        reasons.append("ORGANOID_VAF_AT_OR_BELOW_THRESHOLD")
    effect = info.get("VC", "")
    included = (effect in rule["included_exact_short_effects"]) or bool(
        set(effect.split(",")) & set(rule["included_explicit_verbose_SO_tokens"])
    )
    if (not included):
        if (not effect):
            reasons.append("CONSEQUENCE_MISSING")
        elif ((effect in rule["excluded_short_effects"]) or (
            ("SO:0001624:3_prime_UTR_variant" in effect)
            and ("SO:0001629:splice_site_variant" in effect)
        )):
            reasons.append("CONSEQUENCE_EXCLUDED")
        else:
            reasons.append("CONSEQUENCE_UNKNOWN")
    vd = info.get("VD", "")
    parts = vd.split("|") if (vd) else []
    symbol = parts[0] if (parts) else ""
    transcript = parts[1] if (len(parts) > 1) else ""
    gene = reverse.get(symbol, "")
    if (not vd):
        reasons.append("DEFAULT_ANNOTATION_MISSING")
    elif (not gene):
        reasons.append(
            "GENE_MAPPING_AMBIGUOUS"
            if (symbol in symbol_candidates)
            else "GENE_MAPPING_UNRESOLVED"
        )
    return {
        "eligible": (not reasons),
        "reasons": reasons,
        "gene": gene,
        "symbol": symbol,
        "transcript": transcript,
        "effect": effect,
        "VD": vd,
    }

def header(path):
    samples = None
    normal = []
    with gzip.open(path, "rt", encoding = "utf-8") as stream:
        for (line) in (stream):
            if (line.startswith("##SAMPLE=<ID=NORMAL,")):
                match = re.search(r"SampleName=([^>]+)", line)
                require(match is not None, "Missing normal sample definition")
                normal.append(match.group(1))
            if (line.startswith("#CHROM")):
                samples = line.rstrip().split("\t")[9:]
                break
    require(
        samples is not None and len(samples) == len(set(samples)) == 2,
        "Exactly two source sample columns required",
    )
    require(len(normal) == 1 and normal[0] in samples, "Unique designated normal required")
    return samples, normal[0]

def build_mutations(specifications, models, metadata, rule, base, output):
    mapping_values, candidates = mapping(metadata)
    reverse = {symbol: gene for (gene, symbol) in (mapping_values.items())}
    ids = [model["model_id"] for (model) in (models)]
    pairs = [(row["model_id"], row["kind"]) for (row) in (specifications)]
    require(
        len(pairs) == len(set(pairs)) == 118, "Complete unique SNP and indel sources required"
    )
    require(
        set(pairs) == {(model, kind) for (model) in (ids) for (kind) in (("snp", "indel"))},
        "Mutation source population differs",
    )
    eligible = {model: set() for (model) in (ids)}
    identities = defaultdict(set)
    summaries = []
    for (entry) in (specifications):
        path = resolve(base, entry["path"])
        samples, normal = header(path)
        require(
            samples == entry["sample_columns"] and normal == entry["normal_sample"],
            "Mutation source header identity differs",
        )
        model = entry["model_id"]
        identities[model].add((tuple(samples), normal))
        counts = Counter()
        reasons = Counter()
        ledger = output / (model + "." + entry["kind"] + ".csv.gz")
        with (
            gzip.open(path, "rt", encoding = "utf-8") as stream,
            gzip.open(ledger, "xt", encoding = "utf-8", newline = "") as target,
        ):
            writer = csv.writer(target, lineterminator = "\n")
            writer.writerow(["source_record", "gene", "symbol", "eligible", "reasons"])
            for (line) in (stream):
                if (line.startswith("#")):
                    continue
                counts["source_records"] += 1
                decision = classify(parse(line, samples, normal), rule, reverse, candidates)
                reasons.update(decision["reasons"])
                if (decision["eligible"]):
                    eligible[model].add(decision["gene"])
                    counts["eligible_records"] += 1
                writer.writerow(
                    [
                        counts["source_records"],
                        decision["gene"],
                        decision["symbol"],
                        int(decision["eligible"]),
                        ";".join(decision["reasons"]),
                    ]
                )
        summaries.append(
            {
                "model_id": model,
                "kind": entry["kind"],
                "counts": dict(counts),
                "reasons": dict(reasons),
            }
        )
    require(
        all(len(values) == 1 for (values) in (identities.values())),
        "SNP and indel identities differ",
    )
    normal_donors = defaultdict(set)
    donor_normals = defaultdict(set)
    for (model) in (models):
        normal = next(iter(identities[model["model_id"]]))[1]
        normal_donors[normal].add(model["donor_id"])
        donor_normals[model["donor_id"]].add(normal)
    require(
        len(normal_donors) == 58 and all(len(value) == 1 for (value) in (normal_donors.values())),
        "Normal-DNA donor linkage differs",
    )
    require(
        all(len(value) == 1 for (value) in (donor_normals.values())),
        "Within-donor normal linkage differs",
    )
    with gzip.open(
        output / "mutation_binary.csv.gz", "xt", encoding = "utf-8", newline = ""
    ) as target:
        writer = csv.writer(target, lineterminator = "\n")
        writer.writerow(["gene"] + [model["sample_id"] for (model) in (models)])
        writer.writerows(
            [mapping_values[gene]] + [int(gene in eligible[model]) for (model) in (ids)]
            for (gene) in (sorted(mapping_values))
        )
    save(
        output / "mutation_summary.json",
        {
            "sources": summaries,
            "zero_semantics": "No qualifying source call; variant-only source does not establish wild type",
            "positive_gene_counts": {key: len(value) for (key, value) in (eligible.items())},
        },
    )
    return output / "mutation_binary.csv.gz"
