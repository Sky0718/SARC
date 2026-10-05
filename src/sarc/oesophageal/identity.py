import re
from collections import defaultdict
from .contracts import require
from .mutation import mapping

def stable(value):
    return re.sub(r"\.\d+$", "", value)

def crosswalk(tumour, normal, metadata, original_genes):
    symbols, candidates = mapping(metadata)
    records = []
    for (gene) in (sorted(set(tumour) | set(normal))):
        reasons = []
        if (not gene or len(tumour.get(gene, [])) != 1):
            reasons.append("ABSENT_OR_AMBIGUOUS_TUMOUR")
        if (len(normal.get(gene, [])) != 1):
            reasons.append("ABSENT_OR_AMBIGUOUS_NORMAL")
        if (gene not in symbols):
            reasons.append("NO_UNIQUE_FROZEN_HGNC_MAPPING")
        records.append(
            {
                "stable_ensembl": gene,
                "hgnc_symbol": symbols.get(gene),
                "status": ";".join(reasons) if (reasons) else "ONE_TO_ONE_SHARED",
                "cmp_identities": sorted(tumour.get(gene, [])),
                "gtex_identities": normal.get(gene, []),
                "in_original_7399": symbols.get(gene) in original_genes,
            }
        )
    chosen = [row for (row) in (records) if (row["status"] == "ONE_TO_ONE_SHARED")]
    require(len(chosen) == 39102, "Shared count-gene support differs")
    require(
        len({row["hgnc_symbol"] for (row) in (chosen)} & set(original_genes)) == 7388,
        "Common query support differs",
    )
    return records, chosen
