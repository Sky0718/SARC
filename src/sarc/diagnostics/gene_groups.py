from fractions import Fraction as F

def classify(lower, upper):
    return (
        "STRICTLY_NEGATIVE"
        if (upper < 0)
        else "STRICTLY_POSITIVE"
        if (lower > 0)
        else "IDENTIFIED_ZERO"
        if (lower == upper == 0)
        else "UNRESOLVED"
    )

def group_diagnostic(rows, known, lower, upper):
    ordered = sorted(rows, key = lambda row: (-max(F(), -F(row["upper"])), row["gene"]))
    negative = sorted(rows, key = lambda row: (-max(F(), -F(row["known"])), row["gene"]))
    total_negative = sum((max(F(), -F(row["known"])) for (row) in (rows)), F())
    total_positive = sum((max(F(), F(row["known"])) for (row) in (rows)), F())
    curve, removed, minimum = [], [], None
    for (count) in (range(len(ordered) + 1)):
        if (count):
            row = ordered[count - 1]
            known -= F(row["known"])
            lower -= F(row["lower"])
            upper -= F(row["upper"])
            removed.append(row["gene"])
        curve.append(
            {
                "removed_gene_count": count,
                "last_removed_gene": removed[-1] if (removed) else "",
                "known": str(known),
                "lower": str(lower),
                "upper": str(upper),
                "sign": classify(lower, upper),
            }
        )
        if (minimum is None and upper >= 0):
            minimum = {
                "gene_count": count,
                "genes": list(removed),
                "known": str(known),
                "lower": str(lower),
                "upper": str(upper),
                "sign": classify(lower, upper),
            }
    if (known != 0 or lower != 0 or upper != 0):
        raise ValueError("Complete removal must leave exact zero")
    summary = {
        "genes": len(rows),
        "minimum_complete_genes_to_lose_negative_identification": minimum,
        "fixed_prefixes": {
            str(count): {
                **curve[min(count, len(rows))],
                "genes": [row["gene"] for (row) in (ordered[:count])],
            }
            for (count) in ((1, 5, 10))
        },
        "ordering": [row["gene"] for (row) in (ordered)],
        "negative_net_known": {
            "total_negative_magnitude": str(total_negative),
            "total_positive": str(total_positive),
            "top": {
                str(count): {
                    "genes": [row["gene"] for (row) in (negative[:count])],
                    "mass": str(
                        sum((max(F(), -F(row["known"])) for (row) in (negative[:count])), F())
                    ),
                    "share": str(
                        sum((max(F(), -F(row["known"])) for (row) in (negative[:count])), F())
                        / total_negative
                    )
                    if (total_negative)
                    else None,
                }
                for (count) in ((1, 5, 10))
            },
        },
    }
    return summary, curve
