from collections import defaultdict
from fractions import Fraction

def need(condition, message):
    if (not condition):
        raise ValueError(message)

def degree_rows(rows, domain):
    need(domain and len(domain) == len(set(domain)), "Common-domain duplicates or empty")
    degrees = dict.fromkeys(domain, 0)
    observed, retained, total = set(), 0, 0
    for (row) in (rows):
        need(set(row) == {"gene1", "gene2", "combined_score"}, "Edge field mismatch")
        first, second = row["gene1"], row["gene2"]
        score = Fraction(row["combined_score"])
        need(
            score.denominator == 1 and 700 < score <= 1000,
            "Confidence is not finite integer in (700,1000]",
        )
        need(first != second, "Self edge prohibited")
        pair = tuple(sorted((first, second)))
        need(pair not in observed, "Repeated unordered pair prohibited")
        observed.add(pair)
        total += 1
        if (first in degrees and second in degrees):
            degrees[first] += 1
            degrees[second] += 1
            retained += 1
    need(sum(degrees.values()) == 2 * retained, "Undirected incidence mismatch")
    return degrees, {
        "vertices": len(domain),
        "undirected_edges": retained,
        "isolated_vertices": sum(value == 0 for (value) in (degrees.values())),
        "source_rows": total,
    }

def allocate(query, degrees, budget = 10):
    need(
        len(query) == len(set(query)) and set(query).issubset(degrees),
        "Query is duplicated or outside graph",
    )
    need(type(budget) is int and budget > 0, "Budget invalid")
    groups = defaultdict(list)
    for (gene) in (query):
        need(
            type(degrees[gene]) is int and degrees[gene] >= 0,
            "Degree must be nonnegative integer",
        )
        groups[degrees[gene]].append(gene)
    left, weights, cutoff = Fraction(budget), {}, None
    for (degree) in (sorted(groups, reverse = True)):
        if (not left):
            break
        genes = sorted(groups[degree])
        take = min(left, Fraction(len(genes)))
        weights.update((gene, take / len(genes)) for (gene) in (genes))
        left -= take
        cutoff = degree
    need(
        sum(weights.values(), Fraction()) == min(budget, len(query)),
        "Tie allocation capacity failure",
    )
    return weights, cutoff
