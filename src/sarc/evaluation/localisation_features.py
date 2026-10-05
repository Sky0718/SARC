from collections import Counter, defaultdict
from fractions import Fraction

STRICT = (
    "added_mutant_contact_edges",
    "removed_mutant_contact_edges",
    "mutant_candidates_touched_fraction",
    "mean_mutant_relative_edit_burden",
    "max_mutant_relative_edit_burden",
    "mean_mutant_deg_onehop_edit_burden",
)

def require(value, message):
    if (not value):
        raise ValueError(message)

def plain(value):
    if (isinstance(value, Fraction)):
        return {
            "numerator": value.numerator,
            "denominator": value.denominator,
            "value": float(value),
        }
    if (isinstance(value, dict)):
        return {key: plain(item) for ((key, item)) in (value.items())}
    if (isinstance(value, (tuple, list))):
        return [plain(item) for (item) in (value)]
    return value

def prepare_graph(baseline, changed):
    require(
        all(
            (
                a < b and 700 < score <= 1000
                for (graph) in ((baseline, changed))
                for (((a, b), score)) in (graph.items())
            )
        ),
        "Canonical score>700 STRING graph required",
    )
    require(
        all(
            (
                baseline[edge] == changed[edge]
                for (edge) in (baseline.keys() & changed.keys())
            )
        ),
        "Shared edge score changed",
    )
    added, removed = (
        changed.keys() - baseline.keys(),
        baseline.keys() - changed.keys(),
    )
    incident = defaultdict(set)
    for (a, b) in (added | removed):
        incident[a].add(b)
        incident[b].add(a)
    degrees = Counter((gene for (edge) in (baseline) for (gene) in (edge)))
    return {
        "added": added,
        "removed": removed,
        "incident": incident,
        "degrees": degrees,
    }

def burdens(prepared, candidates):
    return {
        gene: Fraction(
            len(prepared["incident"].get(gene, ())), max(1, prepared["degrees"][gene])
        )
        for (gene) in (candidates)
    }

def strict_block(prepared, candidates, degs):
    if (not candidates):
        return {
            "status": "UNAVAILABLE",
            "reason": "EMPTY_NATIVE_MUTANT_POPULATION",
            "features": dict.fromkeys(STRICT),
            "candidate_count": 0,
        }
    values = burdens(prepared, candidates)
    features = {
        "added_mutant_contact_edges": sum(
            (a in candidates or b in candidates for ((a, b)) in (prepared["added"]))
        ),
        "removed_mutant_contact_edges": sum(
            (a in candidates or b in candidates for ((a, b)) in (prepared["removed"]))
        ),
        "mutant_candidates_touched_fraction": Fraction(
            sum((bool(prepared["incident"].get(gene)) for (gene) in (candidates))),
            len(candidates),
        ),
        "mean_mutant_relative_edit_burden": sum(values.values(), Fraction())
        / len(candidates),
        "max_mutant_relative_edit_burden": max(values.values()),
        "mean_mutant_deg_onehop_edit_burden": sum(
            (
                Fraction(
                    len(prepared["incident"].get(gene, set()) & degs),
                    max(1, prepared["degrees"][gene]),
                )
                for (gene) in (candidates)
            ),
            Fraction(),
        )
        / len(candidates),
    }
    return {
        "status": "AVAILABLE",
        "reason": None,
        "features": features,
        "candidate_count": len(candidates),
    }
