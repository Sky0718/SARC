from collections import Counter, defaultdict
from fractions import Fraction
import math
from ..evaluation.localisation import STRICT, ASSISTED, MASTERS as SEEDS

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
        return {key: plain(item) for (key, item) in (value.items())}
    if (isinstance(value, (tuple, list))):
        return [plain(item) for (item) in (value)]
    return value

def prepare_graph(baseline, changed):
    require(
        all(
            (
                a < b and 700 < score <= 1000
                for (graph) in ((baseline, changed))
                for ((a, b), score) in (graph.items())
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
            (a in candidates or b in candidates for (a, b) in (prepared["added"]))
        ),
        "removed_mutant_contact_edges": sum(
            (a in candidates or b in candidates for (a, b) in (prepared["removed"]))
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

def one_seed(prepared, candidates, ranking):
    genes, scores, status = (ranking["genes"], ranking["scores"], ranking["status"])
    require(
        len(genes) == len(scores) == len(set(genes))
        and all(
            isinstance(gene, str) and gene and gene == gene.strip() for (gene) in (genes)
        )
        and all(math.isfinite(value) for (value) in (scores))
        and not any(a < b for (a, b) in (zip(scores, scores[1:]))),
        "Malformed full native ranking",
    )
    require(
        status in ("SUCCESS", "SUCCESS_EMPTY", "FAILED", "UNAVAILABLE", "ERROR")
        and (status == "SUCCESS") == bool(genes),
        "Original ranking status differs",
    )
    eligible = [
        (gene, score) for (gene, score) in (zip(genes, scores)) if (gene in candidates)
    ]
    diagnostic = {
        "status": status,
        "native_candidate_count": len(eligible),
        "boundary_tie": None,
        "flat_top20_range": None,
    }
    if (status not in ("SUCCESS", "SUCCESS_EMPTY")):
        return {
            "status": "UNAVAILABLE",
            "reason": "OLD_BASELINE_" + status,
            "features": dict.fromkeys(ASSISTED),
            "diagnostic": diagnostic,
        }
    if (len(eligible) < 20):
        return {
            "status": "UNAVAILABLE",
            "reason": "FEWER_THAN_20_NATIVE_CANDIDATES",
            "features": dict.fromkeys(ASSISTED),
            "diagnostic": diagnostic,
        }
    top, following = (eligible[:10], eligible[10:20])
    local = burdens(prepared, [gene for (gene, _) in (top + following)])
    means = [
        sum((local[gene] for (gene, _) in (group)), Fraction()) / 10
        for (group) in ((top, following))
    ]
    exact_scores = [Fraction.from_float(float(score)) for (_, score) in (eligible[:20])]
    span = exact_scores[0] - exact_scores[19]
    diagnostic.update(
        boundary_tie = exact_scores[9] == exact_scores[10], flat_top20_range = span == 0
    )
    features = {
        "top10_touched_fraction": Fraction(
            sum((local[gene] > 0 for (gene, _) in (top))), 10
        ),
        "next10_touched_fraction": Fraction(
            sum((local[gene] > 0 for (gene, _) in (following))), 10
        ),
        "top10_mean_relative_edit_burden": means[0],
        "next10_mean_relative_edit_burden": means[1],
        "top10_minus_next10_relative_edit_burden": means[0] - means[1],
        "boundary_margin_over_top20_range": (exact_scores[9] - exact_scores[10]) / span
        if (span)
        else Fraction(),
    }
    return {
        "status": "AVAILABLE",
        "reason": None,
        "features": features,
        "diagnostic": diagnostic,
    }

def assisted_block(prepared, candidates, by_seed, method):
    expected = list(SEEDS) if (method == "PRODIGY") else [None]
    require(
        list(by_seed) == expected, "Exact complete ordered algorithm seeds required"
    )
    values = [
        {"master_seed": seed, **one_seed(prepared, candidates, by_seed[seed])}
        for (seed) in (expected)
    ]
    ready = all((item["status"] == "AVAILABLE" for (item) in (values)))
    features = (
        {
            key: sum((item["features"][key] for (item) in (values)), Fraction())
            / len(values)
            for (key) in (ASSISTED)
        }
        if (ready)
        else dict.fromkeys(ASSISTED)
    )
    return {
        "status": "AVAILABLE" if (ready) else "UNAVAILABLE",
        "reason": None if (ready) else "INCOMPLETE_REQUIRED_OLD_NATIVE_SEED_FEATURES",
        "features": features,
        "seed_details": values,
        "expected_seed_count": len(expected),
    }
