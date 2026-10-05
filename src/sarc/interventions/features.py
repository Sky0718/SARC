import math
from collections import Counter
from fractions import Fraction
from .original_prediction import METHODS

def prepare_graph(baseline, changed, modules):
    if (any(
        a >= b or not 700 < score <= 1000
        for (((a, b), score)) in (list(baseline.items()) + list(changed.items()))
    )):
        raise ValueError("Canonical thresholded STRING edge input required")
    removed = baseline.keys() - changed.keys()
    added = changed.keys() - baseline.keys()
    if (any(
        baseline[edge] != changed[edge] for (edge) in (baseline.keys() & changed.keys())
    )):
        raise ValueError("B direction must preserve shared-edge confidence exactly")
    edited = sorted(removed | added)
    degree = Counter(gene for (edge) in (baseline) for (gene) in (edge))
    count = len(edited)
    if (any(gene not in modules for (edge) in (edited) for (gene) in (edge))):
        raise ValueError("Older persistent module identity is missing")
    mass = sum(1000 - min(changed[edge], 800) for (edge) in (added)) - sum(
        1000 - min(baseline[edge], 800) for (edge) in (removed)
    )
    values = {
        "log_edit_count": math.log1p(count),
        "between_module_edit_fraction": Fraction(
            sum(modules[a] != modules[b] for ((a, b)) in (edited)), count
        )
        if (count)
        else Fraction(0),
        "mean_endpoint_log_degree": sum(
            math.log1p(degree[gene]) for (edge) in (edited) for (gene) in (edge)
        )
        / (2 * count)
        if (count)
        else 0.0,
        "mean_consumed_cost_change": Fraction(mass, count) if (count) else Fraction(0),
        "added_edges": len(added),
        "removed_edges": len(removed),
    }
    return {"edited": edited, "values": values}

def sample_features(prepared, mutations, degs, method):
    if (method not in METHODS):
        raise ValueError("Unknown inherited method")
    edited = prepared["edited"]
    count = len(edited)
    values = prepared["values"] | {
        "mutated_edit_fraction": Fraction(
            sum(a in mutations or b in mutations for ((a, b)) in (edited)), count
        )
        if (count)
        else Fraction(0),
        "deg_edit_fraction": Fraction(
            sum(a in degs or b in degs for ((a, b)) in (edited)), count
        )
        if (count)
        else Fraction(0),
        "cost_structurally_not_applicable": method != "PRODIGY",
    }
    if (method != "PRODIGY"):
        values["mean_consumed_cost_change"] = None
    return values

def graph_features(baseline, changed, modules, mutations, degs, method):
    return sample_features(
        prepare_graph(baseline, changed, modules), mutations, degs, method
    )
