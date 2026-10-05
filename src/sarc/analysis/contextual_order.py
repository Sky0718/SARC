import itertools
import math
from collections import defaultdict, deque
from fractions import Fraction

def allocation(genes, scores, eligible, k = 10):
    if (genes is None):
        if (scores is not None):
            raise ValueError("Unavailable ranking exposes scores")
        return None
    if (len(genes) != len(scores) or len(set(genes)) != len(genes)):
        raise ValueError("Invalid ranking identities")
    if (any((not math.isfinite(x) for (x) in (scores))) or any(
        (a < b for ((a, b)) in (zip(scores, scores[1:])))
    )):
        raise ValueError("Invalid scores")
    rows = [(g, s) for ((g, s)) in (zip(genes, scores)) if (g in eligible)]
    weights = {g: Fraction(0) for ((g, s)) in (rows)}
    used = 0
    for (score, entries) in (itertools.groupby(rows, key = lambda x: x[1])):
        block = [g for ((g, s)) in (entries)]
        slots = min(len(block), max(0, k - used))
        weights.update({g: Fraction(slots, len(block)) for (g) in (block)})
        used += len(block)
    if (sum(weights.values()) != min(k, len(rows))):
        raise ValueError("Budget mismatch")
    return weights

def constraints(observations):
    strict = defaultdict(set)
    equal = []
    for (model, weights) in (observations.items()):
        if (weights is None):
            continue
        levels = defaultdict(list)
        for (gene, weight) in (weights.items()):
            levels[weight].append(gene)
        for (weight, genes) in (levels.items()):
            if (0 < weight < 1):
                equal.append((model, tuple(sorted(genes))))
        ordered = sorted(levels, reverse = True)
        for (index, high) in (enumerate(ordered)):
            for (low) in (ordered[index + 1 :]):
                for (a) in (levels[high]):
                    for (b) in (levels[low]):
                        strict[a, b].add(model)
    return (dict(strict), equal)

def compatibility(observations):
    strict, equal = constraints(observations)
    genes = sorted(
        {g for (w) in (observations.values()) if (w is not None) for (g) in (w)}
    )
    parent = {g: g for (g) in (genes)}

    def find(g):
        while (parent[g] != g):
            parent[g] = parent[parent[g]]
            g = parent[g]
        return g

    for (model, block) in (equal):
        for (g) in (block[1:]):
            a, b = (find(block[0]), find(g))
            parent[max(a, b)] = min(a, b)
    components = defaultdict(list)
    for (g) in (genes):
        components[find(g)].append(g)
    graph = {c: set() for (c) in (components)}
    edge_example = {}
    for (a, b) in (sorted(strict)):
        u, v = (find(a), find(b))
        graph[u].add(v)
        edge_example.setdefault(
            (u, v), {"better": a, "worse": b, "models": sorted(strict[a, b])}
        )
    incoming = {g: 0 for (g) in (graph)}
    for (u) in (graph):
        for (v) in (graph[u]):
            incoming[v] += 1
    queue = deque(sorted((g for (g) in (graph) if (incoming[g] == 0))))
    order = []
    while (queue):
        u = queue.popleft()
        order.append(u)
        for (v) in (sorted(graph[u])):
            incoming[v] -= 1
            if (incoming[v] == 0):
                queue.append(v)
    remainder = {g for (g) in (graph) if (incoming[g] > 0)}
    cycle = []
    if (remainder):
        previous = {}
        for (u) in (sorted(remainder)):
            for (v) in (sorted(graph[u] & remainder)):
                previous.setdefault(v, u)
        seen, trail = ({}, [])
        u = min(remainder)
        while (u not in seen):
            seen[u] = len(trail)
            trail.append(u)
            u = previous[u]
        backwards = trail[seen[u] :]
        cycle = list(reversed(backwards))
        cycle.append(cycle[0])
    witness = [edge_example[u, v] for ((u, v)) in (zip(cycle, cycle[1:]))]
    return (
        {
            "compatible": not bool(remainder),
            "gene_count": len(genes),
            "equality_component_count": len(components),
            "fractional_cutoff_observations": len(equal),
            "strict_gene_edges": len(strict),
            "strict_observation_constraints": sum(
                (len(v) for (v) in (strict.values()))
            ),
            "contracted_strict_self_edges": sum((a == b for ((a, b)) in (edge_example))),
            "cycle_witness": witness,
            "cycle_equality_components": {
                g: components[g] for (g) in (set(cycle)) if (len(components[g]) > 1)
            },
            "equality_blocks": [{"model": m, "genes": list(g)} for ((m, g)) in (equal)],
            "compatible_order": order if (not remainder) else None,
        },
        strict,
    )

def census(observations, donors):
    result, strict = compatibility(observations)
    model_pairs = set()
    gene_pairs = []
    for ((a, b), forward) in (sorted(strict.items())):
        if (a >= b or (b, a) not in strict):
            continue
        reverse = strict[b, a]
        for (i) in (forward):
            for (j) in (reverse):
                if (i == j):
                    raise ValueError("Self-contradictory observation")
                model_pairs.add(tuple(sorted((i, j))))
        gene_pairs.append(
            {
                "gene_a": a,
                "gene_b": b,
                "a_over_b_models": sorted(forward),
                "b_over_a_models": sorted(reverse),
            }
        )
    opportunities = []
    for (i, j) in (itertools.combinations(sorted(observations), 2)):
        if (observations[i] is not None and observations[j] is not None):
            if (len(observations[i].keys() & observations[j].keys()) >= 2):
                opportunities.append((i, j))
    lower = sum(
        (
            min(len(x["a_over_b_models"]), len(x["b_over_a_models"]))
            for (x) in (gene_pairs)
        )
    )
    result.update(
        {
            "models": len(observations),
            "available_models": sum((w is not None for (w) in (observations.values()))),
            "all_planned_observations_available": all(
                (w is not None for (w) in (observations.values()))
            ),
            "compatibility_scope": "AVAILABLE_OBSERVATIONS_ONLY_IF_ANY_MISSING",
            "unavailable_models": sorted(
                (m for ((m, w)) in (observations.items()) if (w is None))
            ),
            "empty_models": sorted(
                (m for ((m, w)) in (observations.items()) if (w == {}))
            ),
            "uninformative_all_returned_models": sorted(
                (
                    m
                    for ((m, w)) in (observations.items())
                    if (w and all((v == 1 for (v) in (w.values()))))
                )
            ),
            "model_pairs_with_two_shared_genes": len(opportunities),
            "different_donor_pairs_with_two_shared_genes": sum(
                (donors[i] != donors[j] for ((i, j)) in (opportunities))
            ),
            "direct_opposing_gene_pairs": len(gene_pairs),
            "direct_opposing_model_pairs": len(model_pairs),
            "different_donor_direct_opposing_model_pairs": sum(
                (donors[i] != donors[j] for ((i, j)) in (model_pairs))
            ),
            "models_with_direct_opposition": sorted(
                {m for (p) in (model_pairs) for (m) in (p)}
            ),
            "strict_constraint_violation_lower_bound": lower,
            "opposing_gene_pair_witnesses": gene_pairs,
            "per_model_menu_size": {
                m: len(w) if (w is not None) else None
                for ((m, w)) in (observations.items())
            },
        }
    )
    return result
