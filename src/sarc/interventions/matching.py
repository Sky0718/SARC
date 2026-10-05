from .balance import degree_bin
from collections import Counter, defaultdict

CLASSES = ("PB", "PW", "DB", "DW")
CONTROL_SEEDS = (271828, 314159, 161803, 141421, 173205)

def confidence_bin(score):
    return 0 if (score < 750) else 1 if (score < 800) else 2 if (score < 900) else 3

def degree_pair(edge, degrees):
    return tuple(sorted((degree_bin(degrees[edge[0]]), degree_bin(degrees[edge[1]]))))

def removal_stratum(edge, score, degrees):
    return (*degree_pair(edge, degrees), 1000 - min(score, 800), confidence_bin(score))

def jointly_match_removals(baseline, targets, degrees, rng):
    available = defaultdict(list)
    requested = defaultdict(lambda: defaultdict(int))
    for (edge, score) in (baseline.items()):
        available[removal_stratum(edge, score, degrees)].append(edge)
    for (label) in (CLASSES):
        for (edge) in (targets[label]):
            requested[removal_stratum(edge, baseline[edge], degrees)][label] += 1
    selected = {label: set() for (label) in (CLASSES)}
    minimum_slack = None
    for (stratum) in (sorted(requested)):
        needs = requested[stratum]
        count = sum(needs.values())
        pool = sorted(available[stratum])
        if (len(pool) < count):
            raise ValueError("Infeasible exact removal stratum: " + repr(stratum))
        slack = len(pool) - count
        minimum_slack = slack if (minimum_slack is None) else min(minimum_slack, slack)
        chosen = rng.sample(pool, count)
        offset = 0
        for (label) in (CLASSES):
            selected[label].update(chosen[offset : offset + needs[label]])
            offset += needs[label]
    return (
        selected,
        {"strata": len(requested), "minimum_available_minus_required": minimum_slack},
    )

def jointly_match_additions(baseline, observed, targets, degrees, rng):
    nodes = defaultdict(list)
    for (node) in (sorted(degrees)):
        nodes[degree_bin(degrees[node])].append(node)
    excluded_counts = Counter((degree_pair(edge, degrees) for (edge) in (baseline)))
    requested = defaultdict(lambda: defaultdict(list))
    for (label) in (CLASSES):
        for (edge) in (sorted(targets[label])):
            requested[degree_pair(edge, degrees)][label].append(observed[edge])
    selected = {label: {} for (label) in (CLASSES)}
    feasibility = []
    for (stratum) in (sorted(requested)):
        (first, second) = (nodes[index] for (index) in (stratum))
        total_pairs = (
            len(first) * len(second)
            if (stratum[0] != stratum[1])
            else len(first) * (len(first) - 1) // 2
        )
        capacity = total_pairs - excluded_counts[stratum]
        count = sum((len(values) for (values) in (requested[stratum].values())))
        if (capacity < count):
            raise ValueError("Infeasible exact addition stratum: " + repr(stratum))
        chosen = set()
        attempts = 0
        while (len(chosen) < count):
            (a, b) = (rng.choice(first), rng.choice(second))
            attempts += 1
            if (attempts > max(10000, count * 1000)):
                raise RuntimeError(
                    "Addition sampling budget exhausted; preserve failure without relaxing matching"
                )
            if (a == b):
                continue
            edge = tuple(sorted((a, b)))
            if (edge not in baseline):
                chosen.add(edge)
        chosen = sorted(chosen)
        rng.shuffle(chosen)
        offset = 0
        for (label) in (CLASSES):
            scores = list(requested[stratum][label])
            rng.shuffle(scores)
            for (edge, score) in (zip(chosen[offset : offset + len(scores)], scores)):
                selected[label][edge] = score
            offset += len(scores)
        feasibility.append(
            {
                "degree_pair": stratum,
                "required": count,
                "finite_nonedge_capacity": capacity,
                "sampling_attempts": attempts,
            }
        )
    return (selected, feasibility)

def verify_matches(
    baseline,
    observed,
    targets_remove,
    targets_add,
    controls_remove,
    controls_add,
    degrees,
):
    removal_union = set()
    addition_union = set()
    result = {}
    for (label) in (CLASSES):
        removed = controls_remove[label]
        added = controls_add[label]
        if (removal_union & removed or addition_union & added.keys()):
            raise ValueError("Control classes overlap")
        removal_union.update(removed)
        addition_union.update(added)
        target_remove = Counter(
            (
                removal_stratum(edge, baseline[edge], degrees)
                for (edge) in (targets_remove[label])
            )
        )
        control_remove = Counter(
            (removal_stratum(edge, baseline[edge], degrees) for (edge) in (removed))
        )
        target_add = Counter(
            (
                (*degree_pair(edge, degrees), observed[edge])
                for (edge) in (targets_add[label])
            )
        )
        control_add = Counter(
            ((*degree_pair(edge, degrees), score) for ((edge, score)) in (added.items()))
        )
        if (target_remove != control_remove or target_add != control_add):
            raise ValueError("Incorrect exact-stratum match")
        if (not removed <= baseline.keys() or added.keys() & baseline.keys()):
            raise ValueError("Invalid control operation")

        def averages(values, scores):
            edges = list(values)
            return {
                "mean_endpoint_degree": sum(
                    (degrees[node] for (edge) in (edges) for (node) in (edge))
                )
                / (2 * len(edges))
                if (edges)
                else None,
                "mean_confidence": sum((scores[edge] for (edge) in (edges)))
                / len(edges)
                if (edges)
                else None,
            }

        result[label] = {
            "removed": len(removed),
            "added": len(added),
            "exact_degree_cost_confidence_strata": True,
            "original_removed_edges_retained": len(removed & targets_remove[label]),
            "original_added_edges_retained": len(added.keys() & targets_add[label]),
            "real_removal": averages(targets_remove[label], baseline),
            "control_removal": averages(removed, baseline),
            "real_addition": averages(targets_add[label], observed),
            "control_addition": averages(added, added),
        }
    return result
