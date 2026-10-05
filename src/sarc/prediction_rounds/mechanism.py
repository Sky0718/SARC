from dataclasses import dataclass
import io
from types import MappingProxyType
import numpy as np

FEATURE_NAMES = (
    "mean_normalised_absolute_direct_delta",
    "max_normalised_absolute_direct_delta",
    "mean_normalised_absolute_reweighting_delta",
    "max_normalised_absolute_reweighting_delta",
    "mean_normalised_absolute_interaction_delta",
    "max_normalised_absolute_interaction_delta",
    "mean_normalised_absolute_total_delta",
    "max_normalised_absolute_total_delta",
    "similarity_change_fraction",
    "mean_absolute_similarity_change",
    "candidate_entry_fraction",
    "candidate_exit_fraction",
)
RTOL = 1e-10
ATOL = 1e-12

def require(value, message):
    if (not value):
        raise ValueError(message)

def parse_pathways(text):
    require(isinstance(text, str), "Bound historical pathway text required")
    pathways = {}
    for (line) in (io.StringIO(text, newline = None)):
        fields = line.rstrip("\n").split("\t\t\t")
        require(len(fields) >= 2, "Malformed original triple-tab pathway row")
        for (gene) in (set(fields[1].split(",")[:-1])):
            pathways.setdefault(gene, set()).add(fields[0])
    return {gene: frozenset(names) for (gene, names) in (pathways.items())}

def canonical_edges(edges):
    require(hasattr(edges, "items"), "R2 canonical edge-to-confidence mapping required")
    result = {}
    for (pair, confidence) in (edges.items()):
        require(
            isinstance(pair, tuple)
            and len(pair) == 2
            and all((isinstance(gene, str) and gene for (gene) in (pair)))
            and (pair[0] < pair[1]),
            "Canonical nonself undirected edge required",
        )
        require(
            isinstance(confidence, (int, np.integer))
            and (not isinstance(confidence, (bool, np.bool_)))
            and (700 < confidence <= 1000),
            "Original canonical STRING confidence contract differs",
        )
        result[pair] = int(confidence)
    return result

def immutable(array):
    array.setflags(write = False)
    return array

@dataclass(frozen = True)
class State:
    q: np.ndarray
    contacts: np.ndarray
    own_contacts: np.ndarray
    retained_degs: np.ndarray
    similarity: np.ndarray
    support: np.ndarray
    scores: np.ndarray

@dataclass(frozen = True)
class Baseline:
    population: tuple
    mutation_genes: tuple
    deg_genes: tuple
    mutation_index: object
    deg_index: object
    mutation_mask: np.ndarray
    deg_mask: np.ndarray
    pathways: object
    edges: object
    state: State

def accumulate(cache, edges, sign, q, contacts, own_contacts):
    for (left, right) in (edges):
        for (mutation, deg) in (((left, right), (right, left))):
            u, v = (cache.mutation_index.get(mutation), cache.deg_index.get(deg))
            if (u is None or v is None):
                continue
            eligible = cache.mutation_mask[:, u] & cache.deg_mask[:, v]
            if (not np.any(eligible)):
                continue
            weight = len(
                cache.pathways.get(mutation, frozenset())
                & cache.pathways.get(deg, frozenset())
            )
            contacts[u, eligible] += sign
            q[u, eligible] += sign * weight
            own_contacts[eligible, v] += sign

def state_from_counts(cache, q, contacts, own_contacts):
    require(
        np.all(q >= 0) and np.all(contacts >= 0) and np.all(own_contacts >= 0),
        "Input-derived contact subtraction became negative",
    )
    retained = own_contacts > 0
    integer = retained.astype(np.int64)
    intersections = integer @ integer.T
    sizes = integer.sum(axis = 1)
    denominator = np.outer(sizes, sizes)
    similarity = np.divide(
        intersections.astype(np.float64) ** 2,
        denominator,
        out = np.zeros(denominator.shape, dtype = np.float64),
        where = denominator != 0,
    )
    support = cache.mutation_mask & (contacts.sum(axis = 1) > 0)[None, :]
    scores = similarity @ q.T * cache.mutation_mask
    require(
        np.isfinite(scores).all() and np.isfinite(similarity).all(),
        "Nonfinite input-derived arithmetic",
    )
    return State(
        *(
            immutable(value)
            for (value) in ((
                q,
                contacts,
                own_contacts,
                retained,
                similarity,
                support,
                scores,
            ))
        )
    )

def prepare_baseline(baseline_edges, population, mutations, degs, pathways):
    population = tuple(population)
    require(
        population
        and len(population) == len(set(population))
        and all(
            (
                isinstance(sample, str)
                and sample
                and ("_" not in sample)
                and ("." not in sample)
                and (not any((char.isspace() for (char) in (sample))))
                for (sample) in (population)
            )
        ),
        "Unique canonical samples without ambiguous native underscore tokens required",
    )
    require(
        set(mutations) == set(degs) == set(population),
        "Complete same-cohort mutation and DEG populations required",
    )
    mutation_sets = {sample: frozenset(mutations[sample]) for (sample) in (population)}
    deg_sets = {sample: frozenset(degs[sample]) for (sample) in (population)}
    mutation_genes = tuple(sorted(set().union(*mutation_sets.values())))
    deg_genes = tuple(sorted(set().union(*deg_sets.values())))
    require(
        all(
            (
                isinstance(gene, str)
                and gene
                and (not any((char.isspace() for (char) in (gene))))
                for (gene) in (mutation_genes + deg_genes)
            )
        ),
        "Valid existing literal molecular labels required",
    )
    require(
        all(("_" not in gene and "." not in gene for (gene) in (deg_genes))),
        "DEG labels must already follow admitted native normalisation without ambiguous underscore tokens",
    )
    require(
        hasattr(pathways, "items")
        and all(
            (
                isinstance(gene, str)
                and isinstance(names, (set, frozenset))
                and all((isinstance(name, str) for (name) in (names)))
                for (gene, names) in (pathways.items())
            )
        ),
        "Explicit historical pathway membership sets required",
    )
    edges = canonical_edges(baseline_edges)
    mutation_index = {gene: index for (index, gene) in (enumerate(mutation_genes))}
    deg_index = {gene: index for (index, gene) in (enumerate(deg_genes))}
    mutation_mask = np.asarray(
        [
            [gene in mutation_sets[sample] for (gene) in (mutation_genes)]
            for (sample) in (population)
        ],
        dtype = bool,
    ).reshape(len(population), len(mutation_genes))
    deg_mask = np.asarray(
        [[gene in deg_sets[sample] for (gene) in (deg_genes)] for (sample) in (population)],
        dtype = bool,
    ).reshape(len(population), len(deg_genes))
    cache = Baseline(
        population,
        mutation_genes,
        deg_genes,
        MappingProxyType(mutation_index),
        MappingProxyType(deg_index),
        immutable(mutation_mask),
        immutable(deg_mask),
        MappingProxyType({gene: frozenset(names) for (gene, names) in (pathways.items())}),
        MappingProxyType(edges),
        None,
    )
    q = np.zeros((len(mutation_genes), len(population)), dtype = np.int64)
    contacts = np.zeros_like(q)
    own_contacts = np.zeros((len(population), len(deg_genes)), dtype = np.int64)
    accumulate(cache, edges, 1, q, contacts, own_contacts)
    state = state_from_counts(cache, q, contacts, own_contacts)
    return Baseline(
        cache.population,
        cache.mutation_genes,
        cache.deg_genes,
        cache.mutation_index,
        cache.deg_index,
        cache.mutation_mask,
        cache.deg_mask,
        cache.pathways,
        cache.edges,
        state,
    )

def decompose(cache, changed_edges):
    require(
        isinstance(cache, Baseline) and isinstance(cache.state, State),
        "Prepared immutable whole-cohort baseline required",
    )
    changed = canonical_edges(changed_edges)
    require(
        all(
            (
                cache.edges[edge] == changed[edge]
                for (edge) in (cache.edges.keys() & changed.keys())
            )
        ),
        "Shared-edge confidence changed outside the admitted R2 graph contract",
    )
    q, contacts, own_contacts = (
        cache.state.q.copy(),
        cache.state.contacts.copy(),
        cache.state.own_contacts.copy(),
    )
    accumulate(
        cache, cache.edges.keys() - changed.keys(), -1, q, contacts, own_contacts
    )
    accumulate(cache, changed.keys() - cache.edges.keys(), 1, q, contacts, own_contacts)
    after = state_from_counts(cache, q, contacts, own_contacts)
    delta_q = after.q - cache.state.q
    delta_p = after.similarity - cache.state.similarity
    direct = cache.state.similarity @ delta_q.T * cache.mutation_mask
    reweighting = delta_p @ cache.state.q.T * cache.mutation_mask
    interaction = delta_p @ delta_q.T * cache.mutation_mask
    total = direct + reweighting + interaction
    require(
        np.allclose(total, after.scores - cache.state.scores, rtol = RTOL, atol = ATOL),
        "Input-derived S/P/Q algebra differs beyond frozen tolerance",
    )
    return {
        "before": cache.state,
        "after": after,
        "direct": direct,
        "reweighting": reweighting,
        "interaction": interaction,
        "total": total,
        "similarity_delta": delta_p,
    }

def mechanism_features(cache, changed_edges):
    parts = decompose(cache, changed_edges)
    before, after = (parts["before"], parts["after"])
    result = {}
    for (index, sample) in (enumerate(cache.population)):
        old, new = (before.support[index], after.support[index])
        union = old | new
        count = int(np.count_nonzero(union))
        scale = max(1.0, float(before.scores[index, old].sum()))
        values = []
        for (component) in (("direct", "reweighting", "interaction", "total")):
            absolute = np.abs(parts[component][index, union]) / scale
            values.extend(
                (
                    float(absolute.mean()) if (count) else 0.0,
                    float(absolute.max()) if (count) else 0.0,
                )
            )
        similarity = parts["similarity_delta"][index]
        values.extend(
            (
                float(np.count_nonzero(similarity) / len(cache.population)),
                float(np.abs(similarity).mean()),
                float(np.count_nonzero(new & ~old) / max(1, count)),
                float(np.count_nonzero(old & ~new) / max(1, count)),
            )
        )
        require(
            len(values) == len(FEATURE_NAMES)
            and np.isfinite(values).all()
            and all((value >= 0 for (value) in (values))),
            "Exactly twelve finite nonnegative input features required",
        )
        result[sample] = {
            "status": "AVAILABLE",
            "reason": None,
            "features": dict(zip(FEATURE_NAMES, values)),
            "diagnostics": {
                "old_candidate_count": int(np.count_nonzero(old)),
                "new_candidate_count": int(np.count_nonzero(new)),
                "union_candidate_count": count,
                "score_scale": scale,
                "old_retained_deg_count": int(
                    np.count_nonzero(before.retained_degs[index])
                ),
                "new_retained_deg_count": int(
                    np.count_nonzero(after.retained_degs[index])
                ),
            },
        }
    return result
