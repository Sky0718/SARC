from fractions import Fraction

FEATURES = (
    "old_eligible_coverage",
    "new_eligible_coverage",
    "persistent_eligible_coverage",
    "incident_topology_turnover",
    "common_incident_consumed_cost_change",
)

def unique_strings(values):
    values = tuple(values)
    if ((len(values) != len(set(values))) or any(
        (not isinstance(value, str)) or (not value) for (value) in (values)
    )):
        raise ValueError("Identifiers must be nonempty unique strings")
    return values

def network_index(network):
    if (set(network) != {"release", "names", "edges"}):
        raise ValueError("Only input network identities and edges are permitted")
    names = dict(network["names"])
    unique_strings(names)
    unique_strings(names.values())
    edges = {}
    incident = {}
    for (left, right, confidence) in (network["edges"]):
        if (
            (left == right)
            or (left not in names)
            or (right not in names)
            or (type(confidence) is not int)
            or (not (700 < confidence <= 1000))
        ):
            raise ValueError("Invalid frozen physical edge")
        key = tuple(sorted((left, right)))
        if (key in edges):
            raise ValueError("Duplicate undirected edge")
        edges[key] = confidence
        for (node) in (key):
            incident.setdefault(names[node], set()).add(key)
    return {
        "release": network["release"],
        "names": names,
        "edges": edges,
        "incident": incident,
        "nodes": frozenset(incident),
    }

def stable_edge(key, network, stable):
    return tuple(
        sorted(
            (
                ("persistent", node, network["names"][node])
                if (node in stable)
                else (network["release"], node, network["names"][node])
            )
            for (node) in (key)
        )
    )

def cohort_features(samples, old_network, new_network, persistent):
    samples = list(samples)
    if ((not samples) or any(
        set(sample) != {"sample_id", "eligible_ids"} for (sample) in (samples)
    )):
        raise ValueError("Input-only sample projection is required")
    unique_strings(sample["sample_id"] for (sample) in (samples))
    old = network_index(old_network)
    new = network_index(new_network)
    if (old["release"] == new["release"]):
        raise ValueError("A release transition is required")
    stable = {
        node
        for (node) in (set(old["names"]) & set(new["names"]))
        if (old["names"][node] == new["names"][node])
    }
    if (dict(persistent) != {node: old["names"][node] for (node) in (stable)}):
        raise ValueError(
            "Persistent source identity is not the exact unchanged unique mapping"
        )
    stable_names = frozenset(persistent.values())
    records = []
    for (sample) in (samples):
        eligible = frozenset(unique_strings(sample["eligible_ids"]))
        old_keys = set().union(
            *(old["incident"].get(gene, set()) for (gene) in (eligible))
        )
        new_keys = set().union(
            *(new["incident"].get(gene, set()) for (gene) in (eligible))
        )
        old_edges = {
            stable_edge(key, old, stable): old["edges"][key] for (key) in (old_keys)
        }
        new_edges = {
            stable_edge(key, new, stable): new["edges"][key] for (key) in (new_keys)
        }
        union = set(old_edges) | set(new_edges)
        common = set(old_edges) & set(new_edges)
        values = [
            Fraction(len(eligible & nodes), len(eligible)) if (eligible) else None
            for (nodes) in ((old["nodes"], new["nodes"], stable_names))
        ]
        values.append(
            Fraction(len(set(old_edges) ^ set(new_edges)), len(union))
            if (union)
            else None
        )
        values.append(
            sum(
                (
                    Fraction(
                        abs(min(old_edges[key], 800) - min(new_edges[key], 800)), 100
                    )
                    for (key) in (common)
                ),
                Fraction(0),
            )
            / len(common)
            if (common)
            else None
        )
        records.append(
            {
                "sample_id": sample["sample_id"],
                "eligible_count": len(eligible),
                "incident_union_count": len(union),
                "common_incident_count": len(common),
                "values": dict(zip(FEATURES, values)),
            }
        )
    values = {}
    missing = {}
    for (feature) in (FEATURES):
        current = [row["values"][feature] for (row) in (records)]
        missing[feature] = sum(value is None for (value) in (current))
        values[feature] = (
            None if (missing[feature]) else sum(current, Fraction(0)) / len(current)
        )
    return {
        "values": values,
        "missing_sample_counts": missing,
        "sample_count": len(samples),
        "samples": records,
    }

def feature_values(record):
    values = record["values"] if ("values" in record) else record
    if (set(values) != set(FEATURES)):
        raise ValueError("The complete five-feature input projection is required")
    if (any(
        (value is not None)
        and (
            (not isinstance(value, (int, Fraction)))
            or (isinstance(value, bool))
            or (not (0 <= value <= 1))
        )
        for (value) in (values.values())
    )):
        raise ValueError("Input features must be exact values in [0,1] or missing")
    return values

def distance(target, source):
    target = feature_values(target)
    source = feature_values(source)
    differences = []
    for (name) in (FEATURES):
        left, right = target[name], source[name]
        differences.append(
            Fraction(0)
            if ((left is None) and (right is None))
            else Fraction(1)
            if ((left is None) or (right is None))
            else abs(left - right)
        )
    return (
        (sum(differences[:2], Fraction(0)) / 2)
        + (sum(differences[2:], Fraction(0)) / 3)
    ) / 2

def context_weights(target_features, training_features):
    if (not training_features):
        raise ValueError("At least one development training context is required")
    raw = {
        context: (
            Fraction(1) / (1 + distance(target_features, training_features[context]))
        )
        for (context) in (sorted(training_features))
    }
    total = sum(raw.values(), Fraction(0))
    return {context: (value / total) for ((context, value)) in (raw.items())}
