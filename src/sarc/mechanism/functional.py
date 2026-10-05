from collections import Counter
from fractions import Fraction

from ..precision.io import require, write_csv, write_json
from .coefficients import (
    belongs,
    categories,
    category_map,
    contrasts,
    interval,
    signed_coefficients,
    state_names,
)

def unpack_allocations(document):
    roster = document["models"]
    require(document["method"] == "DawnRank", "DawnRank component allocations required")
    require(
        document["master_seeds"] == [104729, 130363, 155921, 196613, 228017],
        "Five deterministic seed references required",
    )
    require(
        len(roster) == 85 and len({row["donor_id"] for (row) in (roster)}) == 83,
        "Full colorectal population required",
    )
    require(
        len({row["sample_id"] for (row) in (roster)}) == 85
        and len({row["model_id"] for (row) in (roster)}) == 85,
        "Nonunique model identity",
    )
    identities = {row["sample_id"]: row for (row) in (roster)}
    allocations = {}
    for (row) in (document["allocations"]):
        if (row.get("budget", 10) != 10):
            continue
        key = row["state"], row["sample_id"]
        require(key not in allocations, "Duplicate state-model allocation")
        require(row.get("numerically_certified") is True, "Uncertified numerical allocation")
        require(
            row.get("stopped_by_threshold") is True, "Numerical stopping threshold not reached"
        )
        require(row["sample_id"] in identities, "Unknown sample identity")
        require(
            all(
                row[name] == identities[row["sample_id"]][name]
                for (name) in (("model_id", "donor_id"))
            ),
            "Allocation identity mismatch",
        )
        weights = {gene: Fraction(value) for (gene, value) in (row["weights"].items())}
        query = identities[row["sample_id"]]["query"]
        require(
            len(query) == len(set(query)) and set(weights) <= set(query),
            "Allocation differs from fixed molecular query",
        )
        require(all(0 < value <= 1 for (value) in (weights.values())), "Invalid allocation weight")
        require(
            sum(weights.values(), Fraction(0)) == 10,
            "Component population requires full ten-slot returns",
        )
        allocations[key] = weights
    expected = {(state, row["sample_id"]) for (state) in (state_names()) for (row) in (roster)}
    require(set(allocations) == expected, "Incomplete eight-state population")
    return roster, allocations

def freeze(document, output):
    roster, allocations = unpack_allocations(document)
    coefficients = []
    for (name, multipliers) in (contrasts().items()):
        for (row) in (roster):
            sample = row["sample_id"]
            values = signed_coefficients(
                allocations, multipliers, [sample], {sample: Fraction(1)}
            )
            coefficients.append(
                {
                    "contrast": name,
                    "sample_id": sample,
                    "coefficients": {
                        gene: value for ((identity, gene), value) in (sorted(values.items()))
                    },
                }
            )
    requested = sorted(
        {(sample, gene) for ((state, sample), weights) in (allocations.items()) for (gene) in (weights)}
    )
    write_json(output / "allocations.json", document)
    write_json(output / "coefficients.json", {"contrasts": contrasts(), "rows": coefficients})
    write_json(
        output / "label_requests.json",
        {"requests": [{"sample_id": sample, "gene": gene} for (sample, gene) in (requested)]},
    )
    return {
        "models": 85,
        "donors": 83,
        "states": 680,
        "contrasts": 22,
        "requested_labels": len(requested),
    }

def selected_labels(document, requested):
    labels = {}
    allowed = ("MEASURED", "ABSENT", "EXPLICIT_NULL", "IDENTITY_CONFLICT", "CONFLICTING_CALLS")
    for (row) in (document["calls"]):
        key = row["sample_id"], row["gene"]
        require(key not in labels, "Duplicate sample-gene label")
        state, value = row["state"], row["value"]
        require(state in allowed, "Unrecognised assay state")
        if (state == "MEASURED"):
            require(type(value) is int and value in (0, 1), "Measured binary value required")
        else:
            require(value is None, "Unknown assay state must have null value")
        labels[key] = state, value
    require(
        set(labels) == set(requested), "Selected labels do not equal frozen allocation support"
    )
    return labels

def biological_groups(roster):
    samples = [row["sample_id"] for (row) in (roster)]
    donors = {row["sample_id"]: row["donor_id"] for (row) in (roster)}
    counts = Counter(donors.values())
    groups = [("model", sample, [sample], {sample: Fraction(1)}) for (sample) in (samples)]
    for (donor) in (sorted(counts)):
        chosen = [sample for (sample) in (samples) if (donors[sample] == donor)]
        groups.append(
            ("donor", donor, chosen, {sample: Fraction(1, len(chosen)) for (sample) in (chosen)})
        )
    groups.append(
        (
            "population",
            "ALL_83_DONORS",
            samples,
            {sample: Fraction(1, len(counts) * counts[donors[sample]]) for (sample) in (samples)},
        )
    )
    return groups

def state_mass(allocations, state, selected, samples, alpha, labels, mapping):
    fields = ("R", "H", "Z", "ABSENT", "NULL", "IDENTITY_CONFLICT", "CONFLICTING_CALLS", "U")
    masses = dict.fromkeys(fields, Fraction(0))
    for (sample) in (samples):
        for (gene, weight) in (allocations[(state, sample)].items()):
            if (not belongs(gene, selected, mapping)):
                continue
            amount = alpha[sample] * weight
            masses["R"] += amount
            label_state, value = labels[(sample, gene)]
            if (label_state == "MEASURED"):
                masses["H" if (value) else "Z"] += amount
            else:
                masses["NULL" if (label_state == "EXPLICIT_NULL") else label_state] += amount
                masses["U"] += amount
    require(
        masses["R"] == masses["H"] + masses["Z"] + masses["U"], "Allocation conservation failed"
    )
    masses["VACANT"] = Fraction(10) - masses["R"] if (selected == "ALL") else None
    return masses

def evaluate(document, label_document, category_document, output):
    roster, allocations = unpack_allocations(document)
    requested = {
        (sample, gene) for ((state, sample), weights) in (allocations.items()) for (gene) in (weights)
    }
    labels = selected_labels(label_document, requested)
    genes = {gene for (sample, gene) in (requested)}
    mapping = category_map(genes, category_document)
    state_rows, contrast_rows = [], []
    for (level, unit, samples, alpha) in (biological_groups(roster)):
        for (state) in (state_names()):
            for (selected) in (categories()):
                masses = state_mass(
                    allocations, state, selected, samples, alpha, labels, mapping
                )
                state_rows.append(
                    {
                        "level": level,
                        "unit": unit,
                        "state": state,
                        "category": selected,
                        **masses,
                    }
                )
        for (name, multipliers) in (contrasts().items()):
            combined = signed_coefficients(allocations, multipliers, samples, alpha)
            for (selected) in (categories()):
                values = {
                    key: value
                    for (key, value) in (combined.items())
                    if (belongs(key[1], selected, mapping))
                }
                contrast_rows.append(
                    {
                        "level": level,
                        "unit": unit,
                        "contrast": name,
                        "category": selected,
                        **interval(values, labels),
                    }
                )
    directions = []
    for (name) in (contrasts()):
        for (selected) in (categories()):
            counts = Counter(
                row["sign"]
                for (row) in (contrast_rows)
                if (
                    row["level"] == "donor"
                    and row["contrast"] == name
                    and row["category"] == selected
                )
            )
            directions.append(
                {
                    "contrast": name,
                    "category": selected,
                    **{
                        key: counts[key]
                        for (key) in (("NEGATIVE", "POSITIVE", "ZERO", "UNRESOLVED"))
                    },
                }
            )
    check_partitions(state_rows, contrast_rows)
    write_csv(output / "state_levels.csv", state_rows)
    write_csv(output / "contrasts.csv", contrast_rows)
    result = {
        "models": 85,
        "donors": 83,
        "states": 680,
        "contrasts": 22,
        "population_levels": [row for (row) in (state_rows) if (row["level"] == "population")],
        "population_contrasts": [
            row for (row) in (contrast_rows) if (row["level"] == "population")
        ],
        "donor_directions": directions,
        "uncertainty": "Finite shared-binary-label bounds, not sampling intervals",
    }
    write_json(output / "results.json", result)
    return result

def check_partitions(state_rows, contrast_rows):
    for (rows, key_name, fields) in ((
        (
            state_rows,
            "state",
            ("R", "H", "Z", "ABSENT", "NULL", "IDENTITY_CONFLICT", "CONFLICTING_CALLS", "U"),
        ),
        (contrast_rows, "contrast", ("R", "H", "lower", "upper")),
    )):
        lookup = {
            (row["level"], row["unit"], row[key_name], row["category"]): row for (row) in (rows)
        }
        for (row) in (rows):
            if (row["category"] != "ALL"):
                continue
            key = row["level"], row["unit"], row[key_name]
            for (field) in (fields):
                total = sum(
                    (
                        lookup[(*key, category)][field]
                        for (category) in ((
                            "BOTH_CORE",
                            "ORGANOID_ONLY",
                            "CELL_ONLY",
                            "NEITHER_CORE",
                        ))
                    ),
                    Fraction(0),
                )
                require(row[field] == total, "Four-way category partition failed")
                require(
                    row[field]
                    == lookup[(*key, "BOTH_CORE")][field]
                    + lookup[(*key, "NOT_BOTH_CORE")][field],
                    "Core complement partition failed",
                )
