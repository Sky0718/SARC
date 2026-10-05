from fractions import Fraction

from .endpoint import require
from .inputs import CATEGORIES, validate_models
from .paired import combine

LEVEL_FIELDS = (
    "returned",
    "positive",
    "zero",
    "unknown",
    "absent",
    "null",
    "identity_conflict",
    "call_conflict",
)
DELTA_FIELDS = ("returned_change", "K", "L", "U", "zero_change", "unknown_change")
UNKNOWN_FIELDS = {
    "ABSENT": "absent",
    "EXPLICIT_NULL": "null",
    "IDENTITY_CONFLICT": "identity_conflict",
    "CONFLICTING_CALLS": "call_conflict",
}

def direction(lower, upper):
    if (upper < 0):
        return "NEGATIVE"
    if (lower > 0):
        return "POSITIVE"
    return "ZERO" if (lower == upper == 0) else "UNRESOLVED"

def contribution(coefficient, call):
    state, label = call
    positive = coefficient if (state == "MEASURED" and label == 1) else Fraction()
    zero = coefficient if (state == "MEASURED" and label == 0) else Fraction()
    unknown = coefficient if (label is None) else Fraction()
    return {
        "returned_change": coefficient,
        "K": positive,
        "L": positive + min(unknown, 0),
        "U": positive + max(unknown, 0),
        "zero_change": zero,
        "unknown_change": unknown,
    }

def finish_delta(row):
    row["direction"] = direction(row["L"], row["U"])
    return row

def summarise(vector):
    values = list(vector.values())
    positive = sum((value for (value) in (values) if (value > 0)), Fraction())
    negative = -sum((value for (value) in (values) if (value < 0)), Fraction())
    output = {
        "units": len(values),
        "positive_count": sum(value > 0 for (value) in (values)),
        "negative_count": sum(value < 0 for (value) in (values)),
        "zero_count": sum(value == 0 for (value) in (values)),
        "net": positive - negative,
        "positive_mass": positive,
        "negative_mass": negative,
        "absolute_mass": positive + negative,
        "minimum": min(values),
        "maximum": max(values),
    }
    for (side, parts, denominator) in ([
        ("positive", [max(value, 0) for (value) in (values)], positive),
        ("negative", [max(-value, 0) for (value) in (values)], negative),
        ("absolute", [abs(value) for (value) in (values)], positive + negative),
    ]):
        ordered = sorted(parts, reverse = True)
        for (number) in ((1, 3, 5, 10)):
            output[f"top_{number}_{side}_share"] = (
                sum(ordered[:number], Fraction()) / denominator if (denominator) else None
            )
    return output

def aggregate_rows(rows, fields, counts, unit):
    buckets = {}
    for (row) in (rows):
        context = tuple(
            (key, row[key])
            for (key) in (("context", "policy", "partition", "category", "state", "contrast"))
            if (key in row)
        )
        key = (row["donor_id"], context) if (unit == "donor") else context
        if (key not in buckets):
            base = {"donor_id": row["donor_id"]} if (unit == "donor") else {}
            buckets[key] = {**base, **dict(context), **{field: Fraction() for (field) in (fields)}}
        weight = Fraction(1, counts[row["donor_id"]])
        if (unit == "population"):
            weight /= len(counts)
        for (field) in (fields):
            buckets[key][field] += weight * row[field]
    result = list(buckets.values())
    for (row) in (result):
        if ("L" in fields):
            finish_delta(row)
        else:
            row["pooled_positive_fraction"] = (
                row["positive"] / row["returned"] if (row["returned"]) else None
            )
    return result

def model_level(base, vector, eligible, calls):
    row = {**base, **{field: Fraction() for (field) in (LEVEL_FIELDS)}}
    for (gene, weight) in (vector.items()):
        if (gene not in eligible):
            continue
        state, label = calls[(base["sample_id"], gene)]
        row["returned"] += weight
        if (state == "MEASURED"):
            row["positive" if (label == 1) else "zero"] += weight
        else:
            row[UNKNOWN_FIELDS[state]] += weight
            row["unknown"] += weight
    row["pooled_positive_fraction"] = (
        row["positive"] / row["returned"] if (row["returned"]) else None
    )
    require(
        (row["returned"] == row["positive"] + row["zero"] + row["unknown"]),
        "Category mass conservation",
    )
    return row

def derive(models, policies, calls, categories, context):
    counts = validate_models(models, context)
    model_levels, model_contrasts, coefficients = [], [], []
    methods = ("PRODIGY", "DawnRank") if (context == "CRC") else ("DawnRank",)
    for (model) in (models):
        mid, sample = model["model_id"], model["sample_id"]
        identity = {key: model[key] for (key) in (("model_id", "sample_id", "donor_id"))}
        identity["context"] = context
        degree_update = combine(
            [(1, policies[mid]["degree"]["G1"]), (-1, policies[mid]["degree"]["G0"])]
        )
        contrasts = {("degree", "degree_update"): degree_update}
        for (method) in (methods):
            joint = combine(
                [(1, policies[mid][method]["G1H1"]), (-1, policies[mid][method]["G0H0"])]
            )
            contrasts[(method, "joint")] = joint
            contrasts[(method, "Theta")] = combine([(1, joint), (-1, degree_update)])
        for (partition, names) in (CATEGORIES.items()):
            for (category) in (names):
                eligible = {
                    gene
                    for (gene) in (model["query"])
                    if (partition == "all" or categories[partition][gene] == category)
                }
                base = {**identity, "partition": partition, "category": category}
                for (policy, states) in (policies[mid].items()):
                    for (state, vector) in (states.items()):
                        model_levels.append(
                            model_level(
                                {**base, "policy": policy, "state": state},
                                vector,
                                eligible,
                                calls,
                            )
                        )
                for ((policy, contrast), vector) in (contrasts.items()):
                    row = {
                        **base,
                        "policy": policy,
                        "contrast": contrast,
                        **{field: Fraction() for (field) in (DELTA_FIELDS)},
                    }
                    for (gene) in (sorted(eligible & set(vector))):
                        coefficient = vector[gene]
                        item = contribution(coefficient, calls[(sample, gene)])
                        for (field) in (DELTA_FIELDS):
                            row[field] += item[field]
                        if (partition == "all"):
                            coefficients.append(
                                {
                                    **identity,
                                    "policy": policy,
                                    "contrast": contrast,
                                    "gene": gene,
                                    "core": categories["core"][gene],
                                    "role": categories["role"][gene],
                                    "label_state": calls[(sample, gene)][0],
                                    "label_value": calls[(sample, gene)][1],
                                    "coefficient": coefficient,
                                    "population_coefficient": coefficient
                                    / (len(counts) * counts[model["donor_id"]]),
                                    **item,
                                }
                            )
                    model_contrasts.append(finish_delta(row))
    donor_levels = aggregate_rows(model_levels, LEVEL_FIELDS, counts, "donor")
    population_levels = aggregate_rows(model_levels, LEVEL_FIELDS, counts, "population")
    donor_contrasts = aggregate_rows(model_contrasts, DELTA_FIELDS, counts, "donor")
    population_contrasts = aggregate_rows(model_contrasts, DELTA_FIELDS, counts, "population")
    oncogenes = sorted(
        {
            gene
            for (model) in (models)
            for (gene) in (model["query"])
            if (categories["role"][gene] == "ONCOGENE")
        }
    )
    gene_rows, concentration = [], []
    contrast_keys = [("degree", "degree_update")] + [
        (method, contrast) for (method) in (methods) for (contrast) in (("joint", "Theta"))
    ]
    for (policy, contrast) in (contrast_keys):
        by_gene = {gene: {field: Fraction() for (field) in (DELTA_FIELDS)} for (gene) in (oncogenes)}
        for (row) in (coefficients):
            if (row["role"] == "ONCOGENE" and (row["policy"], row["contrast"]) == (
                policy,
                contrast,
            )):
                for (field) in (DELTA_FIELDS):
                    by_gene[row["gene"]][field] += row[field] / (
                        len(counts) * counts[row["donor_id"]]
                    )
        gene_rows.extend(
            finish_delta(
                {
                    "context": context,
                    "policy": policy,
                    "gene": gene,
                    "contrast": contrast,
                    **values,
                }
            )
            for (gene, values) in (by_gene.items())
        )
    donor_rows = [
        row
        for (row) in (donor_contrasts)
        if (row["partition"] == "role" and row["category"] == "ONCOGENE")
    ]
    for (policy, contrast) in (contrast_keys):
        for (unit, source, identifier, divisor) in ((
            ("gene", gene_rows, "gene", 1),
            ("donor", donor_rows, "donor_id", len(counts)),
        )):
            for (metric) in (("K", "returned_change", "unknown_change")):
                vector = {
                    row[identifier]: row[metric] / divisor
                    for (row) in (source)
                    if ((row["policy"], row["contrast"]) == (policy, contrast))
                }
                if (vector):
                    concentration.append(
                        {
                            "context": context,
                            "policy": policy,
                            "contrast": contrast,
                            "unit": unit,
                            "metric": metric,
                            **summarise(vector),
                        }
                    )
    return {
        "model_category_levels": model_levels,
        "donor_category_levels": donor_levels,
        "population_category_levels": population_levels,
        "model_category_contrasts": model_contrasts,
        "donor_category_contrasts": donor_contrasts,
        "population_category_contrasts": population_contrasts,
        "model_gene_coefficients": coefficients,
        "oncogene_gene_contributions": gene_rows,
        "oncogene_donor_contributions": donor_rows,
        "oncogene_concentration": concentration,
    }
