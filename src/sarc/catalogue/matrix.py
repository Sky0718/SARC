import math
from fractions import Fraction

from ..transfer import scoring as core
from .constants import (
    CONTEXT,
    ENDPOINTS,
    METHODS,
    PARAMETERS,
    REFERENCES,
    SEEDS,
    group_key,
    network_states,
)

def validate_registry(registry, population):
    if (
        len(population) != 396
        or len(set(population)) != 396
        or any(not isinstance(sample, str) or not sample for (sample) in (population))
    ):
        raise ValueError("The exact complete unique 396-sample population is required")
    expected = {
        (method, network, PARAMETERS[method], seed)
        for (method) in (METHODS)
        for (network) in (network_states(method))
        for (seed) in (SEEDS if (method == "PRODIGY") else (None,))
    }
    observed = []
    for (cell_id, cell) in (registry.items()):
        if (
            cell_id != cell["cell_id"]
            or cell["node"] != "NCS-S06"
            or cell["context"] != CONTEXT
            or cell["sample_order"] != list(population)
            or cell.get("tuning_default_alias") is not False
        ):
            raise ValueError("Foreign cell, altered population or protected retuning")
        observed.append(
            (
                cell["method"],
                cell["network_id"],
                cell["parameter_id"],
                cell["master_seed"],
            )
        )
    if (
        len(registry) != 69
        or len(observed) != len(expected)
        or set(observed) != expected
    ):
        raise ValueError("The complete frozen 69-cell portfolio is required")
    return expected

def validate_references(references, population):
    if (list(references) != list(population)):
        raise ValueError("Reference population differs from the frozen ordered cohort")
    for (sample, row) in (references.items()):
        if (row["sample"] != sample):
            raise ValueError("Reference row identity differs")
        supports = {}
        for (support) in (("common", "native")):
            eligible = row[support + "_eligible_literal_gene_labels"]
            if (len(eligible) != len(set(eligible))):
                raise ValueError("Repeated eligible label")
            supports[support] = set(eligible)
            catalogues = row[support + "_reference_positive_labels"]
            if (set(catalogues) != set(REFERENCES)):
                raise ValueError("The four frozen reference catalogues are required")
            for (values) in (catalogues.values()):
                if (len(values) != len(set(values)) or not set(values).issubset(
                    supports[support]
                )):
                    raise ValueError("Invalid fixed reference eligibility")
        if (not supports["common"].issubset(supports["native"])):
            raise ValueError("Common eligibility exceeds valid original mutations")

def score_cell(cell, raw, references):
    population = cell["sample_order"]
    if (list(raw) != population or list(references) != population):
        raise ValueError("Rankings and references do not retain exact sample order")
    samples = {}
    for (sample) in (population):
        row = raw[sample]
        state = row["method_status"]
        if (state not in ("SUCCESS", "SUCCESS_EMPTY", "FAILED", "UNAVAILABLE")):
            raise ValueError("Unknown method state")
        valid = state in ("SUCCESS", "SUCCESS_EMPTY")
        genes, scores = row["genes"], row["scores"]
        if ((valid and (genes is None or scores is None)) or (
            not valid and (genes is not None or scores is not None)
        )):
            raise ValueError(
                "Missing outputs cannot be replaced with successful empty lists"
            )
        if (valid and (
            bool(genes) != (state == "SUCCESS")
            or len(genes) != len(scores)
            or len(set(genes)) != len(genes)
        )):
            raise ValueError("Native candidate roster is incomplete or duplicated")
        source = references[sample]
        if (valid and not set(genes).issubset(
            source["native_eligible_literal_gene_labels"]
        )):
            raise ValueError(
                "Prediction exceeds frozen valid original mutation eligibility"
            )
        if (valid and any(
            isinstance(score, bool) or not math.isfinite(score) for (score) in (scores)
        )):
            raise ValueError("Invalid native numerical score")
        measurements = {}
        for (support, catalogue, k) in (ENDPOINTS):
            eligible = source[support + "_eligible_literal_gene_labels"]
            positives = source[support + "_reference_positive_labels"][catalogue]
            measurements[(support, catalogue, k)] = core.score_list(
                genes, scores, eligible, positives, k
            )
        samples[sample] = {
            "method_status": state,
            "reason": row.get("reason", ""),
            "genes": genes,
            "scores": scores,
            "measurements": measurements,
        }
    return {"cell": cell, "samples": samples}

def reduce_cells(completed, registry, references):
    population = list(references)
    validate_registry(registry, population)
    validate_references(references, population)
    if (set(completed) != set(registry)):
        raise ValueError("Incomplete frozen method matrix")
    repeated = {}
    for (cell_id, value) in (completed.items()):
        cell = registry[cell_id]
        if (value["cell"] != cell or list(value["samples"]) != population):
            raise ValueError("Scored cell metadata changed")
        for (sample) in (population):
            if (set(value["samples"][sample]["measurements"]) != set(ENDPOINTS)):
                raise ValueError("A required endpoint was dropped")
        key = group_key(cell["method"], cell["network_id"])
        repeated.setdefault(key, {})[cell["master_seed"]] = value
    groups = {}
    for (key, repeats) in (repeated.items()):
        seeds = SEEDS if (key[1] == "PRODIGY") else (None,)
        if (set(repeats) != set(seeds)):
            raise ValueError("Five paired seeds must precede sample reduction")
        measurements = {}
        for (endpoint) in (ENDPOINTS):
            support, reference, k = endpoint
            rows = {}
            for (sample) in (population):
                values = [
                    repeats[seed]["samples"][sample]["measurements"][endpoint]
                    for (seed) in (seeds)
                ]
                source = references[sample]
                eligible = frozenset(source[support + "_eligible_literal_gene_labels"])
                positive = frozenset(
                    source[support + "_reference_positive_labels"][reference]
                )
                if (any(
                    (value["eligible_count"], value["positive_count"], value["k"])
                    != (len(eligible), len(positive), k)
                    for (value) in (values)
                )):
                    raise ValueError("Seed-specific endpoint support changed")
                summary = (
                    core.repeat_summary(values)
                    if (len(seeds) == 5)
                    else dict(values[0])
                )
                hits = summary["hits"]

                def average(name):
                    entries = [value[name] for (value) in (values)]
                    return (
                        None
                        if (any(value is None for (value) in (entries)))
                        else sum(entries, Fraction(0)) / len(entries)
                    )

                summary.update(
                    {
                        "eligible_ids": eligible,
                        "positive_ids": positive,
                        "eligible_count": len(eligible),
                        "positive_count": len(positive),
                        "k": k,
                        "recall": hits / len(positive)
                        if (hits is not None and positive)
                        else None,
                        "repeat_count": len(seeds),
                        "complete_repeat_count": sum(
                            value["hits"] is not None for (value) in (values)
                        ),
                        "native_output_count": average("native_output_count"),
                        "common_output_count": average("common_output_count"),
                        "positive_output_count": average("positive_output_count"),
                    }
                )
                rows[sample] = summary
            measurements[endpoint] = rows
        groups[key] = {
            "population": population,
            "expected_seeds": seeds,
            "measurements": measurements,
            "cohort_summaries": {
                endpoint: core.mean_records(list(rows.values()))
                for ((endpoint, rows)) in (measurements.items())
            },
        }
    return groups

def recall_summary(rows):
    positive = [row for (row) in (rows.values()) if (row["positive_count"] > 0)]
    if (not positive):
        return {
            "full_population": len(rows),
            "positive_eligible_population": 0,
            "value": None,
            "status": "NO_REFERENCE_POSITIVES",
        }
    records = [
        {
            "hits": row["recall"],
            "lower": row["lower"] / row["positive_count"],
            "upper": row["upper"] / row["positive_count"],
        }
        for (row) in (positive)
    ]
    return {
        "full_population": len(rows),
        "positive_eligible_population": len(positive),
        "status": "FIXED_POSITIVE_ELIGIBILITY_SUBSET",
        **core.mean_records(records),
    }
