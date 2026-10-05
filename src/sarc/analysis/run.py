import csv
import gzip
import io
import json
from fractions import Fraction
from pathlib import Path
from . import (
    contextual_budget,
    contextual_order,
    cross_method,
    dawn_factorial,
    domain_bridge,
    factorial_decisions,
    hypothesis_analysis,
    hypothesis_census,
    hypothesis_reader,
    native_context,
    prior_value,
    prodigy_factorial,
    release_consequence,
    robustness,
    source_prior,
    upstream_selected,
)

def object_pairs(pairs):
    output = {}
    for (key, value) in (pairs):
        if (key in output):
            raise ValueError("Duplicate JSON key: " + key)
        output[key] = value
    return output

def reject_constant(value):
    raise ValueError("Nonfinite JSON token: " + value)

def read_input(base, item):
    specification = {"path": item, "format": "json"} if (isinstance(item, str)) else item
    if ("paths" in specification):
        return [read_input(base, value) for (value) in (specification["paths"])]
    if ("mapping" in specification):
        return {key: read_input(base, value) for (key, value) in (specification["mapping"].items())}
    path = (base / specification["path"]).resolve()
    opener = gzip.open if (path.suffix == ".gz") else open
    with opener(path, "rt", encoding = "utf-8-sig") as stream:
        kind = specification.get("format", "json")
        if (kind == "text"):
            return stream.read()
        if (kind in ("csv", "tsv")):
            return list(csv.DictReader(stream, delimiter = "\t" if (kind == "tsv") else ","))
        if (kind == "masked_binary_ledger"):
            return [
                release_consequence.masked_row(line) for (line) in (stream) if (line.strip())
            ]
        if (kind == "jsonl"):
            return [
                json.loads(line, object_pairs_hook = object_pairs, parse_constant = reject_constant)
                for (line) in (stream)
                if (line.strip())
            ]
        if (kind != "json"):
            raise ValueError("Unsupported declared input format")
        result = json.load(
            stream, object_pairs_hook = object_pairs, parse_constant = reject_constant
        )
        for (key) in (specification.get("select", [])):
            result = result[key]
        return result

def serializable(value):
    if (isinstance(value, Fraction)):
        return {"numerator": value.numerator, "denominator": value.denominator}
    if (isinstance(value, dict)):
        if (any(not isinstance(key, str) for (key) in (value))):
            raise ValueError("Output requires explicit string keys or row records")
        return {key: serializable(item) for ((key, item)) in (value.items())}
    if (isinstance(value, (list, tuple))):
        return [serializable(item) for (item) in (value)]
    if (isinstance(value, set)):
        return [serializable(item) for (item) in (sorted(value))]
    return value

def run_operation(config_path, registry):
    config_path = Path(config_path).resolve()
    config = read_input(config_path.parent, config_path.name)
    if (config["operation"] not in registry):
        raise ValueError("Unknown scientific operation")
    inputs = config.get("inputs", {})
    options = config.get("options", {})
    if (set(inputs) & set(options)):
        raise ValueError("An argument cannot be both a file input and an option")
    arguments = {
        name: read_input(config_path.parent, item) for ((name, item)) in (inputs.items())
    }
    arguments.update(options)
    if (config["operation"] == "prodigy_checkpoint_read"):
        arguments["base_directory"] = str(
            (config_path.parent / arguments.get("base_directory", ".")).resolve()
        )
    output = (config_path.parent / config["output"]).resolve()
    if (output.exists()):
        raise FileExistsError("Use a distinct output path; existing results are retained")
    result = registry[config["operation"]](**arguments)
    output.parent.mkdir(parents = True, exist_ok = True)
    with output.open("x", encoding = "utf-8", newline = "\n") as stream:
        json.dump(serializable(result), stream, ensure_ascii = False, allow_nan = False, indent = 2)
        stream.write("\n")
    return {"operation": config["operation"], "output": str(output)}

def contextual_order_analysis(ranks, eligible, models):
    robustness.validate_split(
        {
            "models": models,
            "ordered_donors": sorted({row["study_donor_id"] for (row) in (models)}),
        }
    )
    by_model = {row["sample_ID"]: row for (row) in (models)}
    if (len(by_model) != 85 or len({row["study_donor_id"] for (row) in (models)}) != 83):
        raise ValueError("Complete 85-model, 83-donor population required")
    groups = {}
    for (row) in (ranks):
        method, pool, release, seed, sample = row["key"]
        if (pool != "all34"):
            continue
        observations = groups.setdefault((method, release, seed), {})
        if (sample in observations or sample not in by_model):
            raise ValueError("Repeated or unknown ranking identity")
        menu = eligible["pools"][pool][sample]["native_eligible_literal_gene_labels"]
        observations[sample] = contextual_order.allocation(
            row["genes"], row["scores"], set(menu), 10
        )
    expected = {
        (method, release, seed)
        for (method) in (("DawnRank", "PersonaDrive", "PRODIGY"))
        for (release) in (("native_11_0", "native_11_5", "native_12_0"))
        for (seed) in (
            (104729, 130363, 155921, 196613, 228017) if (method == "PRODIGY") else (None,)
        )
    }
    if (set(groups) != expected or any(
        set(group) != set(by_model) for (group) in (groups.values())
    )):
        raise ValueError("Complete method, release, seed and model grid required")
    output = []
    donors = {sample: row["study_donor_id"] for ((sample, row)) in (by_model.items())}
    for (key, observations) in (sorted(groups.items(), key = lambda item: str(item[0]))):
        for (population) in (("all85", "development", "held_out")):
            chosen = {
                sample: weights
                for ((sample, weights)) in (observations.items())
                if (population == "all85" or by_model[sample]["split"] == population)
            }
            result = contextual_order.census(chosen, donors)
            output.append(
                {
                    "method": key[0],
                    "release": key[1],
                    "seed": key[2],
                    "population": population,
                    **result,
                }
            )
    return {"groups": output}

def robustness_analysis(rows, split, bootstrap, production = True):
    groups, accounting = robustness.accumulate_ledger(rows, split, production)
    return {
        "accounting": accounting,
        "analysis": robustness.analyse(groups, split, bootstrap),
    }

def release_analysis(rows, split, production = True):
    groups, labels, accounting = release_consequence.accumulate(rows, split, production)
    return {
        "accounting": accounting,
        "analysis": release_consequence.analyse(groups, labels, split),
    }

def source_graph(gml, sample, mapping, pathways):
    import networkx as nx

    statistics, interface = source_prior.graph_statistics(
        io.BytesIO(gml.encode("utf-8")),
        sample,
        mapping,
        {gene: set(values) for ((gene, values)) in (pathways.items())},
        nx,
    )
    return {"statistics": statistics, "interface": interface}

def source_model(sample, release, scores, statistics, coefficients, eligibility):
    stream = io.StringIO()
    record = source_prior.model_record(
        sample,
        release,
        scores,
        statistics,
        source_prior.dyadic_row(coefficients),
        eligibility,
        stream,
    )
    return {
        "model": record,
        "genes": [json.loads(line) for (line) in (stream.getvalue().splitlines())],
    }

def source_models(records):
    releases = ("native_11_0", "native_11_5", "native_12_0")
    grouped = {release: [] for (release) in (releases)}
    genes, references = [], {release: {} for (release) in (releases)}
    for (row) in (records):
        if (row["release"] not in grouped):
            raise ValueError("Unexpected source-prior release")
        source_prior.recipient_invariance(
            references[row["release"]], row["statistics"], row["sample"]
        )
        result = source_model(**row)
        grouped[row["release"]].append(result["model"])
        genes.extend(result["genes"])
    if (any(
        {row["sample_id"] for (row) in (grouped[release])}
        != {row["sample_id"] for (row) in (grouped[releases[0]])}
        for (release) in (releases)
    )):
        raise ValueError("Source-prior model populations differ across releases")
    summaries = [
        source_prior.release_summary(release, grouped[release], 85) for (release) in (releases)
    ]
    return {
        "models": [row for (release) in (releases) for (row) in (grouped[release])],
        "genes": genes,
        "summaries": summaries,
    }

def prior_prepare(carriers, models, table5_source, production = True):
    prior_value.validate_models(models, production)
    records = [
        prior_value.carrier_record(row, row["release"], row["sample_id"])
        for (row) in (carriers)
    ]
    expected = {
        (release, model["sample_ID"])
        for (release) in (("native_11_0", "native_11_5", "native_12_0"))
        for (model) in (models)
    }
    observed = {(row["release"], row["sample_ID"]) for (row) in (records)}
    if (len(records) != len(expected) or observed != expected):
        raise ValueError("Full source-prior release and model population required")
    return {
        "prelabel_freeze": True,
        "models": models,
        "records": records,
        "table5_source": table5_source,
        "requested_labels": [
            {"sample_id": sample, "gene": gene}
            for ((sample, gene)) in (prior_value.required_cells(records))
        ],
    }

def prior_evaluate(frozen, selected):
    if (
        frozen.get("prelabel_freeze") is not True
        or frozen["table5_source"] != selected["source_binding"]
    ):
        raise ValueError("Fixed allocation and matching assay source required")
    records = [
        {
            **row,
            "coefficients": None
            if (row["coefficients"] is None)
            else {
                gene: prior_value.fraction(value)
                for ((gene, value)) in (row["coefficients"].items())
            },
        }
        for (row) in (frozen["records"])
    ]
    calls = {}
    for (row) in (selected["calls"]):
        native_context.validate_call(row)
        key = (row["sample_id"], row["gene"])
        if (key in calls):
            raise ValueError("Repeated selected assay identity")
        calls[key] = {
            "mask": {
                "MEASURED": "PRESENT",
                "EXPLICIT_NULL": "EXPLICIT_NULL",
                "ABSENT": "KEY_ABSENT",
            }[row["state"]],
            "call": row["value"],
        }
    if (set(calls) != {
        (row["sample_id"], row["gene"]) for (row) in (frozen["requested_labels"])
    }):
        raise ValueError("Frozen selected-call support differs")
    output = []
    for (release) in (("native_11_0", "native_11_5", "native_12_0")):
        for (subset) in (("held_out", "development", "all")):
            models = [
                row
                for (row) in (frozen["models"])
                if (subset == "all" or row["split"] == subset)
            ]
            ids = {row["sample_ID"] for (row) in (models)}
            chosen = [
                row
                for (row) in (records)
                if (row["release"] == release and row["sample_ID"] in ids)
            ]
            output.append(prior_value.summarise(chosen, models, release, subset, calls))
    return {"summaries": output}

OPERATIONS = {
    "prodigy_freeze": prodigy_factorial.freeze,
    "prodigy_evaluate": prodigy_factorial.evaluate,
    "factorial_decision_freeze": factorial_decisions.freeze,
    "factorial_decision_evaluate": factorial_decisions.evaluate,
    "dawn_freeze": dawn_factorial.build_freeze,
    "dawn_evaluate": dawn_factorial.evaluate,
    "cross_method": cross_method.analysis,
    "bridge_freeze": domain_bridge.prepare_bridge,
    "bridge_evaluate": domain_bridge.finish_from_prepared,
    "native_context_freeze": native_context.prepare,
    "native_context_merge_labels": native_context.merge_selected,
    "native_context_evaluate": native_context.evaluate,
    "contextual_order": contextual_order_analysis,
    "contextual_budgets": contextual_budget.analyse,
    "robustness": robustness_analysis,
    "release_consequence": release_analysis,
    "source_prior_graph": source_graph,
    "source_prior_model": source_model,
    "source_prior": source_models,
    "prior_value_freeze": prior_prepare,
    "prior_value_evaluate": prior_evaluate,
    "prodigy_checkpoint_read": hypothesis_reader.read_rds,
    "prodigy_retained_hypotheses": hypothesis_analysis.analyse,
    "prodigy_hypothesis_census": hypothesis_census.analyse,
    "prodigy_upstream_selected": upstream_selected.analyse,
}

def run(config_path):
    return run_operation(config_path, OPERATIONS)
