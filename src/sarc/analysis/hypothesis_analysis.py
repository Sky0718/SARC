import struct
from collections import Counter
from . import hypothesis_codec as codec
from . import retained_hypotheses as retained
from .contextual_order import allocation

RELEASES = retained.RELEASES
PAIRS = retained.PAIRS
SEEDS = (104729, 130363, 155921, 196613, 228017)
BUDGETS = retained.BUDGETS
MODES = (
    "original_full_pipeline",
    "deterministic_stage_reuse",
    "original_after_preserved_reuse_assertion",
)

def validate_roster(models):
    samples = [row["sample_ID"] for (row) in (models)]
    codec.require(
        len(models) == len(set(samples)) == 85
        and len({row["model_id"] for (row) in (models)}) == 85
        and len({row["study_donor_id"] for (row) in (models)}) == 83,
        "Complete 85-model and 83-donor population required",
    )
    return samples

def source_index(origins):
    result = {}
    for (row) in (origins):
        identity = row["source_id"]
        codec.require(
            isinstance(identity, str) and identity and identity not in result,
            "Unique explicit checkpoint source identities required",
        )
        codec.require(
            row["status"] in ("AVAILABLE", "UNAVAILABLE"),
            "Explicit checkpoint availability required",
        )
        codec.require(
            (row.get("checkpoint") is not None) == (row["status"] == "AVAILABLE"),
            "Checkpoint availability and payload differ",
        )
        if (row["status"] != "AVAILABLE"):
            codec.require(
                isinstance(row.get("reason"), str) and row["reason"],
                "Unavailable source requires its reason",
            )
        result[identity] = row
    return result

def state_index(models, seed_states, origins):
    samples = validate_roster(models)
    expected = {
        (sample, release, seed)
        for (sample) in (samples)
        for (release) in (RELEASES)
        for (seed) in (SEEDS)
    }
    result, used = {}, set()
    for (row) in (seed_states):
        key = (row["sample"], row["release"], row["master_seed"])
        codec.require(
            key in expected and key not in result, "Repeated or unplanned seed state"
        )
        codec.require(
            row["status"] in ("AVAILABLE", "UNAVAILABLE", "FAILED"),
            "Explicit seed availability required",
        )
        origin = row.get("checkpoint_source_id")
        if (origin is not None):
            codec.require(origin in origins, "Checkpoint source identity missing")
            source = origins[origin]
            codec.require(
                source["sample"] == row["sample"]
                and source["release"] == row["release"],
                "Checkpoint source belongs to another model or release",
            )
            codec.require(row["mode"] in MODES, "Unknown checkpoint ancestry mode")
            expected_seed = (
                SEEDS[0]
                if (row["mode"] == "deterministic_stage_reuse")
                else row["master_seed"]
            )
            codec.require(
                source["checkpoint_source_seed"]
                == row["checkpoint_source_seed"]
                == expected_seed,
                "Checkpoint seed ancestry differs",
            )
            used.add(origin)
        if (row["status"] == "AVAILABLE"):
            codec.require(
                origin is not None and origins[origin]["status"] == "AVAILABLE",
                "Available state lacks its available checkpoint",
            )
        else:
            codec.require(
                isinstance(row.get("reason"), str) and row["reason"],
                "Unavailable seed state requires its reason",
            )
        result[key] = row
    codec.require(
        set(result) == expected and used == set(origins),
        "Full seed grid and exact source population required",
    )
    return result

def expected_identity(source):
    return {
        "context": "COLO_ORGANOID_SANGER",
        "sample": source["sample"],
        "network_id": source["release"],
        "normal_pool": "all34",
        "master_seed": source["checkpoint_source_seed"],
        "alpha": 0.05,
    }

def rank_index(models, rankings, eligibility):
    samples = set(validate_roster(models))
    codec.require(set(eligibility) == samples, "Complete model eligibility required")
    for (genes) in (eligibility.values()):
        codec.unique_axis(genes, "Eligible mutation")
    expected = {
        (sample, release, seed)
        for (sample) in (samples)
        for (release) in (RELEASES)
        for (seed) in (SEEDS)
    }
    result = {}
    for (row) in (rankings):
        key = (row["sample"], row["release"], row["seed"])
        codec.require(
            key in expected and key not in result, "Unplanned or repeated ranking"
        )
        genes, scores = row["genes"], row["scores"]
        if (genes is None):
            codec.require(scores is None, "Unavailable ranks require null scores")
        else:
            codec.require(len(genes) == len(scores), "Native score axis differs")
            bits = row.get("scores_little_endian_float64_hex")
            if (bits is not None):
                codec.require(
                    bits == [struct.pack("<d", value).hex() for (value) in (scores)],
                    "Native binary64 ranks differ from recorded score bytes",
                )
        result[key] = {
            "genes": genes,
            "scores": scores,
            "weights": {
                budget: allocation(
                    genes, scores, set(eligibility[row["sample"]]), budget
                )
                for (budget) in (BUDGETS)
            },
        }
    codec.require(set(result) == expected, "All 1275 native ranking states required")
    return result

def derive_states(models, seed_states, origins, rankings, eligibility, results):
    origin_map = source_index(origins)
    states = state_index(models, seed_states, origin_map)
    ranks = rank_index(models, rankings, eligibility)
    decoded, derived, accounting, checks = {}, {}, [], []
    used_results = set()
    for (key, state) in (states.items()):
        sample, release, seed = key
        rank = ranks[key]
        source_id = state.get("checkpoint_source_id")
        value = {
            "status": "UNAVAILABLE_HYPOTHESIS_DERIVATION",
            "weights": rank["weights"],
            "lineage": state,
        }
        reason = state.get("reason")
        if (state["status"] == "FAILED"):
            codec.require(
                rank["genes"] is None and rank["scores"] is None,
                "Failed native state has a rank",
            )
            value["status"] = "UNAVAILABLE_NATIVE_FAILED"
        elif (state["status"] == "AVAILABLE"):
            result_id = state["result_source_id"]
            codec.require(
                isinstance(result_id, str) and result_id in results,
                "Explicit native result source required",
            )
            used_results.add(result_id)
            try:
                native = codec.decode(results[result_id])
                record = native["result"]
                codec.require(
                    codec.scalar(record["sample"]) == sample
                    and codec.scalar(record["network_id"]) == release
                    and codec.scalar(record["master_seed"]) == seed,
                    "Final result model, release or seed identity differs",
                )
                genes = codec.listify(record["genes"])
                scores = (
                    record["scores"].values
                    if (isinstance(record["scores"], codec.Numeric))
                    else []
                )
                codec.require(
                    genes == rank["genes"] and scores == rank["scores"],
                    "Final native result differs from ranking input",
                )
                if (source_id not in decoded):
                    source = origin_map[source_id]
                    decoded[source_id] = codec.checkpoint(
                        codec.decode(source["checkpoint"]), expected_identity(source)
                    )
                checkpoint = decoded[source_id]
                aggregate = codec.aggregate_check(checkpoint, native)
                codec.require(
                    set(genes) <= set(checkpoint["genes"]),
                    "Emitted target absent from checkpoint mutation axis",
                )
                checks.append(
                    {
                        "sample": sample,
                        "release": release,
                        "master_seed": seed,
                        "checkpoint_source_id": source_id,
                        "result_source_id": result_id,
                        "checks": aggregate,
                    }
                )
                value.update(status = "AVAILABLE", checkpoint = checkpoint)
            except (KeyError, ValueError, TypeError, struct.error) as error:
                reason = type(error).__name__ + ": " + str(error)
        derived[key] = value
        accounting.append(
            {
                "sample": sample,
                "release": release,
                "seed": seed,
                "status": value["status"],
                "reason": reason,
                "checkpoint_source_id": source_id,
                "checkpoint_source_seed": state.get("checkpoint_source_seed"),
                "mode": state.get("mode"),
            }
        )
    codec.require(
        used_results == set(results), "Extra unreferenced native result payload"
    )
    return derived, accounting, checks

def analyse(models, rankings, eligibility, origins, seed_states, results, pathways):
    states, accounting, aggregate_checks = derive_states(
        models, seed_states, origins, rankings, eligibility, results
    )
    dictionary = codec.pathway_dictionary(codec.decode(pathways))
    weights = retained.model_weights(models)
    records, contrasts, rectangles = [], [], []
    for (model) in (models):
        sample = model["sample_ID"]
        cache = {}
        for (old_release, new_release) in (PAIRS):
            for (seed) in (SEEDS):
                old, new = (
                    states[sample, old_release, seed],
                    states[sample, new_release, seed],
                )
                prefix = "|".join((sample, old_release, new_release, str(seed)))
                if (old["status"] == new["status"] == "AVAILABLE"):
                    rectangles.append(
                        dict(
                            retained.common_rectangle(
                                old["checkpoint"], new["checkpoint"]
                            ),
                            identity = prefix,
                        )
                    )

                def contrast(gene):
                    key = (
                        old_release,
                        new_release,
                        old["lineage"]["checkpoint_source_id"],
                        new["lineage"]["checkpoint_source_id"],
                        gene,
                    )
                    if (key not in cache):
                        value = retained.gene_contrast(
                            old["checkpoint"], new["checkpoint"], gene, dictionary
                        )
                        value["contrast_id"] = prefix + "|" + gene
                        cache[key] = value
                        contrasts.append(value)
                    return cache[key]

                for (budget) in (BUDGETS):
                    records.append(
                        retained.pair_record(
                            old,
                            new,
                            budget,
                            {
                                "sample": sample,
                                "model_id": model["model_id"],
                                "donor": model["study_donor_id"],
                                "model_weight": weights[sample],
                                "release_pair": [old_release, new_release],
                                "seed": seed,
                            },
                            contrast,
                        )
                    )
    codec.require(
        len(accounting) == 1275 and len(records) == 5100,
        "Finite diagnostic population incomplete",
    )
    return {
        "seed_accounting": accounting,
        "aggregate_checks": aggregate_checks,
        "pair_seed_budget": records,
        "gene_contrasts": contrasts,
        "common_rectangles": rectangles,
        "summaries": retained.summarize(records, weights),
        "counts": {
            "models": 85,
            "donors": 83,
            "tasks": 255,
            "seed_states": 1275,
            "pair_seed_budget_records": 5100,
            "seed_status_counts": dict(
                Counter(row["status"] for (row) in (accounting))
            ),
        },
        "functional_labels_used": False,
    }
