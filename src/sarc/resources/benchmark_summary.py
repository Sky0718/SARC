from collections import defaultdict
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from .comparators import STRATEGIES
from .io import check, read_json, read_rows, write_json
from .benchmark_scenarios import EVALUATION_SEEDS, scenario_specs

METRICS = (
    "correct_reported_answer_fraction",
    "reported_answer_coverage",
    "ineligible_answer_fraction",
    "correct_abstention_fraction",
    "support_denominator_mismatch_fraction",
    "false_alert_fraction",
    "missed_event_fraction",
)

def opportunity(row, metric):
    eligible = row["oracle_eligible"]
    answered = row["answer_state"] == "ANSWERED"
    if (metric == "correct_reported_answer_fraction"):
        return int(row["correct_answer"]), int(eligible and answered)
    if (metric == "reported_answer_coverage"):
        return int(eligible and answered), int(eligible)
    if (metric == "ineligible_answer_fraction"):
        return int(not eligible and answered), int(not eligible)
    if (metric == "correct_abstention_fraction"):
        return int(not eligible and row["answer_state"] == "ABSTAINED"), int(
            not eligible
        )
    if (metric == "support_denominator_mismatch_fraction"):
        return int(answered and row["denominator_mismatch"]), int(answered)
    event_task = row["task"] in {
        "fixed_leading_set_change",
        "native_leading_set_change",
    }
    if (metric == "false_alert_fraction"):
        eligible = (
            event_task and eligible and answered and row["truth_boolean"] is False
        )
        return int(eligible and row["answer_boolean"] is True), int(eligible)
    if (metric == "missed_event_fraction"):
        eligible = event_task and eligible and answered and row["truth_boolean"] is True
        return int(eligible and row["answer_boolean"] is False), int(eligible)
    raise ValueError("Unregistered metric")

def closed_fraction(numerator, denominator):
    if (denominator == 0):
        return {
            "value": None,
            "numerator": None,
            "denominator": None,
            "support_denominator": 0,
            "reason": "NO_ELIGIBLE_OPPORTUNITIES",
        }
    return {
        "value": numerator / denominator,
        "numerator": numerator,
        "denominator": denominator,
        "support_denominator": denominator,
        "reason": None,
    }

def descriptor_index(population):
    required = {"disease_id", "size_stratum", "leader_stratum", "native_leader_margin"}
    check(
        population and all(required <= set(row) for (row) in (population)),
        "Complete sampled-panel descriptors are required",
    )
    identifiers = [row["disease_id"] for (row) in (population)]
    check(len(set(identifiers)) == len(identifiers), "Duplicate sampled panel")
    check(
        all(
            row["size_stratum"] in ("30_99", "100_499", "500_plus")
            and row["leader_stratum"] in ("unique", "tied")
            for (row) in (population)
        ),
        "Unknown baseline stratum",
    )
    return (
        identifiers,
        {panel: index for ((index, panel)) in (enumerate(identifiers))},
        {row["disease_id"]: row for (row) in (population)},
    )

def aggregate_tasks(source, population):
    panel_ids, panel_indices, descriptors = descriptor_index(population)
    groups = {}
    group_panels = defaultdict(set)
    paired = {}
    contrast_counts = defaultdict(lambda: np.zeros((len(panel_ids), 5), dtype = np.int64))
    identities = set()
    observed_panels = set()
    instances = {}
    specifications = defaultdict(set)
    task_names = {
        "fixed_leading_set_change",
        "native_leading_set_change",
        "persistent_versus_roster_membership_description",
        "mapping_ambiguity_requires_abstention",
    }
    for (row) in (read_rows(source)):
        check(row["phase"] in ("development", "evaluation"), "Unknown benchmark phase")
        if (row["phase"] != "evaluation"):
            continue
        check(
            row["source_panel_id"] in descriptors
            and row["strategy"] in STRATEGIES
            and row["task"] in task_names,
            "Unknown evaluation panel, strategy or task",
        )
        metadata = (row["source_panel_id"], row["scenario"], row["level"], row["seed"])
        if (row["instance_id"] in instances):
            check(
                instances[row["instance_id"]] == metadata,
                "Instance metadata differs across task rows",
            )
        else:
            check(
                metadata[1:] not in specifications[metadata[0]],
                "Repeated evaluation scenario cell",
            )
            instances[row["instance_id"]] = metadata
            specifications[metadata[0]].add(metadata[1:])
        key = (row["instance_id"], row["strategy"], row["task"])
        check(
            key not in identities and "answer_json" in row,
            "Duplicate task or missing exact answer",
        )
        identities.add(key)
        observed_panels.add(row["source_panel_id"])
        panel = descriptors[row["source_panel_id"]]
        strata = (
            ("overall", "all"),
            ("scenario", row["scenario"]),
            ("size", panel["size_stratum"]),
            ("leader", panel["leader_stratum"]),
        )
        for (dimension, stratum) in (strata):
            group_panels[(dimension, stratum)].add(row["source_panel_id"])
        if (row["strategy"] in ("fixed_and_native", "support_explicit")):
            pair_key = (row["instance_id"], row["task"])
            pair = paired.setdefault(pair_key, {})
            check(row["strategy"] not in pair, "Duplicate paired strategy")
            pair[row["strategy"]] = row
            if (len(pair) == 2):
                combined, full = pair["fixed_and_native"], pair["support_explicit"]
                check(
                    all(
                        combined[field] == full[field]
                        for (field) in (
                            (
                                "source_panel_id",
                                "oracle_eligible",
                                "expected_support_denominator",
                                "scenario",
                                "task",
                            )
                        )
                    ),
                    "Paired benchmark identities differ",
                )
                index = panel_indices[combined["source_panel_id"]]
                common = (
                    combined["oracle_eligible"]
                    and combined["answer_state"] == "ANSWERED"
                    and full["answer_state"] == "ANSWERED"
                )
                equal = (
                    common
                    and combined["answer_json"] == full["answer_json"]
                    and combined["reported_support_denominator"]
                    == full["reported_support_denominator"]
                )
                contrast_counts[row["task"]][index] += (
                    int(equal),
                    int(common),
                    int(combined["oracle_eligible"]),
                    int(
                        combined["oracle_eligible"]
                        and combined["answer_state"] == "ANSWERED"
                    ),
                    int(full["oracle_eligible"] and full["answer_state"] == "ANSWERED"),
                )
                del paired[pair_key]
        for (metric) in (METRICS):
            numerator, denominator = opportunity(row, metric)
            for (dimension, stratum) in (strata):
                key = (dimension, stratum, row["strategy"], row["task"], metric)
                if (key not in groups):
                    groups[key] = np.zeros((len(panel_ids), 2), dtype = np.int64)
                groups[key][panel_indices[row["source_panel_id"]]] += (
                    numerator,
                    denominator,
                )
    check(
        not paired and observed_panels == set(panel_ids),
        "Unpaired strategies or incomplete evaluation population",
    )
    expected = set(scenario_specs(EVALUATION_SEEDS))
    check(
        all(specifications[panel] == expected for (panel) in (panel_ids)),
        "Complete evaluation scenario and seed design is required",
    )
    check(
        len(identities) == len(instances) * len(STRATEGIES) * len(task_names),
        "Incomplete strategy by task design",
    )
    return panel_ids, groups, group_panels, contrast_counts, instances

def comparative_records(panel_ids, groups, group_panels, bootstrap):
    records = []
    for (key, counts) in (sorted(groups.items())):
        dimension, stratum, strategy, task, metric = key
        numerator, denominator = (int(value) for (value) in (counts.sum(axis = 0)))
        record = {
            "resource": "OpenTargets_controlled",
            "baseline_release": "26.06",
            "followup_release": "controlled",
            "unit": "source_panel_cluster_and_task_instance",
            "dimension": dimension,
            "stratum": stratum,
            "strategy": strategy,
            "task": task,
            "metric": metric,
            **closed_fraction(numerator, denominator),
            "source_panel_count": len(group_panels[(dimension, stratum)]),
            "baseline_sample_panel_count": len(panel_ids),
            "bootstrap_lower": None,
            "bootstrap_upper": None,
            "bootstrap_evaluable_resamples": None,
        }
        if (dimension == "overall" and denominator > 0):
            replicates = bootstrap @ counts
            valid = replicates[:, 1] > 0
            values = replicates[valid, 0] / replicates[valid, 1]
            record["bootstrap_lower"], record["bootstrap_upper"] = (
                float(value)
                for (value) in (np.quantile(values, [0.025, 0.975], method = "linear"))
            )
            record["bootstrap_evaluable_resamples"] = int(valid.sum())
        records.append(record)
    return records

def paired_contrasts(panel_ids, contrast_counts, bootstrap):
    contrasts = []
    for (task, counts) in (sorted(contrast_counts.items())):
        equal, common, eligible, combined_covered, full_covered = (
            int(value) for (value) in (counts.sum(axis = 0))
        )
        records_bootstrap = bootstrap @ counts
        valid = records_bootstrap[:, 2] > 0
        differences = (
            records_bootstrap[valid, 4] - records_bootstrap[valid, 3]
        ) / records_bootstrap[valid, 2]
        contrasts.append(
            {
                "task": task,
                "first_strategy": "fixed_and_native",
                "second_strategy": "support_explicit",
                "common_answer_agreement": closed_fraction(equal, common),
                "common_answer_coverage": closed_fraction(common, eligible),
                "coverage_difference_full_minus_combined": (
                    full_covered - combined_covered
                )
                / eligible
                if (eligible)
                else None,
                "coverage_difference_bootstrap_lower": float(
                    np.quantile(differences, 0.025)
                )
                if (len(differences))
                else None,
                "coverage_difference_bootstrap_upper": float(
                    np.quantile(differences, 0.975)
                )
                if (len(differences))
                else None,
                "source_panel_count": len(panel_ids),
                "common_task_instances": common,
            }
        )
    return {
        "contrasts": contrasts,
        "bootstrap_resamples": 2000,
        "bootstrap_seed": 960625,
        "cluster": "source_panel",
        "interpretation": "equivalent_answers_are_retained; no_global_winner_score",
    }

def endpoint_rows(record, descriptor):
    row = record["comparison"]
    common = {
        "instance_id": record["instance_id"],
        "source_panel_id": row["panel_id"],
        "scenario": record["scenario"],
        "level": record["level"],
        "seed": record["seed"],
        "size_stratum": descriptor["size_stratum"],
        "leader_stratum": descriptor["leader_stratum"],
        "baseline_leader_margin": descriptor["native_leader_margin"],
    }
    for (view) in (("native", "fixed")):
        for (label, endpoint) in ([("leader", row[view]["leader"])] + [
            ("top_" + key, value) for ((key, value)) in (row[view]["top_k"].items())
        ]):
            yield {
                **common,
                "view": view,
                "metric": label + "_set_changed",
                "value": float(endpoint["value"]) if (endpoint["estimable"]) else None,
                "reason": endpoint["reason"],
                "support_denominator": endpoint["support_denominator"],
                "persistence": None,
                "tie_variant": None,
            }
        rbo = row[view]["rbo"]
        for (variant) in (("w", "b")):
            for (persistence) in (("0.9", "0.8", "0.95")):
                for (measure) in (("ext", "min", "max", "res")):
                    value = (
                        rbo["value"][variant][persistence][measure]
                        if (rbo["estimable"])
                        else None
                    )
                    yield {
                        **common,
                        "view": view,
                        "metric": "rank_biased_overlap_" + measure,
                        "value": value,
                        "reason": rbo["reason"],
                        "support_denominator": rbo["support_denominator"],
                        "persistence": float(persistence),
                        "tie_variant": variant,
                    }
    for (metric) in (("kendall_tau_b", "spearman")):
        value = row["fixed"]["rank"]
        yield {
            **common,
            "view": "fixed",
            "metric": metric,
            "value": value.get(metric),
            "reason": value["reason"],
            "support_denominator": value["support_denominator"],
            "persistence": None,
            "tie_variant": None,
        }
    for (metric) in ((
        "mean_signed_change",
        "mean_absolute_change",
        "median_absolute_change",
        "q90_absolute_change",
        "maximum_absolute_change",
    )):
        value = row["fixed"]["score"]
        yield {
            **common,
            "view": "fixed",
            "metric": metric,
            "value": value.get(metric),
            "reason": value["reason"],
            "support_denominator": value["support_denominator"],
            "persistence": None,
            "tie_variant": None,
        }

def write_endpoints(source, output, population, instances):
    descriptors = {row["disease_id"]: row for (row) in (population)}
    fields = [
        (name, pa.string())
        for (name) in (("instance_id", "source_panel_id", "scenario"))
    ]
    fields += [("level", pa.float64()), ("seed", pa.int64())]
    fields += [(name, pa.string()) for (name) in (("size_stratum", "leader_stratum"))]
    fields += [("baseline_leader_margin", pa.float64())]
    fields += [(name, pa.string()) for (name) in (("view", "metric"))]
    fields += [
        ("value", pa.float64()),
        ("reason", pa.string()),
        ("support_denominator", pa.int64()),
        ("persistence", pa.float64()),
        ("tie_variant", pa.string()),
    ]
    schema = pa.schema(fields)
    pending = []
    observed = set()
    with pq.ParquetWriter(
        output / "endpoint_values.parquet", schema, compression = "zstd"
    ) as writer:
        for (record) in (read_rows(source)):
            check(
                record["phase"] in ("development", "evaluation"),
                "Unknown benchmark phase",
            )
            if (record["phase"] != "evaluation"):
                continue
            panel = record["comparison"]["panel_id"]
            check(
                panel in descriptors and record["instance_id"] not in observed,
                "Unexpected or repeated endpoint instance",
            )
            check(
                instances.get(record["instance_id"])
                == (panel, record["scenario"], record["level"], record["seed"]),
                "Endpoint and task identities differ",
            )
            observed.add(record["instance_id"])
            pending.extend(endpoint_rows(record, descriptors[panel]))
            if (len(pending) >= 10000):
                writer.write_table(pa.Table.from_pylist(pending, schema = schema))
                pending.clear()
        if (pending):
            writer.write_table(pa.Table.from_pylist(pending, schema = schema))
    check(observed == set(instances), "Incomplete evaluation endpoints")

def summarise_benchmark(config_path):
    config = read_json(config_path)
    root = Path(config_path).resolve().parent
    population = list(read_rows(root / config["panels"]))
    panel_ids, groups, group_panels, contrast_counts, instances = aggregate_tasks(
        root / config["tasks"], population
    )
    bootstrap = np.random.default_rng(960625).multinomial(
        len(panel_ids), np.full(len(panel_ids), 1 / len(panel_ids)), size = 2000
    )
    output = (root / config["output"]).resolve()
    output.mkdir(parents = True, exist_ok = False)
    records = comparative_records(panel_ids, groups, group_panels, bootstrap)
    pq.write_table(
        pa.Table.from_pylist(records),
        output / "comparative_summary.parquet",
        compression = "zstd",
    )
    write_json(
        output / "strong_baseline_contrasts.json",
        paired_contrasts(panel_ids, contrast_counts, bootstrap),
    )
    write_endpoints(root / config["comparisons"], output, population, instances)
