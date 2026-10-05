import gzip
import json
from dataclasses import asdict
from pathlib import Path

from . import benchmark_run as benchmark
from . import benchmark_scenarios as scenarios
from .comparators import compare_panels
from .contracts import MappingProposal, PanelSnapshot
from .io import check, read_json, read_rows, write_json, write_parquet
from .personadrive_comparison import compare_cohort

def resource_paths(config_path):
    config = read_json(config_path)
    root = Path(config_path).resolve().parent
    return config, root

def full_path(root, value):
    path = Path(value)
    return path if (path.is_absolute()) else root / path

def panel(record):
    return PanelSnapshot(
        record["resource"],
        record["release"],
        record["panel_id"],
        tuple(tuple(row) for (row) in (record["members"])),
        record.get("identity_decidable", True),
    )

def compare_resources(config_path):
    config, root = resource_paths(config_path)
    pairs = read_json(full_path(root, config["pairs"]))
    output = full_path(root, config["output"])
    output.parent.mkdir(parents = True, exist_ok = True)
    with output.open("w", encoding = "utf-8") as handle:
        for (record) in (pairs):
            mappings = [
                MappingProposal(**row) for (row) in (record.get("mappings", []))
            ]
            result = compare_panels(
                panel(record["baseline"]),
                panel(record["followup"]),
                mapping_proposals = mappings,
                **config.get("comparison", {}),
            )
            handle.write(json.dumps(result, allow_nan = False) + "\n")

def compare_personadrive(config_path):
    config, root = resource_paths(config_path)
    baseline = read_json(full_path(root, config["baseline"]))
    followup = read_json(full_path(root, config["followup"]))
    result = compare_cohort(baseline, followup, config["samples"], config["settings"])
    write_json(full_path(root, config["output"]), result)

def controlled_benchmark(config_path):
    config, root = resource_paths(config_path)
    baselines = read_json(full_path(root, config["baselines"]))
    design = list(read_rows(full_path(root, config["design"])))
    check(
        len({row["instance_id"] for (row) in (design)}) == len(design),
        "Duplicate scenario identity",
    )
    check(
        {int(row["seed"]) for (row) in (design)}
        == set(scenarios.DEVELOPMENT_SEEDS + scenarios.EVALUATION_SEEDS),
        "All five fixed scenario seeds are required",
    )
    check(
        {row["scenario"] for (row) in (design)} == set(scenarios.SCENARIO_CLASSES),
        "All ten scenario classes are required",
    )
    expected = set(
        scenarios.scenario_specs(
            scenarios.DEVELOPMENT_SEEDS + scenarios.EVALUATION_SEEDS
        )
    )
    observed = {}
    for (row) in (design):
        level = None if (row["level"] in (None, "", "null")) else float(row["level"])
        observed.setdefault(row["panel_id"], set()).add(
            (row["scenario"], level, int(row["seed"]))
        )
    check(
        set(observed) == set(baselines)
        and all(value == expected for (value) in (observed.values())),
        "Complete panel by scenario by seed design is required",
    )
    check(
        len(design) == len(baselines) * len(expected),
        "Repeated panel by scenario by seed design cell",
    )
    output = full_path(root, config["output"])
    output.mkdir(parents = True, exist_ok = False)
    tasks = []
    with (
        gzip.open(
            output / "scenarios.jsonl.gz", "wt", encoding = "utf-8"
        ) as source_handle,
        gzip.open(
            output / "comparisons.jsonl.gz", "wt", encoding = "utf-8"
        ) as result_handle,
    ):
        for (row) in (design):
            seed = int(row["seed"])
            level = (
                None if (row["level"] in (None, "", "null")) else float(row["level"])
            )
            scenario = scenarios.create_scenario(
                row["panel_id"],
                baselines[row["panel_id"]],
                row["scenario"],
                level,
                seed,
                int(row["random_seed"]),
                row["identifier_prefix"],
            )
            left, right, mappings = benchmark.scenario_inputs(scenario)
            report = compare_panels(left, right, mapping_proposals = mappings)
            truth = benchmark.expected_answers(scenario)
            benchmark.evaluate_truth(report, truth)
            phase = (
                "development" if (seed in scenarios.DEVELOPMENT_SEEDS) else "evaluation"
            )
            tasks.extend(
                record | {"level": level, "seed": seed}
                for (record) in (
                    benchmark.task_records(
                        row["instance_id"],
                        row["panel_id"],
                        row["scenario"],
                        phase,
                        report,
                        truth,
                    )
                )
            )
            source_handle.write(
                json.dumps(
                    {"instance_id": row["instance_id"], **asdict(scenario)},
                    allow_nan = False,
                )
                + "\n"
            )
            result_handle.write(
                json.dumps(
                    {
                        "instance_id": row["instance_id"],
                        "phase": phase,
                        "scenario": row["scenario"],
                        "level": level,
                        "seed": seed,
                        "comparison": report,
                        "truth": truth,
                    },
                    allow_nan = False,
                )
                + "\n"
            )
    write_parquet(output / "task_answers.parquet", tasks)
