from pathlib import Path

from ..resources.io import check, read_json, write_json
from . import background, design, residuals

def inputs(config_path):
    config = read_json(config_path)
    root = Path(config_path).resolve().parent
    records = {
        name: read_json(root / value) for ((name, value)) in (config["inputs"].items())
    }
    return config, root, records

def assess_design(config_path):
    config, root, records = inputs(config_path)
    rows = records["rows"]
    expected = {
        (cohort, transition, direction, label, seed)
        for (cohort) in (design.COHORTS)
        for (transition) in (design.TRANSITIONS)
        for (direction) in (design.DIRECTIONS)
        for (label) in (design.CLASSES)
        for (seed) in (design.SEEDS)
    }
    identities = [
        (
            row["cohort"],
            row["transition"],
            row["direction"],
            row["class_id"],
            row["control_seed"],
        )
        for (row) in (rows)
    ]
    check(
        len(identities) == len(set(identities)) and set(identities) == expected,
        "Complete matched-control design is required",
    )
    rows.sort(
        key = lambda row: (
            list(design.COHORTS).index(row["cohort"]),
            design.TRANSITIONS.index(row["transition"]),
            design.DIRECTIONS.index(row["direction"]),
            design.CLASSES.index(row["class_id"]),
            design.SEEDS.index(row["control_seed"]),
        )
    )
    prepared = []
    for (row) in (rows):
        prepared.append(
            row
            | {
                "stratum": "::".join(
                    (
                        row["cohort"],
                        row["transition"],
                        row["direction"],
                        row["class_id"],
                    )
                ),
                "covariates": design.residuals(row["balance"], row.get("components")),
            }
        )
    result = {
        "rows": prepared,
        "cohort_designs": {
            cohort: design.design(
                [row for (row) in (prepared) if (row["cohort"] == cohort)]
            )
            for (cohort) in (design.COHORTS)
        },
    }
    write_json(root / config["output"], result)

def project_residuals(config_path):
    config, root, records = inputs(config_path)
    write_json(
        root / config["output"],
        residuals.project_all(records["assessment"], records["scalars"]),
    )

def background_contrasts(config_path):
    config, root, records = inputs(config_path)
    output = []
    for (group) in (records["groups"]):
        cohort = group["cohort"]
        check(cohort in background.COHORTS, "Unknown cohort")
        check(
            group["transition"] in background.TRANSITIONS
            and group["class_id"] in background.CLASSES,
            "Unknown transition or edit class",
        )
        check(
            len(group["insert"]) == len(group["retract"]) == background.COHORTS[cohort],
            "Full cohort is required",
        )
        for (direction) in (background.DIRECTIONS):
            check(
                all(
                    row["cohort"] == cohort
                    and row["transition"] == group["transition"]
                    and row["class_id"] == group["class_id"]
                    and row["direction"] == direction
                    for (row) in (group[direction])
                ),
                "Background group and row identities differ",
            )
        old = {row["sample_id"]: row for (row) in (group["insert"])}
        new = {row["sample_id"]: row for (row) in (group["retract"])}
        check(
            len(old) == background.COHORTS[cohort] and set(old) == set(new),
            "Paired sample roster differs",
        )
        for (budget) in (background.BUDGETS):
            paired = [
                background.paired_case(old[sample], new[sample], budget)
                for (sample) in (old)
            ]
            output.append(
                {
                    "cohort": cohort,
                    "transition": group["transition"],
                    "class_id": group["class_id"],
                    "k": budget,
                    "rows": paired,
                    "summary": background.summarise(paired),
                }
            )
    write_json(root / config["output"], output)
