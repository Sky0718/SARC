from collections import Counter, defaultdict

from ..analysis.run import run_operation
from .comparison import compare, project
from .cost_join import COUNTS, build, require

def describe(catalogue):
    require(
        len(catalogue["rows"]) == catalogue["total_rows"] == 32268,
        "The complete integration catalogue is required",
    )
    rows = [project(row) for (row) in (catalogue["rows"])]
    require(
        len({row["catalog_id"] for (row) in (rows)}) == len(rows),
        "Duplicate catalogue identity",
    )
    counts = dict(Counter(row["family"] for (row) in (rows)))
    require(counts == COUNTS, "Complete branch and policy families are required")
    panels = defaultdict(list)
    for (row) in (rows):
        panels[row["panel_id"]].append(row)
    require(len(panels) == 10338, "Scientific comparison contracts differ")
    for (group) in (panels.values()):
        compare(group)
    return {
        "rows": rows,
        "counts": counts,
        "panels": len(panels),
        "classification_counts": dict(
            Counter(row["descriptive_point"]["status"] for (row) in (rows))
        ),
        "cost_frontier_claimed": False,
    }

def join_costs(descriptive, non_b_lineage, b_packets, source_ids):
    require(
        len(non_b_lineage["rows"]) == 32088 and len(non_b_lineage["scenarios"]) == 379,
        "Complete physical-work metadata are required",
    )
    result = build(descriptive, non_b_lineage, b_packets, source_ids)
    require(
        result["counts"] == COUNTS and result["panels"] == 10338,
        "Incomplete final integration",
    )
    return result

def run(config_path):
    return run_operation(config_path, {"describe": describe, "join_costs": join_costs})
