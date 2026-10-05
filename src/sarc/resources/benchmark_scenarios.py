import math
import random
import statistics
from dataclasses import dataclass

SCENARIO_CLASSES = (
    "unchanged",
    "row_permutation",
    "monotone_score_rescaling",
    "persistent_score_perturbation",
    "persistent_leader_replacement",
    "roster_entry",
    "roster_exit",
    "exact_reversible_identifier_rename",
    "ambiguous_mapping",
    "mixed_score_roster_mapping",
)
LEVEL_FREE = frozenset(
    {
        "unchanged",
        "row_permutation",
        "monotone_score_rescaling",
        "persistent_leader_replacement",
    }
)
LEVELS = (0.01, 0.05, 0.2)
DEVELOPMENT_SEEDS = (1103, 2207)
EVALUATION_SEEDS = (3301, 4409, 5501)

@dataclass(frozen = True)
class Scenario:
    panel_id: str
    scenario: str
    level: float | None
    seed: int
    baseline: tuple[tuple[str, float], ...]
    followup: tuple[tuple[str, str, float], ...]
    mappings: tuple[tuple[str, str], ...]
    ledger: tuple[dict, ...]
    identity_decidable: bool

def deterministic_generator(random_seed):
    if (not isinstance(random_seed, int)):
        raise ValueError("The archived integer random seed is required")
    return random.Random(random_seed)

def selected_members(values, count, generator):
    return generator.sample(sorted(values), min(count, len(values)))

def create_scenario(
    panel_id, scores, scenario, level, seed, random_seed, identifier_prefix
):
    if (
        scenario not in SCENARIO_CLASSES
        or (scenario in LEVEL_FREE and level is not None)
        or (scenario not in LEVEL_FREE and level not in LEVELS)
    ):
        raise ValueError(
            "scenario and intervention level do not match the frozen registry"
        )
    if (not scores or any(
        (
            not math.isfinite(value) or not 0 < value <= 1
            for (value) in (scores.values())
        )
    )):
        raise ValueError(
            "scenario baselines require observed positive Open Targets scores"
        )
    generator = deterministic_generator(random_seed)
    baseline = tuple(sorted(scores.items()))
    values = dict(baseline)
    mappings = []
    ledger = []
    amount = max(1, math.ceil((level or 0) * len(baseline)))
    prefix = identifier_prefix
    if (not isinstance(prefix, str) or not prefix):
        raise ValueError("The archived synthetic-identifier prefix is required")
    if (scenario == "monotone_score_rescaling"):
        values = {member: value * 0.5 for ((member, value)) in (values.items())}
        ledger.append({"operation": "multiply_all_scores", "factor": 0.5})
    if (scenario in {"persistent_score_perturbation", "mixed_score_roster_mapping"}):
        chosen = selected_members(values, amount, generator)
        for (index, member) in (enumerate(chosen)):
            sign = 1 if (index % 2 == 0) else -1
            values[member] = min(
                1.0, max(math.nextafter(0.0, 1.0), values[member] * (1 + sign * level))
            )
        ledger.append(
            {
                "operation": "score_perturbation",
                "members": chosen,
                "balanced_sign_order": "positive_then_negative",
                "fraction": level,
            }
        )
    if (scenario == "persistent_leader_replacement"):
        maximum = max(values.values())
        leaders = sorted(
            (member for ((member, value)) in (values.items()) if (value == maximum))
        )
        nonleaders = sorted(
            (member for ((member, value)) in (values.items()) if (value < maximum))
        )
        if (nonleaders):
            (first, second) = (leaders[0], nonleaders[0])
            (values[first], values[second]) = (values[second], values[first])
            ledger.append(
                {"operation": "leader_score_swap", "members": [first, second]}
            )
        else:
            values = {
                member: value if (member == leaders[0]) else value * 0.5
                for ((member, value)) in (values.items())
            }
            ledger.append(
                {"operation": "all_tied_keep_one", "retained_leader": leaders[0]}
            )
    if (scenario in {"roster_entry", "mixed_score_roster_mapping"}):
        (maximum, median) = (max(values.values()), statistics.median(values.values()))
        added = []
        for (index) in (range(amount)):
            member = f"entry_{prefix}_{index:06d}"
            values[member] = maximum if (index < math.ceil(amount / 2)) else median
            added.append(member)
        ledger.append(
            {
                "operation": "entry",
                "members": added,
                "maximum_count": math.ceil(amount / 2),
                "maximum_score": maximum,
                "median_score": median,
            }
        )
    if (scenario in {"roster_exit", "mixed_score_roster_mapping"}):
        removed = selected_members(values, amount, generator)
        for (member) in (removed):
            del values[member]
        ledger.append({"operation": "exit", "members": removed})
    if (scenario in {"exact_reversible_identifier_rename", "mixed_score_roster_mapping"}):
        candidates = {
            member: value for ((member, value)) in (values.items()) if (member in scores)
        }
        renamed = selected_members(candidates, amount, generator)
        for (index, member) in (enumerate(renamed)):
            replacement = f"rename_{prefix}_{index:06d}"
            values[replacement] = values.pop(member)
            mappings.append((member, replacement))
        ledger.append({"operation": "reversible_rename", "mappings": list(mappings)})
    identity_decidable = True
    rows = [(member, member, value) for ((member, value)) in (values.items())]
    if (scenario == "ambiguous_mapping"):
        count = min(len(values) - len(values) % 2, max(2, amount + amount % 2))
        chosen = selected_members(values, count, generator)
        replacements = {}
        for (index, member) in (enumerate(chosen)):
            replacement = f"merged_{prefix}_{index // 2:06d}"
            replacements[member] = replacement
            mappings.append((member, replacement))
        rows = [
            (member, replacements.get(member, member), value)
            for ((member, value)) in (values.items())
        ]
        ledger.append(
            {
                "operation": "ambiguous_pairs",
                "mappings": list(mappings),
                "raw_rows_preserved": True,
            }
        )
        identity_decidable = False
    if (scenario == "row_permutation"):
        generator.shuffle(rows)
        ledger.append(
            {"operation": "row_permutation", "row_order": [row[0] for (row) in (rows)]}
        )
    if (scenario == "unchanged"):
        ledger.append({"operation": "unchanged"})
    return Scenario(
        panel_id,
        scenario,
        level,
        seed,
        baseline,
        tuple(rows),
        tuple(mappings),
        tuple(ledger),
        identity_decidable,
    )

def scenario_specs(seeds):
    for (scenario) in (SCENARIO_CLASSES):
        for (level) in ((None,) if (scenario in LEVEL_FREE) else LEVELS):
            for (seed) in (seeds):
                yield (scenario, level, seed)
