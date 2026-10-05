import itertools
import math
import struct
from fractions import Fraction
from . import representation_design as mechanism

SEEDS = mechanism.SEEDS
REPRESENTATIONS = mechanism.REPRESENTATIONS
REPEATS = mechanism.REPEATS
LABELS = tuple(itertools.product(REPRESENTATIONS, REPEATS))
KS = (1, 5, 10, 20)
SUCCESSES = {"SUCCESS", "SUCCESS_EMPTY"}
STATUSES = SUCCESSES | {"FAILED", "UNAVAILABLE"}

def _finite(value):
    return (
        isinstance(value, (float, int))
        and (not isinstance(value, bool))
        and math.isfinite(value)
    )

def checked_returned(observation):
    status = observation["status"]
    if (status not in STATUSES):
        raise ValueError("Unexpected observation status")
    if (status not in SUCCESSES):
        return None
    rows = observation["returned"]
    if (not isinstance(rows, list)):
        raise ValueError("Native returned list is required")
    genes, scores = ([], [])
    for (row) in (rows):
        gene, score = (row["gene"], row["score"])
        if (not isinstance(gene, str) or not gene or gene.strip() != gene):
            raise ValueError("Invalid native gene identity")
        if (not _finite(score)):
            raise ValueError("Native score must be finite binary64")
        raw = bytes.fromhex(row["float64_le_hex"])
        if (len(raw) != 8):
            raise ValueError("Expected exact binary64 score carrier")
        decoded = struct.unpack("<d", raw)[0]
        if (not math.isfinite(decoded) or struct.pack("<d", float(score)) != raw):
            raise ValueError("JSON score differs from lossless native score")
        if ("decimal17" in row and float(row["decimal17"]) != decoded):
            raise ValueError("Decimal carrier disagrees with native score")
        genes.append(gene)
        scores.append(decoded)
    if (len(set(genes)) != len(genes)):
        raise ValueError("Duplicate native gene identity")
    if (any((left < right for ((left, right)) in (zip(scores, scores[1:]))))):
        raise ValueError("Native scores are not descending; no silent resort")
    if ((status == "SUCCESS_EMPTY") != (not rows)):
        raise ValueError("Successful-empty status and returned rows disagree")
    if ("genes" in observation and observation["genes"] != genes):
        raise ValueError("Redundant native gene carrier differs")
    if ("scores" in observation and observation["scores"] != scores):
        raise ValueError("Redundant native score carrier differs")
    return list(zip(genes, scores))

def fractional_membership(rows, k):
    if (type(k) is not int or k < 1):
        raise ValueError("Positive integer k required")
    membership = {}
    start = 0
    remaining = k
    while (start < len(rows)):
        end = start + 1
        while (end < len(rows) and rows[end][1] == rows[start][1]):
            end += 1
        weight = Fraction(min(max(remaining, 0), end - start), end - start)
        for (gene, _) in (rows[start:end]):
            membership[gene] = weight
        remaining -= end - start
        start = end
    assert (sum(membership.values(), Fraction()) == min(k, len(rows)))
    return membership

def exact_l1(left, right, k):
    a, b = (fractional_membership(left, k), fractional_membership(right, k))
    distance = sum(
        (
            abs(a.get(gene, Fraction()) - b.get(gene, Fraction()))
            for (gene) in (a.keys() | b.keys())
        ),
        Fraction(),
    )
    return (distance, sum(a.values(), Fraction()), sum(b.values(), Fraction()))

def rational(value):
    return f"{value.numerator}/{value.denominator}"

def converted(value, exact = False):
    if (isinstance(value, Fraction)):
        return rational(value) if (exact) else float(value)
    if (isinstance(value, dict)):
        return {key: converted(item, exact) for ((key, item)) in (value.items())}
    if (isinstance(value, (tuple, list))):
        return [converted(item, exact) for (item) in (value)]
    return value

def pair_metrics(left, right):
    if ((left["case_id"], left["master_seed"]) != (
        right["case_id"],
        right["master_seed"],
    )):
        raise ValueError("Only same case/master observations may be paired")
    a, b = (checked_returned(left), checked_returned(right))
    left_rep = left.get("physical_representation", left["representation"])
    right_rep = right.get("physical_representation", right["representation"])
    result = {
        "case_id": left["case_id"],
        "master_seed": left["master_seed"],
        "left": [left["representation"], left["repeat"]],
        "right": [right["representation"], right["repeat"]],
        "left_logical_task_id": left.get("logical_task_id", left.get("task_id")),
        "right_logical_task_id": right.get("logical_task_id", right.get("task_id")),
        "left_status": left["status"],
        "right_status": right["status"],
        "status": "AVAILABLE" if (a is not None and b is not None) else "UNAVAILABLE",
        "alias_pair": left_rep == right_rep and left["repeat"] == right["repeat"],
        "representation_alias": left_rep == right_rep
        and left["representation"] != right["representation"],
        "top10_fractional_L1": None,
        "topk": {},
    }
    secondary = (
        "order_equal",
        "positional_differences",
        "shared_positions",
        "unmatched_tail_positions",
        "support_added",
        "support_removed",
        "support_jaccard",
        "score_equal",
        "max_abs",
        "max_symmetric_relative",
        "shared_support_count",
    )
    result.update({field: None for (field) in (secondary)})
    if (a is None or b is None):
        result["unavailable_reasons"] = [
            {"side": side, "status": obs["status"], "reason": obs.get("reason")}
            for ((side, obs)) in ((("left", left), ("right", right)))
            if (obs["status"] not in SUCCESSES)
        ]
        for (k) in (KS):
            result["topk"][str(k)] = {
                key: None
                for (key) in (
                    (
                        "fractional_L1",
                        "fractional_L1_exact",
                        "left_membership_sum",
                        "right_membership_sum",
                    )
                )
            }
        return (result, None)
    exact_primary = None
    for (k) in (KS):
        distance, mass_a, mass_b = exact_l1(a, b, k)
        result["topk"][str(k)] = {
            "fractional_L1": float(distance),
            "fractional_L1_exact": rational(distance),
            "left_membership_sum": float(mass_a),
            "right_membership_sum": float(mass_b),
        }
        if (k == 10):
            exact_primary = distance
            result["top10_fractional_L1"] = float(distance)
    genes_a, genes_b = ([gene for ((gene, _)) in (a)], [gene for ((gene, _)) in (b)])
    map_a, map_b = (dict(a), dict(b))
    support_a, support_b = (set(map_a), set(map_b))
    shared = support_a & support_b
    union = support_a | support_b
    absolute = [abs(map_a[gene] - map_b[gene]) for (gene) in (shared)]
    relative = []
    for (gene) in (shared):
        x, y = (map_a[gene], map_b[gene])
        scale = max(abs(x), abs(y))
        relative.append(0.0 if (scale == 0) else abs(x / scale - y / scale))
    if (any((not math.isfinite(value) for (value) in (absolute + relative)))):
        raise ValueError("Score difference overflows reporting precision")
    result.update(
        {
            "order_equal": genes_a == genes_b,
            "positional_differences": sum(
                (x != y for ((x, y)) in (zip(genes_a, genes_b)))
            )
            + abs(len(a) - len(b)),
            "shared_positions": min(len(a), len(b)),
            "unmatched_tail_positions": abs(len(a) - len(b)),
            "support_added": [gene for (gene) in (genes_b) if (gene not in support_a)],
            "support_removed": [
                gene for (gene) in (genes_a) if (gene not in support_b)
            ],
            "support_jaccard": len(shared) / len(union) if (union) else None,
            "score_equal": support_a == support_b
            and all((map_a[gene] == map_b[gene] for (gene) in (shared))),
            "max_abs": max(absolute) if (shared) else None,
            "max_symmetric_relative": max(relative) if (shared) else None,
            "shared_support_count": len(shared),
        }
    )
    return (result, exact_primary)

def five_seed_summary(values):
    if (len(values) != 5):
        raise ValueError("Five planned master values required")
    count = sum((value is not None for (value) in (values)))
    return {
        "complete_five_seed_mean": sum(values, Fraction()) / 5
        if (count == 5)
        else None,
        "available_seed_count": count,
        "all_seed_values": values,
    }

def analyse(observation_document, case_document):
    if (observation_document.get("schema") != "item3_normalised_lossless_v1"):
        raise ValueError("Wrong lossless input schema")
    cases = case_document["cases"]
    case_ids = [case["case_id"] for (case) in (cases)]
    if (len(cases) != 8 or len(set(case_ids)) != 8):
        raise ValueError("Exactly eight frozen case-graph units required")
    observations = observation_document["observations"]
    if (len(observations) != 240):
        raise ValueError("All 240 logical observations, including failures, required")
    by_key = {}
    case_map = {case["case_id"]: case for (case) in (cases)}
    for (observation) in (observations):
        key = tuple(
            (
                observation[field]
                for (field) in (("case_id", "master_seed", "representation", "repeat"))
            )
        )
        if (
            key[0] not in case_map
            or key[1] not in SEEDS
            or key[2] not in REPRESENTATIONS
            or (key[3] not in REPEATS)
        ):
            raise ValueError("Unexpected observation identity")
        if (key in by_key):
            raise ValueError("Duplicate logical observation")
        aliases = case_map[key[0]]["aliases"]
        if (set(aliases) != set(REPRESENTATIONS) or any(
            (aliases[aliases[name]] != aliases[name] for (name) in (REPRESENTATIONS))
        )):
            raise ValueError("Malformed frozen representation aliases")
        if (observation.get("physical_representation", key[2]) != aliases[key[2]]):
            raise ValueError("Physical representation differs from frozen alias map")
        checked_returned(observation)
        by_key[key] = observation
    if (len(by_key) != 8 * len(SEEDS) * len(LABELS)):
        raise ValueError("Incomplete planned identity grid")
    for (case_id, seed, representation, repeat) in (by_key):
        observation = by_key[case_id, seed, representation, repeat]
        physical = case_map[case_id]["aliases"][representation]
        reference = by_key[case_id, seed, physical, repeat]
        if ((observation["status"], observation.get("returned")) != (
            reference["status"],
            reference.get("returned"),
        )):
            raise ValueError(
                "Alias observations disagree with their physical observation"
            )
    pairs, blocks, case_results = ([], [], [])
    for (case) in (cases):
        case_blocks = []
        for (seed) in (SEEDS):
            block_pairs = []
            for (left, right) in (itertools.combinations(LABELS, 2)):
                pair, exact_primary = pair_metrics(
                    by_key[(case["case_id"], seed) + tuple(left)],
                    by_key[(case["case_id"], seed) + tuple(right)],
                )
                pairs.append(pair)
                block_pairs.append(
                    {"left": left, "right": right, "top10_fractional_L1": exact_primary}
                )
            exact_summary = mechanism.block_summary(block_pairs)
            available = sum(
                (pair["top10_fractional_L1"] is not None for (pair) in (block_pairs))
            )
            block = {
                "case_id": case["case_id"],
                "master_seed": seed,
                "summary": exact_summary,
                "planned_pair_count": 15,
                "available_pair_count": available,
                "all_pairs_available": available == 15,
            }
            case_blocks.append(block)
            blocks.append(
                {
                    **converted(block),
                    "exact_summary": converted(exact_summary, exact = True),
                }
            )
        exact_case = mechanism.case_summary(case_blocks)
        within = {
            representation: five_seed_summary(
                [
                    block["summary"]["within_representation_repeat_distances"][
                        representation
                    ]
                    for (block) in (case_blocks)
                ]
            )
            for (representation) in (REPRESENTATIONS)
        }
        cross = {
            representation: {
                metric: five_seed_summary(
                    [
                        block["summary"]["contrasts"][representation][metric]
                        for (block) in (case_blocks)
                    ]
                )
                for (metric) in (("cross_mean", "mean_within_repeat"))
            }
            for (representation) in (REPRESENTATIONS[1:])
        }
        case_results.append(
            {
                "case_id": case["case_id"],
                "metadata": case,
                "summary": converted(exact_case),
                "exact_summary": converted(exact_case, exact = True),
                "within_representation_repeat_distances": converted(within),
                "within_representation_repeat_distances_exact": converted(
                    within, exact = True
                ),
                "cross_representation": converted(cross),
                "cross_representation_exact": converted(cross, exact = True),
            }
        )
    status_counts = {
        status: sum((row["status"] == status for (row) in (observations)))
        for (status) in (sorted(STATUSES))
    }
    return {
        "pairs.json": {
            "schema": "item3_pair_metrics_v1",
            "primary_k": 10,
            "secondary_k": [1, 5, 20],
            "planned_pairs": 600,
            "pairs": pairs,
        },
        "blocks.json": {
            "schema": "item3_block_summaries_v1",
            "planned_blocks": 40,
            "blocks": blocks,
        },
        "case_summaries.json": {
            "schema": "item3_case_summaries_v1",
            "planned_cases": 8,
            "cases": case_results,
        },
        "SUMMARY.json": {
            "schema": "item3_descriptive_summary_v1",
            "logical_observations": 240,
            "status_counts": status_counts,
            "planned_pairs": 600,
            "available_pairs": sum(
                (pair["status"] == "AVAILABLE" for (pair) in (pairs))
            ),
            "planned_blocks": 40,
            "complete_blocks": sum(
                (block["all_pairs_available"] for (block) in (blocks))
            ),
            "planned_cases": 8,
            "biological_samples": 2,
            "primary_k": 10,
            "secondary_k": [1, 5, 20],
            "biological_inference": False,
            "all_seed_values_retained": True,
            "alias_pairs": sum((pair["alias_pair"] for (pair) in (pairs))),
            "symmetric_relative_definition": "abs(a-b)/max(abs(a),abs(b)); both zero gives zero; shared support only",
            "interpretation": "FINITE_DESCRIPTIVE_REPRESENTATION_VERSUS_REPEAT_NOT_CAUSAL_OR_BIOLOGICAL_SUPERIORITY",
        },
    }
