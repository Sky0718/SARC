import itertools
import math
from fractions import Fraction

METHODS = ("DawnRank", "PRODIGY", "PersonaDrive")
SEEDS = (104729, 130363, 155921, 196613, 228017)
POLICIES = ("STABLE_ONLY_3", "FIXED_V2", "RC")
COMPARATORS = ("LATEST_PREFIX_R", "RC_PREFIX_R", "THINNING_EXPECTATION")
METRICS = ("returned", "assessed", "unassessed", "supported", "vacancies")

def require(condition, message):
    if (not condition):
        raise ValueError(message)

def seeds(method):
    return SEEDS if (method == "PRODIGY") else (None,)

def mean(values):
    return (
        None
        if (not values or any((value is None for (value) in (values))))
        else sum(values, Fraction(0)) / len(values)
    )

def roster(split, production = True):
    rows = split["models"]
    models = [row["sample_ID"] for (row) in (rows)]
    donors = split["ordered_donors"]
    require(
        len(models) == len(set(models)) and len(donors) == len(set(donors)),
        "DUPLICATE_ROSTER",
    )
    require(
        set(donors) == {row["study_donor_id"] for (row) in (rows)},
        "DONOR_ROSTER_CHANGED",
    )
    require(
        all((row["split"] in ("held_out", "development") for (row) in (rows))),
        "SPLIT_CHANGED",
    )
    require(
        all(
            (
                len(
                    {
                        row["split"]
                        for (row) in (rows)
                        if (row["study_donor_id"] == donor)
                    }
                )
                == 1
                for (donor) in (donors)
            )
        ),
        "DONOR_CROSSES_SPLIT",
    )
    if (production):
        require((len(rows), len(donors)) == (85, 83), "FULL_ROSTER_REQUIRED")
        require(
            [
                (
                    sum((row["split"] == label for (row) in (rows))),
                    len(
                        {
                            row["study_donor_id"]
                            for (row) in (rows)
                            if (row["split"] == label)
                        }
                    ),
                )
                for (label) in (("held_out", "development"))
            ]
            == [(59, 58), (26, 25)],
            "FULL_SPLITS_REQUIRED",
        )
    return rows

def load_ledger(
    rows, split, rank_binding, eligibility_binding, checkpoint, production = True
):
    models = {row["sample_ID"]: row for (row) in (roster(split, production))}
    groups, labels, count = ({}, {sample: {} for (sample) in (models)}, 0)
    for (row) in (rows):
        count += 1
        if (count % 20000 == 0):
            checkpoint()
        if ((row["normal_reference_variant"], row["support"], row["offered_k"]) != (
            "all34",
            "native",
            10,
        ) or row["policy_id"] not in POLICIES):
            continue
        sample, method, seed, policy = (
            row["model_id"],
            row["method"],
            row["seed"],
            row["policy_id"],
        )
        require(
            sample in models and method in METHODS and (seed in seeds(method)),
            "SELECTED_IDENTITY_MISMATCH",
        )
        require(
            row["donor_id"] == models[sample]["study_donor_id"]
            and row["split"] == models[sample]["split"],
            "LEDGER_ROSTER_MISMATCH",
        )
        require(
            row["rank_binding"] == rank_binding
            and row["input_binding"] == eligibility_binding,
            "LEDGER_INPUT_CHAIN_CHANGED",
        )
        key = (sample, method, seed, policy)
        state = row["prediction_state"]
        require(
            state in ("unavailable", "available_empty", "available_nonempty"),
            "PREDICTION_STATE_INVALID",
        )
        current = groups.setdefault(
            key, {"state": state, "weights": {}, "sentinel_count": 0}
        )
        require(current["state"] == state, "MIXED_PREDICTION_STATE")
        gene, call = (row["gene_id"], row["binary_call"])
        if (state != "available_nonempty"):
            require(
                gene is None
                and call is None
                and (row["weight_numerator"] is None)
                and (row["weight_denominator"] is None)
                and (row["assay_state"] == "not_applicable"),
                "BAD_SENTINEL",
            )
            current["sentinel_count"] += 1
            require(current["sentinel_count"] == 1, "DUPLICATE_SENTINEL")
            continue
        require(
            isinstance(gene, str) and gene and (gene not in current["weights"]),
            "DUPLICATE_OR_EMPTY_GENE",
        )
        n, d = (row["weight_numerator"], row["weight_denominator"])
        require(type(n) is int and type(d) is int and (d > 0), "BAD_RATIONAL_WEIGHT")
        weight = Fraction(n, d)
        require(0 < weight <= 1, "WEIGHT_OUTSIDE_UNIT_INTERVAL")
        require(
            row["assay_state"] == "eligible_measured"
            and call in (0, 1)
            or (row["assay_state"] == "not_measured" and call is None),
            "ASSAY_STATE_INVALID",
        )
        call = None if (call is None) else int(call)
        require(
            gene not in labels[sample] or labels[sample][gene] == call,
            "SHARED_ASSAY_LABEL_CONTRADICTION",
        )
        labels[sample][gene] = call
        current["weights"][gene] = weight
    expected = {
        (sample, method, seed, policy)
        for (sample) in (models)
        for (method) in (METHODS)
        for (seed) in (seeds(method))
        for (policy) in (POLICIES)
    }
    require(set(groups) == expected, "FULL_SEED_POLICY_GRID_REQUIRED")
    for (group) in (groups.values()):
        require(
            sum(group["weights"].values(), Fraction(0)) <= 10, "OFFERED_BUDGET_CHANGED"
        )
        require(
            (group["state"] == "available_nonempty") == bool(group["weights"]),
            "WEIGHT_STATE_MISMATCH",
        )
    if (production):
        require(count == 667365 and len(groups) == 1785, "FULL_LEDGER_COUNTS_CHANGED")
    return (groups, labels, count)

def rank_index(rows, split, production = True):
    result = {}
    for (row) in (rows):
        key = tuple(row["key"])
        require(key not in result, "DUPLICATE_RANK_STATE")
        genes, scores = (row["genes"], row["scores"])
        if (genes is None):
            require(scores is None, "UNAVAILABLE_RANK_SCORES")
        else:
            require(
                scores is not None
                and len(genes) == len(scores)
                and (len(set(genes)) == len(genes)),
                "RANK_SCHEMA_INVALID",
            )
            require(
                all((isinstance(gene, str) and gene for (gene) in (genes))),
                "RANK_GENE_INVALID",
            )
            require(
                all(
                    (
                        type(score) in (int, float) and math.isfinite(score)
                        for (score) in (scores)
                    )
                )
                and all((a >= b for ((a, b)) in (zip(scores, scores[1:])))),
                "RANK_SCORE_INVALID",
            )
        result[key] = (genes, scores)
    expected = {
        (method, pool, "native_" + release, seed, model["sample_ID"])
        for (model) in (roster(split, production))
        for (method) in (METHODS)
        for (seed) in (seeds(method))
        for (pool) in (("all34", "distal16", "proximal18"))
        for (release) in (("11_0", "11_5", "12_0"))
    }
    require(set(result) == expected, "FULL_RANK_GRID_REQUIRED")
    return result

def prefix(weights, scores, target):
    require(isinstance(target, Fraction) and target >= 0, "INVALID_TARGET_MASS")
    if (weights is None):
        return (None, "COMPARATOR_UNAVAILABLE")
    if (sum(weights.values(), Fraction(0)) < target):
        return (None, "COMPARATOR_MASS_INSUFFICIENT")
    if (target == 0):
        return ({}, "AVAILABLE_ZERO_MASS")
    require(set(weights) <= set(scores), "RANK_MEMBER_MISSING")
    ordered = sorted(weights, key = lambda gene: (-scores[gene], gene))
    result, remaining = ({}, target)
    for (score, members) in (itertools.groupby(ordered, key = lambda gene: scores[gene])):
        block = list(members)
        require(
            len({weights[gene] for (gene) in (block)}) == 1,
            "INCONSISTENT_ACCEPTED_TIE_CAP",
        )
        capacity = sum((weights[gene] for (gene) in (block)), Fraction(0))
        taken = min(remaining, capacity)
        if (taken):
            result.update(((gene, taken / len(block)) for (gene) in (block)))
            remaining -= taken
        if (remaining == 0):
            break
    require(
        remaining == 0
        and sum(result.values(), Fraction(0)) == target
        and all((value <= weights[gene] for ((gene, value)) in (result.items()))),
        "PREFIX_MASS_OR_CAP_FAILURE",
    )
    return (result, "AVAILABLE")

def metrics(weights, labels):
    if (weights is None):
        return dict.fromkeys(METRICS)
    require(set(weights) <= set(labels), "ACCEPTED_ASSAY_MEMBERSHIP_MISSING")
    returned = sum(weights.values(), Fraction(0))
    assessed = sum(
        (weight for ((gene, weight)) in (weights.items()) if (labels[gene] is not None)),
        Fraction(0),
    )
    supported = sum(
        (
            weight * labels[gene]
            for ((gene, weight)) in (weights.items())
            if (labels[gene] is not None)
        ),
        Fraction(0),
    )
    return {
        "returned": returned,
        "assessed": assessed,
        "unassessed": returned - assessed,
        "supported": supported,
        "vacancies": 10 - returned,
    }

def average_weights(rows):
    if (any((row is None for (row) in (rows)))):
        return None
    return {
        gene: sum((row.get(gene, Fraction(0)) for (row) in (rows)), Fraction(0))
        / len(rows)
        for (gene) in (sorted(set().union(*(set(row) for (row) in (rows)))))
    }

def sharp(left, right, labels):
    if (left is None or right is None):
        return {
            "observed_difference": None,
            "lower": None,
            "upper": None,
            "coefficients": None,
        }
    known, lower_extra, upper_extra, witness = (
        Fraction(0),
        Fraction(0),
        Fraction(0),
        [],
    )
    for (gene) in (sorted(set(left) | set(right))):
        require(gene in labels, "SHARED_LABEL_MISSING")
        a, b, call = (
            left.get(gene, Fraction(0)),
            right.get(gene, Fraction(0)),
            labels[gene],
        )
        delta = a - b
        if (call is None):
            lower_extra += min(delta, 0)
            upper_extra += max(delta, 0)
        else:
            known += delta * call
        witness.append(
            {
                "gene_id": gene,
                "stable_weight": a,
                "comparator_weight": b,
                "coefficient": delta,
                "binary_call": call,
                "assay_state": "not_measured"
                if (call is None)
                else "eligible_measured",
            }
        )
    return {
        "observed_difference": known,
        "lower": known + lower_extra,
        "upper": known + upper_extra,
        "coefficients": witness,
    }

def summarize(rows):
    return {
        "full_count": len(rows),
        "available_count": sum(
            (row["observed_difference"] is not None for (row) in (rows))
        ),
        "observed_difference": mean([row["observed_difference"] for (row) in (rows)]),
        "lower": mean([row["lower"] for (row) in (rows)]),
        "upper": mean([row["upper"] for (row) in (rows)]),
        "stable": {
            key: mean([row["stable"][key] for (row) in (rows)]) for (key) in (METRICS)
        },
        "comparator": {
            key: mean([row["comparator"][key] for (row) in (rows)])
            for (key) in (METRICS)
        },
    }

def quantile(values, q):
    values = sorted(values)
    position = (len(values) - 1) * q
    lower, upper = (math.floor(position), math.ceil(position))
    return values[lower] + (values[upper] - values[lower]) * (position - lower)

def uncertainty(vectors, indices, production = True):
    from .matched_uncertainty import uncertainty as exact_uncertainty

    return exact_uncertainty(vectors, indices, production)

def analyze(
    groups,
    labels,
    ranks,
    eligibility,
    split,
    bootstrap,
    consensus,
    checkpoint,
    production = True,
):
    models = roster(split, production)
    seed_rows, model_rows = ([], [])
    for (model) in (models):
        checkpoint()
        sample, donor = (model["sample_ID"], model["study_donor_id"])
        eligible = eligibility["pools"]["all34"][sample]
        native, common = (
            set(eligible["native_eligible_literal_gene_labels"]),
            set(eligible["common_eligible_literal_gene_labels"]),
        )
        for (method) in (METHODS):
            collected = {name: [] for (name) in (COMPARATORS)}
            for (seed) in (seeds(method)):
                policy_weights = {
                    policy: None
                    if (groups[sample, method, seed, policy]["state"] == "unavailable")
                    else groups[sample, method, seed, policy]["weights"]
                    for (policy) in (POLICIES)
                }
                stable, latest, rc = (policy_weights[policy] for (policy) in (POLICIES))
                require(
                    all(
                        (
                            weights is None or set(weights) <= axis
                            for ((weights, axis)) in (
                                ((stable, native), (latest, native), (rc, common))
                            )
                        )
                    ),
                    "INHERITED_ELIGIBILITY_CHANGED",
                )
                target = None if (stable is None) else sum(stable.values(), Fraction(0))
                raw = [
                    ranks[method, "all34", "native_" + release, seed, sample]
                    for (release) in (("11_0", "11_5", "12_0"))
                ]
                latest_scores = {} if (raw[2][0] is None) else dict(zip(*raw[2]))
                consensus_result = consensus(raw, common)
                rc_scores = (
                    {}
                    if (consensus_result["genes"] is None)
                    else dict(
                        zip(consensus_result["genes"], consensus_result["scores"])
                    )
                )
                for (comparator) in (COMPARATORS):
                    q = None
                    if (stable is None):
                        right, status = (None, "STABLE_UNAVAILABLE")
                    elif (comparator == "THINNING_EXPECTATION"):
                        if (latest is None):
                            right, status = (None, "COMPARATOR_UNAVAILABLE")
                        else:
                            mass = sum(latest.values(), Fraction(0))
                            require(target <= mass, "STABLE_NOT_WITHIN_LATEST_MASS")
                            q = target / mass if (mass) else Fraction(0)
                            right, status = (
                                {
                                    gene: q * weight
                                    for ((gene, weight)) in (latest.items())
                                    if (q * weight)
                                },
                                "AVAILABLE_EXPECTATION",
                            )
                    else:
                        right, status = prefix(
                            latest if (comparator == "LATEST_PREFIX_R") else rc,
                            latest_scores
                            if (comparator == "LATEST_PREFIX_R")
                            else rc_scores,
                            target,
                        )
                    base = {
                        "model_id": sample,
                        "donor_id": donor,
                        "split": model["split"],
                        "method": method,
                        "seed": seed,
                        "comparison": comparator,
                        "status": status,
                        "target_mass": target,
                        "thinning_q": q,
                        "stable_weights": stable,
                        "comparator_weights": right,
                        "stable": metrics(stable, labels[sample]),
                        "comparator": metrics(right, labels[sample]),
                    }
                    seed_rows.append(base)
                    collected[comparator].append(base)
            for (comparator, rows) in (collected.items()):
                stable = average_weights([row["stable_weights"] for (row) in (rows)])
                right = average_weights([row["comparator_weights"] for (row) in (rows)])
                bound = sharp(stable, right, labels[sample])
                model_rows.append(
                    {
                        "model_id": sample,
                        "donor_id": donor,
                        "split": model["split"],
                        "method": method,
                        "comparison": comparator,
                        "seed_count": len(rows),
                        "seed_statuses": [row["status"] for (row) in (rows)],
                        "stable": metrics(stable, labels[sample]),
                        "comparator": metrics(right, labels[sample]),
                    }
                    | bound
                )
    donor_rows, summaries = ([], [])
    for (donor) in (split["ordered_donors"]):
        for (method) in (METHODS):
            for (comparator) in (COMPARATORS):
                rows = [
                    row
                    for (row) in (model_rows)
                    if (
                        (row["donor_id"], row["method"], row["comparison"])
                        == (donor, method, comparator)
                    )
                ]
                require(
                    rows and len({row["split"] for (row) in (rows)}) == 1,
                    "DONOR_GROUP_INVALID",
                )
                donor_rows.append(
                    {
                        "donor_id": donor,
                        "split": rows[0]["split"],
                        "method": method,
                        "comparison": comparator,
                        "model_ids": [row["model_id"] for (row) in (rows)],
                    }
                    | summarize(rows)
                )
    for (label) in (("held_out", "development", "all")):
        for (method) in (METHODS):
            for (comparator) in (COMPARATORS):
                rows = [
                    row
                    for (row) in (donor_rows)
                    if (
                        row["method"] == method
                        and row["comparison"] == comparator
                        and (label == "all" or row["split"] == label)
                    )
                ]
                summaries.append(
                    {
                        "stratum": label,
                        "method": method,
                        "comparison": comparator,
                        "weighting": "equal_models_within_donor_then_equal_donors",
                        "bounds_kind": "sharp_shared_model_gene_identification_not_confidence_interval",
                    }
                    | summarize(rows)
                )
    draw = bootstrap["splits"]["held_out"]
    held_donors = [
        donor
        for (donor) in (split["ordered_donors"])
        if (
            any(
                (
                    row["study_donor_id"] == donor and row["split"] == "held_out"
                    for (row) in (models)
                )
            )
        )
    ]
    require(draw["donor_order"] == held_donors, "PAIRED_DONOR_ORDER_CHANGED")
    donor_lookup = {
        (row["donor_id"], row["method"], row["comparison"]): row
        for (row) in (donor_rows)
    }
    vectors = {
        method + "::" + comparator: [
            donor_lookup[donor, method, comparator]["observed_difference"]
            for (donor) in (held_donors)
        ]
        for (method) in (METHODS)
        for (comparator) in (COMPARATORS[:2])
    }
    intervals = uncertainty(vectors, draw["indices"], production)
    if (production):
        require(
            (len(seed_rows), len(model_rows), len(donor_rows), len(summaries))
            == (1785, 765, 747, 27),
            "NEW_OUTPUT_SCOPE_CHANGED",
        )
    return {
        "status": "POSTHOC_NEW_EQUAL_RETURN_AND_SHARED_LABEL_DIAGNOSTICS_NOT_ACCEPTANCE",
        "normal_reference_variant": "all34",
        "support": "native_with_inherited_RC_common_axis",
        "offered_k": 10,
        "original_primary_unchanged": True,
        "independent_confirmation": False,
        "seed_comparisons": seed_rows,
        "model_comparisons": model_rows,
        "donor_comparisons": donor_rows,
        "summaries": summaries,
        "heldout_uncertainty": intervals,
    }
