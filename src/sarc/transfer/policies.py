import itertools
import math
from fractions import Fraction

METHODS = ("DawnRank", "PRODIGY", "PersonaDrive")
VERSIONS = ("11_0", "11_5", "12_0")
PAIRS = (
    ("PersonaDrive", "DawnRank"),
    ("PersonaDrive", "PRODIGY"),
    ("DawnRank", "PRODIGY"),
)
REFERENCES = (
    "NCG6_primary_all",
    "NCG6_known_subgroup",
    "NCG6_candidate_subgroup",
    "CancerMine2019_secondary",
)
ENDPOINTS = tuple(
    (
        (support, reference, k)
        for (support) in (("common", "native"))
        for (reference) in (REFERENCES)
        for (k) in ((10, 1, 5, 20))
    )
)
PRIMARY = ("common", "NCG6_primary_all", 10)
SEEDS = (104729, 130363, 155921, 196613, 228017)
PARAMETERS = {"DawnRank": "mu_3", "PRODIGY": "alpha_0.05", "PersonaDrive": "original"}

def arithmetic(first, second, sign = -1):
    return None if (first is None or second is None) else first + sign * second

def mean(values):
    return (
        None
        if (not values or any((value is None for (value) in (values))))
        else sum(values) / len(values)
    )

def tie_order(candidate):
    (method, version) = candidate
    return (VERSIONS.index(version) if (version in VERSIONS) else len(VERSIONS), method)

def ordered(scores):
    available = [
        candidate for ((candidate, value)) in (scores.items()) if (value is not None)
    ]
    return sorted(
        available, key = lambda candidate: (-scores[candidate], tie_order(candidate))
    )

def consensus(rankings, eligible):
    if (len(rankings) != 3):
        raise ValueError("Three release lists are required")
    eligible = set(eligible)
    losses = []
    for (genes, scores) in (rankings):
        if (genes is None):
            if (scores is not None):
                raise ValueError("Unavailable list exposes scores")
            return {
                "status": "UNAVAILABLE",
                "genes": None,
                "scores": None,
                "rank_losses": None,
            }
        if (scores is None or len(genes) != len(scores) or len(set(genes)) != len(genes)):
            raise ValueError("Malformed ranked list")
        if (any((not math.isfinite(score) for (score) in (scores))) or any(
            (left < right for ((left, right)) in (zip(scores, scores[1:])))
        )):
            raise ValueError("Scores must be finite and decreasing")
        retained = [
            (gene, score)
            for ((gene, score)) in (zip(genes, scores))
            if (gene in eligible)
        ]
        n = len(retained)
        used = 0
        release = {}
        for (score, entries) in (itertools.groupby(retained, key = lambda pair: pair[1])):
            group = list(entries)
            loss = Fraction(2 * used + len(group), 2 * n)
            release.update(((gene, loss) for ((gene, value)) in (group)))
            used += len(group)
        losses.append(release)
    genes = set().union(*(set(values) for (values) in (losses)))
    combined = {
        gene: sum((values.get(gene, Fraction(1)) for (values) in (losses)), Fraction(0))
        / 3
        for (gene) in (genes)
    }
    result = sorted(genes, key = lambda gene: (combined[gene], gene))
    return {
        "status": "SUCCESS" if (result) else "SUCCESS_EMPTY",
        "genes": result,
        "scores": [-combined[gene] for (gene) in (result)],
        "rank_losses": combined,
    }

def select_policies(development_primary, methods = METHODS):
    methods = tuple(methods)
    if (methods not in (METHODS, ("DawnRank", "PersonaDrive"))):
        raise ValueError(
            "Only the full primary or explicitly approved binary transfer portfolio is admitted"
        )
    required = {
        (method, version) for (method) in (methods) for (version) in ((*VERSIONS, "RC"))
    }
    if (set(development_primary) != required):
        raise ValueError("Complete three-method policy menu is required")
    policies = {
        "V" + str(index): {method: (method, version) for (method) in (methods)}
        for ((index, version)) in (enumerate(VERSIONS))
    }
    chosen = {}
    for (method) in (methods):
        candidates = ordered(
            {
                (method, version): development_primary[method, version]
                for (version) in (VERSIONS)
            }
        )
        chosen[method] = candidates[0] if (candidates) else None
    policies["MVS"] = chosen
    policies["RC"] = {method: (method, "RC") for (method) in (methods)}
    leaders = {}
    for (policy, mapping) in (policies.items()):
        candidates = ordered(
            {
                candidate: development_primary[candidate]
                for (candidate) in (mapping.values())
                if (candidate is not None)
            }
        )
        leaders[policy] = candidates
    joint = ordered(
        {
            candidate: value
            for ((candidate, value)) in (development_primary.items())
            if (candidate[1] in VERSIONS)
        }
    )
    if (bool(joint) != bool(leaders["MVS"]) or (joint and joint[0] != leaders["MVS"][0])):
        raise ValueError("Joint and per-method search winners differ")
    return {
        "policies": policies,
        "development_order": leaders,
        "JVS_choice": joint[0] if (joint) else None,
        "aliases": {"latest": "V2", "JVS_winner": "MVS_winner"},
    }

def binary_transfer_export(development_primary):
    methods = ("DawnRank", "PersonaDrive")
    selected = {
        candidate: value
        for ((candidate, value)) in (development_primary.items())
        if (candidate[0] in methods)
    }
    return {
        "role": "SECONDARY_BINARY_APPLICABILITY_PORTFOLIO_FOR_XE",
        "methods": methods,
        "primary_three_method_result_unchanged": True,
        "not_substitution_for_PRODIGY": True,
        "choice_inputs": "STRING_DEVELOPMENT_PRIMARY_ONLY",
        "selection": select_policies(selected, methods),
    }

def portfolio_loss(values, chosen):
    if (
        chosen is None
        or chosen not in values
        or any((value is None for (value) in (values.values())))
    ):
        return None
    return max(values.values()) - values[chosen]

def evaluate_means(development, target, exact_development_primary = None):
    if (set(development) != set(ENDPOINTS) or set(target) != set(ENDPOINTS)):
        raise ValueError("All 32 endpoint projections are required")
    choices = select_policies(
        development[PRIMARY]
        if (exact_development_primary is None)
        else exact_development_primary
    )
    metrics = {}
    for (endpoint) in (ENDPOINTS):
        prefix = "::".join(map(str, endpoint))
        (source, destination) = (development[endpoint], target[endpoint])
        if (set(source) != set(development[PRIMARY]) or set(destination) != set(source)):
            raise ValueError("Projection-specific candidate menu changed")
        native = {
            candidate: value
            for ((candidate, value)) in (destination.items())
            if (candidate[1] in VERSIONS)
        }
        for (policy, mapping) in (choices["policies"].items()):
            local = {}
            for (method, candidate) in (mapping.items()):
                left = source.get(candidate)
                right = destination.get(candidate)
                local[method] = (left, right)
                metrics[prefix + "::" + policy + "::" + method + "::development"] = left
                metrics[prefix + "::" + policy + "::" + method + "::target"] = right
            candidates = choices["development_order"][policy]
            choice = candidates[0] if (candidates) else None
            menu = {
                candidate: destination[candidate]
                for (candidate) in (mapping.values())
                if (candidate is not None)
            }
            metrics[prefix + "::" + policy + "::loss_same_policy"] = (
                portfolio_loss(menu, choice) if (len(menu) == 3) else None
            )
            metrics[prefix + "::" + policy + "::loss_nine_native"] = (
                portfolio_loss(native, choice) if (policy != "RC") else None
            )
            if (policy == "RC"):
                metrics[prefix + "::RC::loss_enlarged_twelve"] = portfolio_loss(
                    destination, choice
                )
            for (first, second) in (PAIRS):
                old_gap = arithmetic(local[first][0], local[second][0])
                new_gap = arithmetic(local[first][1], local[second][1])
                pair = prefix + "::" + policy + "::" + first + "-" + second
                metrics[pair + "::development_gap"] = old_gap
                metrics[pair + "::target_gap"] = new_gap
                metrics[pair + "::gap_transfer_difference"] = arithmetic(
                    old_gap, new_gap
                )
                metrics[pair + "::nominal_reversal"] = (
                    None
                    if (old_gap is None or new_gap is None)
                    else int(old_gap * new_gap < 0)
                )
                for (delta) in ((Fraction(1, 2), Fraction(1), Fraction(2))):
                    value = (
                        None
                        if (old_gap is None or new_gap is None)
                        else int(
                            old_gap > delta
                            and new_gap < -delta
                            or (old_gap < -delta and new_gap > delta)
                        )
                    )
                    metrics[pair + "::meaningful_reversal_" + str(delta)] = value
        ranking = choices["development_order"]["MVS"]
        (winner, runner) = ranking[:2] if (len(ranking) >= 2) else (None, None)
        apparent = arithmetic(source.get(winner), source.get(runner))
        transferred = arithmetic(destination.get(winner), destination.get(runner))
        metrics[prefix + "::optimism::apparent_selected_advantage"] = apparent
        metrics[prefix + "::optimism::transferred_same_choices"] = transferred
        metrics[prefix + "::optimism::difference"] = arithmetic(apparent, transferred)
        for (version) in (VERSIONS):
            first = (winner[0], version) if (winner) else None
            second = (runner[0], version) if (runner) else None
            fixed_source = arithmetic(source.get(first), source.get(second))
            fixed_target = arithmetic(destination.get(first), destination.get(second))
            development_gain = arithmetic(apparent, fixed_source)
            target_gain = arithmetic(transferred, fixed_target)
            base = prefix + "::excess::" + version
            metrics[base + "::fixed_development_gap"] = fixed_source
            metrics[base + "::fixed_target_gap"] = fixed_target
            metrics[base + "::development_gain"] = development_gain
            metrics[base + "::target_gain"] = target_gain
            metrics[base + "::optimism"] = arithmetic(development_gain, target_gain)
    return {"choices": choices, "metrics": metrics}

def ordinary_parameter_diagnostic(training, target):
    if (set(training) != set(ENDPOINTS) or set(target) != set(ENDPOINTS)):
        raise ValueError("All endpoint projections must be retained")
    grids = {
        "DawnRank": ("mu_1", "mu_3", "mu_10"),
        "PRODIGY": ("alpha_0.01", "alpha_0.05", "alpha_0.1"),
    }
    result = []
    for (method, grid) in (grids.items()):
        version_menu = [
            (method, version, PARAMETERS[method]) for (version) in (VERSIONS)
        ]
        parameter_menu = [(method, "12_0", parameter) for (parameter) in (grid)]
        for (arm, menu) in ((("version", version_menu), ("parameter", parameter_menu))):
            if (any(
                (
                    candidate not in training[endpoint]
                    or candidate not in target[endpoint]
                    for (endpoint) in (ENDPOINTS)
                    for (candidate) in (menu)
                )
            )):
                raise ValueError("Missing planned ordinary-budget candidate")
            available = [
                candidate
                for (candidate) in (menu)
                if (training[PRIMARY][candidate] is not None)
            ]

            def order(candidate):
                if (arm == "version"):
                    return (VERSIONS.index(candidate[1]), "")
                return (int(candidate[2] != PARAMETERS[method]), candidate[2])

            selected = (
                min(
                    available,
                    key = lambda candidate: (
                        -training[PRIMARY][candidate],
                        order(candidate),
                    ),
                )
                if (available)
                else None
            )
            for (endpoint) in (ENDPOINTS):
                source = training[endpoint].get(selected)
                destination = target[endpoint].get(selected)
                baseline = (method, "12_0", PARAMETERS[method])
                source_gain = arithmetic(source, training[endpoint].get(baseline))
                target_gain = arithmetic(destination, target[endpoint].get(baseline))
                loss = portfolio_loss(
                    {candidate: target[endpoint][candidate] for (candidate) in (menu)},
                    selected,
                )
                result.append(
                    {
                        "method": method,
                        "arm": arm,
                        "endpoint": endpoint,
                        "choice": selected,
                        "opportunities": 3,
                        "master_seed": 104729 if (method == "PRODIGY") else None,
                        "development": source,
                        "target": destination,
                        "development_gain": source_gain,
                        "target_gain": target_gain,
                        "gain_transfer_difference": arithmetic(
                            source_gain, target_gain
                        ),
                        "loss_same_three_trial_menu": loss,
                        "equal_cpu_cost_claim": False,
                        "TCGA_parameter_transfer": False,
                    }
                )
    return result
