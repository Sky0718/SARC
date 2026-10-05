from fractions import Fraction

from .features import context_weights

METHODS = ("DawnRank", "PRODIGY", "PersonaDrive")
CONTEXTS = ("COAD_CCLE", "LUAD_CCLE")
TRANSITIONS = ("11_0_to_11_5", "11_5_to_12_0")
NATIVE = ("native_11_0", "native_11_5", "native_12_0")
RHOS = (Fraction(0), Fraction(1, 2), Fraction(1))
FAMILIES = ("strong", "candidate")
PRIMARY = ("common", "NCG6_primary_all", 10)

def mean(values):
    values = tuple(values)
    return (
        None
        if ((not values) or any(value is None for (value) in (values)))
        else sum(values, Fraction(0)) / len(values)
    )

def exact(value):
    if ((not isinstance(value, (int, Fraction))) or isinstance(value, bool)):
        raise ValueError("Exact rational performance is required")
    return Fraction(value)

def choose_scores(scores):
    if (set(scores) != set(METHODS)):
        raise ValueError("The fixed complete three-method portfolio is required")
    if (any(value is None for (value) in (scores.values()))):
        return {
            "choice": None,
            "status": "UNAVAILABLE_TRAINING_PORTFOLIO",
            "scores": scores,
        }
    for (value) in (scores.values()):
        exact(value)
    best = max(scores.values())
    return {
        "choice": next(method for (method) in (METHODS) if (scores[method] == best)),
        "status": "AVAILABLE",
        "scores": scores,
    }

def choose_rule(training_estimates, training_features, target_features, family, rho):
    if (
        (family not in FAMILIES)
        or (rho not in RHOS)
        or (set(training_estimates) != set(training_features))
        or (not training_estimates)
    ):
        raise ValueError("Invalid frozen rule family, parameter or training contexts")
    if (any(context not in CONTEXTS for (context) in (training_estimates))):
        raise ValueError("Training is restricted to the two development contexts")
    weights = context_weights(target_features, training_features)
    missing = False
    for (estimates) in (training_estimates.values()):
        if (set(estimates) != set(METHODS)):
            raise ValueError("Training cannot remove a method")
        for (states) in (estimates.values()):
            if (set(states) != {"a", "b", "c", "d"}):
                raise ValueError("Training needs all four actual input states")
            for (value) in (states.values()):
                if (value is None):
                    missing = True
                else:
                    exact(value)
    if (missing):
        return choose_scores(dict.fromkeys(METHODS)) | {
            "weights": weights,
            "family": family,
            "rho": rho,
            "estimates": None,
        }
    weighted = {
        method: {
            state: sum(
                (
                    weights[context] * training_estimates[context][method][state]
                    for (context) in (weights)
                ),
                Fraction(0),
            )
            for (state) in (("a", "b", "c", "d"))
        }
        for (method) in (METHODS)
    }
    scores = {}
    for (method, states) in (weighted.items()):
        a, b, c, d = (states[state] for (state) in (("a", "b", "c", "d")))
        risk = (min(a, b) + min(c, d)) / 2 if (family == "strong") else min(a, b, c, d)
        scores[method] = ((1 - rho) * b) + (rho * risk)
    return choose_scores(scores) | {
        "weights": weights,
        "family": family,
        "rho": rho,
        "estimates": weighted,
    }

def portfolio_loss(population, rows, chosen):
    population = tuple(population)
    if (
        (not population)
        or (len(population) != len(set(population)))
        or (set(rows) != set(METHODS))
        or (chosen not in METHODS and chosen is not None)
    ):
        raise ValueError("Invalid fixed portfolio loss domain")
    for (method) in (METHODS):
        if (tuple(rows[method]) != population):
            raise ValueError("Complete frozen ordered population is required")
    for (sample) in (population):
        signatures = {
            (
                frozenset(rows[method][sample]["eligible_ids"]),
                frozenset(rows[method][sample]["positive_ids"]),
                rows[method][sample]["k"],
            )
            for (method) in (METHODS)
        }
        if (len(signatures) != 1):
            raise ValueError("Method supports or task references differ")
    estimates = {}
    for (method) in (METHODS):
        current = [rows[method][sample] for (sample) in (population)]
        for (row) in (current):
            lower, upper = exact(row["lower"]), exact(row["upper"])
            if (
                (lower < 0)
                or (lower > upper)
                or (upper > min(row["k"], len(row["positive_ids"])))
            ):
                raise ValueError("Invalid attainable score bounds")
            if ((row["hits"] is not None) and (
                (exact(row["hits"]) != lower) or (lower != upper)
            )):
                raise ValueError("Successful performance and bounds disagree")
        estimates[method] = {
            "value": mean(row["hits"] for (row) in (current)),
            "lower": mean(row["lower"] for (row) in (current)),
            "upper": mean(row["upper"] for (row) in (current)),
            "failed_count": sum(row["hits"] is None for (row) in (current)),
        }
    complete = all(item["value"] is not None for (item) in (estimates.values()))
    if (chosen is None):
        return {
            "value": None,
            "lower": Fraction(0),
            "upper": max(item["upper"] for (item) in (estimates.values())),
            "selected_failed_count": None,
            "population_count": len(population),
            "portfolio_size": 3,
            "complete": complete,
            "estimates": estimates,
            "status": "UNAVAILABLE_TRAINING_PORTFOLIO",
        }
    alternatives = [method for (method) in (METHODS) if (method != chosen)]
    lower = max(
        Fraction(0),
        max(estimates[method]["lower"] for (method) in (alternatives))
        - estimates[chosen]["upper"],
    )
    upper = max(
        Fraction(0),
        max(estimates[method]["upper"] for (method) in (alternatives))
        - estimates[chosen]["lower"],
    )
    value = (
        max(item["value"] for (item) in (estimates.values()))
        - estimates[chosen]["value"]
        if (complete)
        else None
    )
    return {
        "value": value,
        "lower": lower,
        "upper": upper,
        "selected_failed_count": estimates[chosen]["failed_count"],
        "population_count": len(population),
        "portfolio_size": 3,
        "complete": complete,
        "estimates": estimates,
        "status": ("AVAILABLE" if (complete) else "UNAVAILABLE_REQUIRED_TEST_SCORE"),
    }

def expected_random_loss(population, rows):
    losses = [portfolio_loss(population, rows, method) for (method) in (METHODS)]
    return {
        "value": mean(item["value"] for (item) in (losses)),
        "lower": mean(item["lower"] for (item) in (losses)),
        "upper": mean(item["upper"] for (item) in (losses)),
        "expected_selected_failed_count": mean(
            item["selected_failed_count"] for (item) in (losses)
        ),
        "population_count": len(population),
        "portfolio_size": 3,
        "method_losses": dict(zip(METHODS, losses)),
        "setting": "EXACT_UNIFORM_EXPECTATION",
    }

def oracle_loss(population, rows):
    validation = portfolio_loss(population, rows, METHODS[0])
    choice = choose_scores(
        {method: item["value"] for ((method, item)) in (validation["estimates"].items())}
    )
    result = portfolio_loss(population, rows, choice["choice"])
    return {
        "choice": choice["choice"],
        "loss": result,
        "setting": "POSTHOC_FULL_PORTFOLIO_COST_ORACLE",
        "information": "All evaluation method outcomes; never used to select a rule",
    }

def state_networks(transition):
    if (transition not in TRANSITIONS):
        raise ValueError("Unregistered transition")
    index = TRANSITIONS.index(transition)
    return dict(
        zip(
            ("a", "b", "c", "d"),
            (
                NATIVE[index],
                NATIVE[index + 1],
                transition + "__persistent_old",
                transition + "__persistent_new",
            ),
        )
    )

def training_estimates(dataset, contexts, transition):
    if ((not contexts) or any(context not in CONTEXTS for (context) in (contexts))):
        raise ValueError("Invalid development training contexts")
    return {
        context: {
            method: {
                state: dataset["cohort_means"][(context, method, network, PRIMARY)]
                for ((state, network)) in (state_networks(transition).items())
            }
            for (method) in (METHODS)
        }
        for (context) in (contexts)
    }

def simple_rule(dataset, contexts, transition, baseline):
    if (baseline not in ("development_single_best", "previous_release_best")):
        raise ValueError("Unregistered simple baseline")
    networks = (
        NATIVE
        if (baseline == "development_single_best")
        else (state_networks(transition)["a"],)
    )
    return choose_scores(
        {
            method: mean(
                dataset["cohort_means"][(context, method, network, PRIMARY)]
                for (context) in (contexts)
                for (network) in (networks)
            )
            for (method) in (METHODS)
        }
    ) | {"baseline": baseline}

def endpoint_rows(dataset, context, network, endpoint):
    return {
        method: dataset["rows"][(context, method, network, endpoint)]
        for (method) in (METHODS)
    }

def portfolio_summary(population, rows):
    summaries = {}
    for (method) in (METHODS):
        current = [rows[method][sample] for (sample) in (population)]
        with_positives = [row for (row) in (current) if (row["positive_ids"])]
        summaries[method] = {
            "population_count": len(population),
            "successful_count": sum(row["hits"] is not None for (row) in (current)),
            "no_reference_positive_count": len(current) - len(with_positives),
            "eligible_positive_recall": mean(
                row["recall"] for (row) in (with_positives)
            ),
            "recall_denominator": len(with_positives),
            "native_output_count_mean": mean(
                row["native_output_count"] for (row) in (current)
            ),
            "eligible_output_count_mean": mean(
                row["common_output_count"] for (row) in (current)
            ),
            "reference_positive_output_count_mean": mean(
                row["positive_output_count"] for (row) in (current)
            ),
            "eligible_count_mean": mean(
                len(row["eligible_ids"]) for (row) in (current)
            ),
            "positive_count_mean": mean(
                len(row["positive_ids"]) for (row) in (current)
            ),
        }
    return summaries

def fitted_parameter(trials):
    if (tuple(trials) != RHOS):
        raise ValueError("All three frozen opportunities must be retained in order")
    means = {rho: mean(trials[rho]) for (rho) in (RHOS)}
    available = {rho: value for ((rho, value)) in (means.items()) if (value is not None)}
    if (not available):
        return {
            "rho": Fraction(1, 2),
            "status": "DEFAULT_NOT_TUNED",
            "losses": [{"rho": rho, "value": means[rho]} for (rho) in (RHOS)],
        }
    best = min(available.values())
    preferred = (Fraction(1, 2), Fraction(0), Fraction(1))
    return {
        "rho": next(rho for (rho) in (preferred) if (available.get(rho) == best)),
        "status": "DEVELOPMENT_SELECTED",
        "losses": [{"rho": rho, "value": means[rho]} for (rho) in (RHOS)],
    }

def evaluate_development(dataset, feature_records):
    if ((tuple(dataset["contexts"]) != CONTEXTS) or (
        set(feature_records) != set(CONTEXTS)
    )):
        raise ValueError("Exactly two development contexts are required")
    for (context) in (CONTEXTS):
        if (set(feature_records[context]) != set(TRANSITIONS)):
            raise ValueError("Input features must cover both transitions")
    trials = {family: {rho: [] for (rho) in (RHOS)} for (family) in (FAMILIES)}
    predictions = []
    baselines = []
    for (held_context) in (CONTEXTS):
        train = tuple(context for (context) in (CONTEXTS) if (context != held_context))
        for (transition) in (TRANSITIONS):
            new_network = state_networks(transition)["b"]
            population = dataset["populations"][held_context]
            primary_rows = endpoint_rows(dataset, held_context, new_network, PRIMARY)
            estimates = training_estimates(dataset, train, transition)
            features = {
                context: feature_records[context][transition] for (context) in (train)
            }
            for (family) in (FAMILIES):
                for (rho) in (RHOS):
                    prediction = choose_rule(
                        estimates,
                        features,
                        feature_records[held_context][transition],
                        family,
                        rho,
                    )
                    loss = portfolio_loss(
                        population, primary_rows, prediction["choice"]
                    )
                    trials[family][rho].append(loss["value"])
                    predictions.append(
                        {
                            "held_context": held_context,
                            "training_contexts": list(train),
                            "transition": transition,
                            "family": family,
                            "rho": rho,
                            "prediction": prediction,
                            "primary_loss": loss,
                        }
                    )
            for (baseline) in (("development_single_best", "previous_release_best")):
                prediction = simple_rule(dataset, train, transition, baseline)
                baselines.append(
                    {
                        "held_context": held_context,
                        "transition": transition,
                        "baseline": baseline,
                        "prediction": prediction,
                    }
                )
            baselines.extend(
                [
                    {
                        "held_context": held_context,
                        "transition": transition,
                        "baseline": "uniform_random",
                    },
                    {
                        "held_context": held_context,
                        "transition": transition,
                        "baseline": "posthoc_oracle",
                    },
                ]
            )
    fitted = {family: fitted_parameter(trials[family]) for (family) in (FAMILIES)}
    endpoints = []
    for (entry) in (predictions + baselines):
        if (("family" in entry) and (entry["rho"] != fitted[entry["family"]]["rho"])):
            continue
        context, transition = entry["held_context"], entry["transition"]
        population = dataset["populations"][context]
        for (endpoint) in (dataset["endpoints"]):
            rows = endpoint_rows(
                dataset, context, state_networks(transition)["b"], endpoint
            )
            strategy = entry.get("family", entry.get("baseline"))
            loss = (
                expected_random_loss(population, rows)
                if (strategy == "uniform_random")
                else oracle_loss(population, rows)
                if (strategy == "posthoc_oracle")
                else portfolio_loss(population, rows, entry["prediction"]["choice"])
            )
            endpoints.append(
                {
                    "held_context": context,
                    "transition": transition,
                    "strategy": strategy,
                    "endpoint": list(endpoint),
                    "rho": entry.get("rho"),
                    "choice": entry.get("prediction", {}).get("choice"),
                    "loss": loss,
                    "outcome_coverage": portfolio_summary(population, rows),
                }
            )
    final = {
        "training_contexts": list(CONTEXTS),
        "families": fitted,
        "features": feature_records,
        "transitions": {
            transition: {
                "estimates": training_estimates(dataset, CONTEXTS, transition),
                "simple_rules": {
                    baseline: simple_rule(dataset, CONTEXTS, transition, baseline)
                    for (baseline) in (
                        ("development_single_best", "previous_release_best")
                    )
                },
            }
            for (transition) in (TRANSITIONS)
        },
        "setting": "First analysis of an unseen complete cohort; no target outcomes or positive labels",
    }
    return {
        "schema": "s05_selection_development_v1",
        "predictions": predictions,
        "baselines": baselines,
        "fitted": fitted,
        "endpoints": endpoints,
        "final_model": final,
        "independent_generalisation_claim": False,
        "decision_count": 4,
    }

def predict_frozen(model, target_features, transition, family):
    current = model["transitions"][transition]
    features = {
        context: model["features"][context][transition]
        for (context) in (model["training_contexts"])
    }
    return choose_rule(
        current["estimates"],
        features,
        target_features,
        family,
        model["families"][family]["rho"],
    )

def conditional_selection_draw(population, rows, chosen, indices):
    population = tuple(population)
    indices = tuple(indices)
    if ((len(indices) != len(population)) or any(
        (type(index) is not int) or (not (0 <= index < len(population)))
        for (index) in (indices)
    )):
        raise ValueError("A complete paired sample-index draw is required")
    portfolio_loss(population, rows, chosen)
    repeated = tuple(str(index) for (index) in (range(len(indices))))
    drawn = {
        method: {
            label: rows[method][population[position]]
            for ((label, position)) in (zip(repeated, indices))
        }
        for (method) in (METHODS)
    }
    return portfolio_loss(repeated, drawn, chosen)
