import random
from fractions import Fraction

import numpy as np

from ..transfer.policies import (
    ENDPOINTS,
    PRIMARY,
    evaluate_means,
    mean,
    ordinary_parameter_diagnostic,
)

class MeanBank:
    def __init__(self, contexts, populations, biological_units):
        if ((set(contexts) != set(populations)) or (
            set(contexts) != set(biological_units)
        )):
            raise ValueError("Context inventory mismatch")
        self.contexts = contexts
        self.populations = populations
        self.biological_units = biological_units
        self.keys = None
        self.arrays = {}
        self.available = {}
        observed_units = set()
        for (context, groups) in (contexts.items()):
            population = populations[context]
            units = biological_units[context]
            if (
                (len(population) != len(set(population)))
                or (not population)
                or (len(units) != len(population))
                or (len(units) != len(set(units)))
            ):
                raise ValueError(
                    "Repeated or missing biological units require a grouped design, not independent rows"
                )
            if (observed_units & set(units)):
                raise ValueError(
                    "Overlapping biological units across contexts require root design reconciliation"
                )
            observed_units.update(units)
            keys = [
                (endpoint, candidate)
                for (endpoint) in (ENDPOINTS)
                for (candidate) in (sorted(groups[endpoint]))
            ]
            if (self.keys is None):
                self.keys = keys
            if (keys != self.keys):
                raise ValueError("Candidate or endpoint inventory differs")
            columns = [groups[endpoint][candidate] for ((endpoint, candidate)) in (keys)]
            if (any(len(values) != len(population) for (values) in (columns))):
                raise ValueError("Sample denominator differs")
            if (any(
                (value is not None) and (not isinstance(value, Fraction))
                for (values) in (columns)
                for (value) in (values)
            )):
                raise ValueError("Exact rational endpoint values are required")
            self.arrays[context] = np.asarray(
                [
                    [0.0 if (value is None) else float(value) for (value) in (values)]
                    for (values) in (columns)
                ],
                dtype = float,
            ).T
            self.available[context] = np.asarray(
                [
                    all(value is not None for (value) in (values))
                    for (values) in (columns)
                ]
            )

    def means(self, context, counts = None):
        n = len(self.populations[context])
        result = {endpoint: {} for (endpoint) in (ENDPOINTS)}
        if (counts is None):
            for (endpoint, candidate) in (self.keys):
                result[endpoint][candidate] = mean(
                    self.contexts[context][endpoint][candidate]
                )
            return result
        if (
            (len(counts) != n)
            or (sum(counts) != n)
            or any((not isinstance(count, int)) or (count < 0) for (count) in (counts))
        ):
            raise ValueError(
                "Resampling weights differ from the full biological population"
            )
        values = np.asarray(counts, dtype = float) @ self.arrays[context] / n
        for (index, (endpoint, candidate)) in (enumerate(self.keys)):
            if (not self.available[context][index]):
                value = None
            elif (endpoint == PRIMARY):
                source = self.contexts[context][endpoint][candidate]
                value = (
                    sum(
                        (
                            source[position] * count
                            for ((position, count)) in (enumerate(counts))
                            if (count)
                        ),
                        Fraction(0),
                    )
                    / n
                )
            else:
                value = float(values[index])
            result[endpoint][candidate] = value
        return result

def combine(context_means, contexts):
    if (not contexts):
        raise ValueError("No development contexts")
    first = context_means[contexts[0]]
    return {
        endpoint: {
            candidate: mean(
                [
                    context_means[context][endpoint][candidate]
                    for (context) in (contexts)
                ]
            )
            for (candidate) in (first[endpoint])
        }
        for (endpoint) in (ENDPOINTS)
    }

def percentile(values, probability):
    if (not values):
        return None
    position = (len(values) - 1) * probability
    lower = int(position)
    upper = min((lower + 1), (len(values) - 1))
    return values[lower] + ((position - lower) * (values[upper] - values[lower]))

def diagnostic_metrics(records):
    fields = (
        "development",
        "target",
        "development_gain",
        "target_gain",
        "gain_transfer_difference",
        "loss_same_three_trial_menu",
    )
    return {
        "::".join(map(str, (*row["endpoint"], row["method"], row["arm"], field))): row[
            field
        ]
        for (row) in (records)
        for (field) in (fields)
    }

def evaluate(bank, designs, draws = 2000, seed = 20260927, ordinary_bank = None):
    if (draws < 2):
        raise ValueError("At least two synthetic or 2000 production draws are required")
    observed = {context: bank.means(context) for (context) in (bank.contexts)}
    points = {}
    samples = {}
    choice_counts = {}
    ordinary_points, ordinary_samples = {}, {}
    ordinary_designs = [
        design
        for (design) in (designs)
        if (
            ordinary_bank
            and all(
                context in ordinary_bank.contexts
                for (context) in ((*design[1], design[2]))
            )
        )
    ]
    if (ordinary_bank):
        if (set(ordinary_bank.contexts) != {"COAD_CCLE", "LUAD_CCLE"}):
            raise ValueError(
                "Ordinary-parameter diagnostic must not include target contexts"
            )
        for (context) in (ordinary_bank.contexts):
            if ((ordinary_bank.populations[context] != bank.populations[context]) or (
                ordinary_bank.biological_units[context]
                != bank.biological_units[context]
            )):
                raise ValueError(
                    "Both diagnostic arms must use the identical paired development populations"
                )
        ordinary_observed = {
            context: ordinary_bank.means(context)
            for (context) in (ordinary_bank.contexts)
        }
        for (name, development, target) in (ordinary_designs):
            ordinary_points[name] = ordinary_parameter_diagnostic(
                combine(ordinary_observed, development), ordinary_observed[target]
            )
            ordinary_samples[name] = {
                key: [] for (key) in (diagnostic_metrics(ordinary_points[name]))
            }
        ordinary_points["combined_development_resubstitution_not_transfer"] = (
            ordinary_parameter_diagnostic(
                combine(ordinary_observed, tuple(sorted(ordinary_observed))),
                combine(ordinary_observed, tuple(sorted(ordinary_observed))),
            )
        )
    for (name, development, target) in (designs):
        if ((target in development) or any(
            context not in bank.contexts for (context) in ((*development, target))
        )):
            raise ValueError("Invalid development/target split")
        points[name] = evaluate_means(combine(observed, development), observed[target])
        samples[name] = {key: [] for (key) in (points[name]["metrics"])}
        choice_counts[name] = {}
    rng = random.Random(seed)
    for (draw) in (range(draws)):
        current, ordinary_current = {}, {}
        for (context) in (sorted(bank.contexts)):
            n = len(bank.populations[context])
            counts = [0] * n
            for (position) in (range(n)):
                counts[rng.randrange(n)] += 1
            current[context] = bank.means(context, counts)
            if (ordinary_bank and (context in ordinary_bank.contexts)):
                ordinary_current[context] = ordinary_bank.means(context, counts)
        for (name, development, target) in (designs):
            result = evaluate_means(combine(current, development), current[target])
            choice = str(
                (result["choices"]["JVS_choice"], result["choices"]["policies"]["MVS"])
            )
            choice_counts[name][choice] = choice_counts[name].get(choice, 0) + 1
            for (key, value) in (result["metrics"].items()):
                if (value is not None):
                    samples[name][key].append(float(value))
        for (name, development, target) in (ordinary_designs):
            values = ordinary_parameter_diagnostic(
                combine(ordinary_current, development), ordinary_current[target]
            )
            for (key, value) in (diagnostic_metrics(values).items()):
                if (value is not None):
                    ordinary_samples[name][key].append(float(value))
    intervals = {}
    for (name, records) in (samples.items()):
        intervals[name] = {}
        for (key, values) in (records.items()):
            values.sort()
            intervals[name][key] = {
                "draws": draws,
                "available_draws": len(values),
                "lower_95": percentile(values, 0.025),
                "upper_95": percentile(values, 0.975),
                "status": (
                    "CONDITIONAL_FIXED_CONTEXT_DEVELOPMENT_RESELECTION"
                    if (len(values) == draws)
                    else "UNAVAILABLE_DRAWS_RETAINED"
                ),
            }
    ordinary_intervals = {}
    for (name, records) in (ordinary_samples.items()):
        ordinary_intervals[name] = {}
        for (key, values) in (records.items()):
            values.sort()
            ordinary_intervals[name][key] = {
                "draws": draws,
                "available_draws": len(values),
                "lower_95": percentile(values, 0.025),
                "upper_95": percentile(values, 0.975),
            }
    return {
        "points": points,
        "intervals": intervals,
        "choice_frequencies": choice_counts,
        "ordinary_parameter_points": ordinary_points,
        "ordinary_parameter_intervals": ordinary_intervals,
        "draws": draws,
        "seed": seed,
        "development_reselection_in_each_draw": True,
        "target_used_for_choice": False,
        "paired_all_policies": True,
        "original_exposed_TCGA_not_confirmatory": True,
    }
