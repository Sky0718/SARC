import math
from fractions import Fraction

import numpy as np

from ..transfer import scoring as core
from . import review_rules as rules
from .constants import METHODS

class PairedBootstrap:
    def __init__(self, population):
        self.population = tuple(population)
        if (not self.population or len(set(self.population)) != len(self.population)):
            raise ValueError("A unique full ordered population is required")
        self.indices = core.bootstrap_indices(len(population), 2000, 20260926)
        self.counts = np.asarray(
            [np.bincount(draw, minlength = len(population)) for (draw) in (self.indices)],
            dtype = np.int64,
        )
        if (not np.all(self.counts.sum(axis = 1) == len(population))):
            raise ValueError("Paired bootstrap indices do not close")
        self.float_counts = self.counts.astype(np.float64)

    def sum_vectors(self, vectors):
        matrix = np.asarray(vectors, dtype = np.float64)
        if (matrix.ndim == 1):
            matrix = matrix[:, None]
        if (matrix.shape[0] != len(self.population) or not np.isfinite(matrix).all()):
            raise ValueError("Bootstrap vector differs from the complete sample axis")
        return self.float_counts @ matrix

    def integer_vectors(self, vectors):
        matrix = np.asarray(vectors, dtype = np.int64)
        if (matrix.ndim == 1):
            matrix = matrix[:, None]
        if (
            matrix.shape[0] != len(self.population)
            or np.any(matrix < 0)
            or np.any(matrix > 1)
        ):
            raise ValueError("Exact event counting requires full binary vectors")
        result = self.float_counts @ matrix.astype(np.float64)
        if (
            np.any(result != np.floor(result))
            or np.any(result < 0)
            or np.any(result > len(self.population))
        ):
            raise ValueError("Event sums must remain exact bounded integers")
        return result.astype(np.int64)

    def interval(self, values, status = "CONDITIONAL_FIXED_DECISION_PERCENTILE"):
        values = np.asarray(values, dtype = np.float64)
        if (values.shape != (2000,)):
            raise ValueError("Exactly 2,000 frozen bootstrap draws are required")
        available = values[np.isfinite(values)]
        complete = len(available) == 2000
        ordered = sorted(available.tolist()) if (complete) else []
        return {
            "available_draws": len(available),
            "draws": 2000,
            "lower_95": core.percentile(ordered, 0.025) if (complete) else None,
            "upper_95": core.percentile(ordered, 0.975) if (complete) else None,
            "status": status
            if (complete)
            else "UNAVAILABLE_ZERO_DENOMINATOR_OR_MISSING_DRAW",
        }

    def contrast_interval(self, result):
        rows = result["sample_contrasts"]
        if (tuple(row["sample_id"] for (row) in (rows)) != self.population):
            raise ValueError("Contrast ordering differs from frozen paired positions")
        values = [row["hits"] for (row) in (rows)]
        if (any(value is None for (value) in (values))):
            return {
                "lower_95": None,
                "upper_95": None,
                "status": "UNAVAILABLE_FULL_POPULATION_CONTRAST",
            }
        means = self.sum_vectors([float(value) for (value) in (values)])[:, 0] / len(
            values
        )
        result = self.interval(means, "CONDITIONAL_DESCRIPTIVE_PERCENTILE")
        return {
            name: result[name]
            for (name) in (("lower_95", "upper_95", "status", "draws"))
        }

    def selection_intervals(self, rows, choices):
        if (set(rows) != set(METHODS) or any(
            tuple(rows[method]) != self.population for (method) in (METHODS)
        )):
            raise ValueError(
                "Fixed complete portfolio and paired sample order required"
            )
        means = {}
        lower = {}
        upper = {}
        complete = np.ones(2000, dtype = bool)
        for (method) in (METHODS):
            records = [rows[method][sample] for (sample) in (self.population)]
            failures = self.integer_vectors(
                [int(row["hits"] is None) for (row) in (records)]
            )[:, 0]
            value = self.sum_vectors(
                [
                    float(row["hits"]) if (row["hits"] is not None) else 0.0
                    for (row) in (records)
                ]
            )[:, 0] / len(records)
            value[failures != 0] = np.nan
            means[method] = value
            complete &= failures == 0
            lower[method] = self.sum_vectors(
                [float(row["lower"]) for (row) in (records)]
            )[:, 0] / len(records)
            upper[method] = self.sum_vectors(
                [float(row["upper"]) for (row) in (records)]
            )[:, 0] / len(records)
        best = np.max(
            np.column_stack([means[method] for (method) in (METHODS)]), axis = 1
        )
        losses = {}
        for (method) in (METHODS):
            alternatives = [name for (name) in (METHODS) if (name != method)]
            losses[method] = {
                "value": best - means[method],
                "lower": np.maximum(
                    0,
                    np.max(
                        np.column_stack([lower[name] for (name) in (alternatives)]),
                        axis = 1,
                    )
                    - upper[method],
                ),
                "upper": np.maximum(
                    0,
                    np.max(
                        np.column_stack([upper[name] for (name) in (alternatives)]),
                        axis = 1,
                    )
                    - lower[method],
                ),
            }
        result = {}
        for (strategy, choice) in (choices.items()):
            if (strategy == "uniform_random"):
                vectors = {
                    field: sum(losses[method][field] for (method) in (METHODS)) / 3
                    for (field) in (("value", "lower", "upper"))
                }
            elif (strategy == "posthoc_oracle"):
                vectors = {
                    "value": np.where(complete, 0.0, np.nan),
                    "lower": np.zeros(2000),
                    "upper": np.where(
                        complete,
                        0.0,
                        np.max(np.column_stack(list(upper.values())), axis = 1),
                    ),
                }
            elif (choice is None):
                vectors = {
                    "value": np.full(2000, np.nan),
                    "lower": np.zeros(2000),
                    "upper": np.max(np.column_stack(list(upper.values())), axis = 1),
                }
            elif (choice in METHODS):
                vectors = losses[choice]
            else:
                raise ValueError("Unknown frozen chosen method")
            result[strategy] = {
                "intervals": {
                    field: self.interval(vector)
                    for ((field, vector)) in (vectors.items())
                },
                "draw_count": 2000,
                "seed": 20260926,
                "parameters_refitted": False,
                "chosen_method_recomputed": False,
                "posthoc_oracle_recomputed_for_loss_only": True,
                "status": "CONDITIONAL_SAMPLE_COMPOSITION_NOT_TRAINED_RULE_GENERALISATION",
            }
        return result

    def review_intervals(self, rows, allocations, delta):
        if (tuple(row["sample_id"] for (row) in (rows)) != self.population):
            raise ValueError(
                "Review interval cannot reallocate duplicate sample identities"
            )
        events = [rules.event_state(row, delta) for (row) in (rows)]
        known = np.asarray(
            [event["event"] is not None for (event) in (events)], dtype = np.int64
        )
        event = np.asarray(
            [event["event"] is True for (event) in (events)], dtype = np.int64
        )
        lower = np.asarray(
            [value["event_lower"] for (value) in (events)], dtype = np.int64
        )
        upper = np.asarray(
            [value["event_upper"] for (value) in (events)], dtype = np.int64
        )
        actionable = np.asarray([row["actionable"] for (row) in (rows)], dtype = np.int64)
        vectors = [known, event, lower, upper, actionable]
        positions = {}
        for (name, allocation) in (allocations.items()):
            selected = frozenset(allocation["allocated"])
            if (not selected.issubset(self.population) or any(
                not row["actionable"]
                for (row) in (rows)
                if (row["sample_id"] in selected)
            )):
                raise ValueError(
                    "Frozen allocation includes unavailable or foreign sample"
                )
            inspect = np.asarray(
                [row["sample_id"] in selected for (row) in (rows)], dtype = np.int64
            )
            positions[name] = len(vectors)
            vectors.extend(
                (
                    inspect,
                    inspect * event,
                    inspect * lower,
                    inspect * upper,
                    inspect * (1 - known),
                    inspect * known * (1 - event),
                )
            )
        sums = self.integer_vectors(np.column_stack(vectors))
        known_draw, events_draw, lower_draw, upper_draw, actionable_draw = (
            sums[:, index] for (index) in (range(5))
        )
        n = len(self.population)

        def ratio(numerator, denominator):
            with np.errstate(divide = "ignore", invalid = "ignore"):
                return np.where(denominator != 0, numerator / denominator, np.nan)

        result = {}
        for (name, offset) in (positions.items()):
            (
                inspected,
                captured,
                captured_lower,
                captured_upper,
                unavailable_inspections,
                non_events,
            ) = (sums[:, offset + index] for (index) in (range(6)))
            counts = {
                "event_known": known_draw,
                "events_known": events_draw,
                "non_events_known": known_draw - events_draw,
                "events_unavailable": n - known_draw,
                "captured_known": captured,
                "missed_known": events_draw - captured,
                "non_event_inspections": non_events,
                "unavailable_inspections": unavailable_inspections,
                "events_lower": lower_draw,
                "events_upper": upper_draw,
                "captured_lower": captured_lower,
                "captured_upper": captured_upper,
                "missed_lower": lower_draw - captured_lower,
                "missed_upper": upper_draw - captured_upper,
            }
            vectors = {
                "captured_yield": np.where(
                    unavailable_inspections == 0, captured / n, np.nan
                ),
                "captured_yield_lower": captured_lower / n,
                "captured_yield_upper": captured_upper / n,
                "precision": np.where(
                    unavailable_inspections == 0, ratio(captured, inspected), np.nan
                ),
                "known_event_capture_rate": ratio(captured, events_draw),
                "full_event_capture_rate": np.where(
                    known_draw == n, ratio(captured, events_draw), np.nan
                ),
                "full_event_capture_lower": ratio(
                    captured_lower, captured_lower + counts["missed_upper"]
                ),
                "full_event_capture_upper": ratio(
                    captured_upper, captured_upper + counts["missed_lower"]
                ),
                "full_event_rate": np.where(known_draw == n, events_draw / n, np.nan),
            }
            result[name] = {
                "status": "CONDITIONAL_SAMPLE_COMPOSITION_NOT_TRAINED_RULE_GENERALISATION",
                "draw_count": 2000,
                "seed": 20260926,
                "allocation_recomputed": False,
                "parameters_refitted": False,
                "allocated_sample_ids": list(allocations[name]["allocated"]),
                "intervals": {
                    field: self.interval(
                        value, "CONDITIONAL_FIXED_ALLOCATION_PERCENTILE"
                    )
                    for ((field, value)) in ((vectors | counts).items())
                },
            }
        return result
