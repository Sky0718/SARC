from pathlib import Path

import numpy as np
from scipy.sparse import csc_matrix

from .io import binary_matrix, read_axes, read_csv, read_json, require, resolve, write_json

def positive_minimum(values):
    selected = values[values > 0]
    return float(np.min(selected)) if (len(selected)) else 1.0

def gamma(count):
    value = np.asarray(count) * (2.0**-53)
    require(
        np.isfinite(value).all() and (value >= 0).all() and (value < 0.5).all(),
        "Roundoff count outside certificate domain",
    )
    return value / (1 - value)

def bound(matrix, damping, q, values):
    require(matrix.shape == (len(values), len(values)), "Operator dimensions differ")
    require(len(damping) == len(values) and len(q) == len(values), "Vector dimensions differ")
    require(
        all(np.isfinite(array).all() for (array) in ((matrix.data, damping, q, values))),
        "Nonfinite numerical state",
    )
    require(
        (matrix.data >= 0).all() and (q >= 0).all() and (values >= 0).all(), "Negative state"
    )
    require(
        (damping >= 0).all() and (damping < 1).all() and q.sum() > 0,
        "Invalid damping or restart mass",
    )
    row_counts = np.diff(matrix.indptr)
    col_counts = np.diff(matrix.tocsc().indptr)
    signal = damping * matrix.dot(values)
    restart = (1 - damping) * q
    residual = signal + restart - values
    allowance = gamma(row_counts + 3) / (1 - gamma(row_counts + 3)) * np.abs(signal)
    allowance += gamma(2) / (1 - gamma(2)) * np.abs(restart)
    allowance += gamma(3) / (1 - gamma(3)) * (np.abs(signal) + np.abs(restart) + np.abs(values))
    roundoff = np.sum(allowance) / (1 - gamma(len(values)))
    norms = matrix.multiply(damping[:, None]).tocsc().sum(axis = 0)
    upper = float(np.max(norms)) / (1 - gamma(int(max(col_counts)) + 3))
    minimum = min(
        positive_minimum(matrix.data) * positive_minimum(values) * positive_minimum(damping),
        positive_minimum(1 - damping) * positive_minimum(q),
    )
    require(
        np.isfinite(residual).all() and upper < 1 and minimum > np.finfo(float).tiny,
        "Arithmetic certificate unavailable",
    )
    error = (np.sum(np.abs(residual)) / (1 - gamma(len(values))) + roundoff) / (1 - upper)
    require(np.isfinite(error), "Nonfinite certificate")
    return float(error), float(upper)

def independent_iterations(matrix, damping, q, epsilon, limit):
    current = q.copy()
    restart = (1 - damping[:, None]) * q
    active = np.ones(q.shape[1], dtype = bool)
    iterations = np.zeros(q.shape[1], dtype = int)
    for (iteration) in (range(1, limit + 1)):
        columns = np.flatnonzero(active)
        if (not len(columns)):
            break
        updated = damping[:, None] * matrix.dot(current[:, columns]) + restart[:, columns]
        movement = np.sqrt(np.sum((updated - current[:, columns]) ** 2, axis = 0))
        require(
            np.isfinite(updated).all() and np.isfinite(movement).all(),
            "Nonfinite independent recurrence",
        )
        current[:, columns] = updated
        iterations[columns] = iteration
        stopped = columns[movement < epsilon]
        active[stopped] = False
    return current, iterations, ~active

def verify_directory(directory):
    directory = Path(directory)
    manifest = read_json(directory / "scores.json")
    genes, samples = read_axes(directory)
    n, m = len(genes), len(samples)
    query = binary_matrix(directory / "mutation.bin", n, m)
    require(np.isin(query, [0, 1]).all(), "Invalid query")
    rows = read_csv(directory / "cutoff_certificates.csv")
    diagnostics = {(row["state"], row["sample_id"], int(row["budget"])): row for (row) in (rows)}
    require(len(diagnostics) == len(rows), "Duplicate certificate")
    budgets = sorted({int(row["budget"]) for (row) in (rows)})
    expected = {
        (state, sample, budget)
        for (state) in (manifest["states"])
        for (sample) in (samples)
        for (budget) in (budgets)
    }
    require(set(diagnostics) == expected, "Incomplete certificate population")
    checks = []
    for (state, specification) in (manifest["states"].items()):
        graph = "G" + str(specification["a"])
        matrix = csc_matrix(
            (
                np.fromfile(directory / (graph + "_P_values.bin"), dtype = "<f8"),
                np.fromfile(directory / (graph + "_P_rows.bin"), dtype = "<i4"),
                np.fromfile(directory / (graph + "_P_pointers.bin"), dtype = "<i4"),
            ),
            shape = (n, n),
        ).tocsr()
        damping = np.fromfile(
            directory / ("G" + str(specification["b"]) + "_damping.bin"), dtype = "<f8"
        )
        q = binary_matrix(directory / ("H" + str(specification["h"]) + "_q.bin"), n, m)
        actual = binary_matrix(directory / (state + "_scores.bin"), n, m)
        current, iterations, stopped = independent_iterations(
            matrix, damping, q, manifest["epsilon"], manifest["max_iterations"]
        )
        for (index, sample) in (enumerate(samples)):
            difference = np.abs(current[:, index] - actual[:, index])
            require(
                np.all(difference <= 1e-12 + 1e-10 * np.abs(actual[:, index])),
                "Independent full-vector mismatch",
            )
            eligible = np.flatnonzero(query[:, index])
            first = eligible[np.argsort(-actual[eligible, index], kind = "stable")]
            second = eligible[np.argsort(-current[eligible, index], kind = "stable")]
            require(np.array_equal(first, second), "Independent full-query order mismatch")
            error, upper = bound(matrix, damping, q[:, index], actual[:, index])
            for (budget) in (budgets):
                record = diagnostics[(state, sample, budget)]
                require(
                    iterations[index] == int(record["iterations"]),
                    "Stopping iteration mismatch",
                )
                require(
                    bool(stopped[index]) == (record["stopped_by_threshold"] == "TRUE"),
                    "Stopping status mismatch",
                )
                stored = float(record["error_bound"])
                require(
                    abs(error - stored) <= max(1e-18, 1e-4 * stored),
                    "Independent residual-bound mismatch",
                )
                require(
                    abs(upper - float(record["operator_l1_upper"])) <= 1e-14,
                    "Independent operator norm mismatch",
                )
                gap = (
                    float(actual[first[budget - 1], index] - actual[first[budget], index])
                    if (len(first) > budget)
                    else None
                )
                certified = gap is None or gap > 2 * max(error, stored)
                require(
                    certified == (record["certified"] == "TRUE"),
                    "Cutoff certificate status mismatch",
                )
                checks.append(
                    {
                        "state": state,
                        "sample_id": sample,
                        "budget": budget,
                        "iterations": int(iterations[index]),
                        "score_bits_identical": bool(
                            np.array_equal(
                                current[:, index].view(np.uint64),
                                actual[:, index].view(np.uint64),
                            )
                        ),
                        "maximum_score_difference": float(np.max(difference)),
                        "independent_error_bound": error,
                        "stored_error_bound": stored,
                        "certified": certified,
                        "full_query_order_identical": True,
                    }
                )
    return {
        "status": "PASS",
        "vectors": len(samples) * len(manifest["states"]),
        "budgets": budgets,
        "checks": checks,
        "all_score_vectors_bit_identical": all(row["score_bits_identical"] for (row) in (checks)),
    }

def run(config_path):
    config_path = Path(config_path).resolve()
    config = read_json(config_path)
    result = verify_directory(resolve(config_path.parent, config["directory"]))
    output = resolve(config_path.parent, config["output"])
    output.parent.mkdir(parents = True, exist_ok = True)
    write_json(output, result)
    return result
