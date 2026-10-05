import itertools
import math
import numpy as np

COHORTS = {"COAD_CCLE": 36, "LUAD_CCLE": 36, "COAD_TCGA": 396}
TRANSITIONS = ("11_0_to_11_5", "11_5_to_12_0")
DIRECTIONS = ("insert", "retract")
CLASSES = ("PB", "PW", "DB", "DW", "PB_PW")
SEEDS = (271828, 314159, 173205)
AXES = (
    ("removed", "lower"),
    ("removed", "higher"),
    ("added", "lower"),
    ("added", "higher"),
)

def require(condition, message):
    if (not condition):
        raise ValueError(message)

def finite(value):
    require(
        type(value) in (int, float) and math.isfinite(value),
        "Finite numeric metadata required",
    )
    return float(value)

def residuals(balance, components = None):
    values = []
    for (direction, role) in (AXES):
        item = balance["continuous_log_degree_balance"][direction][role]
        require(
            item["state"] == "DEFINED" and item["pass"] is True and (item["n"] > 0),
            "Missing/failed matching covariate",
        )
        values.append(
            finite(item["real_mean_log_degree"])
            - finite(item["control_mean_log_degree"])
        )
    for (direction) in (("removed", "added")):
        entries = balance["continuous_log_degree_balance"][direction]
        require(
            entries["pooled"]["n"]
            == 2 * entries["lower"]["n"]
            == 2 * entries["higher"]["n"],
            "Endpoint multiplicity differs",
        )
        for (side) in (("real", "control")):
            field = side + "_mean_log_degree"
            pooled = finite(entries["pooled"][field])
            expected = (
                finite(entries["lower"][field]) + finite(entries["higher"][field])
            ) / 2
            require(
                abs(pooled - expected)
                <= 1e-12 + 1e-10 * max(abs(pooled), abs(expected)),
                "Pooled algebra differs",
            )
    if (balance.get("exact_degree_cost_confidence_strata") is not True):
        require(
            components is not None
            and len(components) == 2
            and all(
                (
                    item.get("exact_degree_cost_confidence_strata") is True
                    for (item) in (components)
                )
            ),
            "Exact matching strata not accepted",
        )
        for (direction) in (("removed", "added")):
            require(
                balance[direction] == sum((item[direction] for (item) in (components))),
                "PB_PW disjoint component count differs",
            )
            for (role) in (("lower", "higher", "pooled")):
                combined = balance["continuous_log_degree_balance"][direction][role]
                parts = [
                    item["continuous_log_degree_balance"][direction][role]
                    for (item) in (components)
                ]
                require(
                    combined["n"] == sum((item["n"] for (item) in (parts)))
                    and all(
                        (
                            item["pass"] is True and item["state"] == "DEFINED"
                            for (item) in (parts)
                        )
                    ),
                    "PB_PW component support differs",
                )
                for (side) in (("real", "control")):
                    key = side + "_mean_log_degree"
                    expected = (
                        sum((finite(item[key]) * item["n"] for (item) in (parts)))
                        / combined["n"]
                    )
                    actual = finite(combined[key])
                    require(
                        abs(actual - expected)
                        <= 1e-12 + 1e-10 * max(abs(actual), abs(expected)),
                        "PB_PW component weighted mean differs",
                    )
    return values

def svd_rank(matrix):
    matrix = np.asarray(matrix, dtype = float)
    values = np.linalg.svd(matrix, compute_uv = False)
    threshold = (
        max(matrix.shape) * np.finfo(float).eps * values[0] if (len(values)) else 0.0
    )
    rank = int(np.sum(values > threshold))
    return {
        "rows": matrix.shape[0],
        "columns": matrix.shape[1],
        "singular_values": values.tolist(),
        "tolerance": float(threshold),
        "rank": rank,
        "condition_number": float(values[0] / values[-1])
        if (len(values) and values[-1] > 0)
        else None,
    }

def hull_zero(points):
    points = np.asarray(points, dtype = float)
    require(
        points.ndim == 2 and points.shape[1] == 4 and np.isfinite(points).all(),
        "Four finite covariates required",
    )
    best = None
    for (size) in (range(1, len(points) + 1)):
        for (subset) in (itertools.combinations(range(len(points)), size)):
            selected = points[list(subset)]
            if (size == 1):
                weights = np.ones(1)
            else:
                differences = (selected[:-1] - selected[-1]).T
                first = np.linalg.lstsq(differences, -selected[-1], rcond = None)[0]
                weights = np.concatenate((first, [1.0 - first.sum()]))
            if (weights.min() < -1e-12):
                continue
            weights = np.maximum(weights, 0.0)
            weights /= weights.sum()
            distance = float(np.linalg.norm(weights @ selected))
            if (best is None or distance < best["distance"]):
                full = np.zeros(len(points))
                full[list(subset)] = weights
                best = {"distance": distance, "weights": full.tolist()}
    require(best is not None, "Convex hull candidate missing")
    tolerance = 1e-12 + 1e-10 * max(1.0, float(np.abs(points).max()))
    return best | {
        "tolerance": tolerance,
        "inside_convex_hull": best["distance"] <= tolerance,
        "coordinatewise_zero_inside": bool(
            ((points.min(axis = 0) <= 0) & (points.max(axis = 0) >= 0)).all()
        ),
    }

def design(rows):
    require(bool(rows), "Empty design")
    keys = [row["stratum"] for (row) in (rows)]
    strata = list(dict.fromkeys(keys))
    x = np.asarray([row["covariates"] for (row) in (rows)], dtype = float)
    require(
        x.shape == (len(rows), 4) and np.isfinite(x).all(), "Invalid covariate matrix"
    )
    centered = x.copy()
    onehot = np.zeros((len(rows), len(strata)))
    support = []
    for (column, key) in (enumerate(strata)):
        mask = np.asarray([value == key for (value) in (keys)])
        require(mask.sum() == 3, "Exactly three controls per stratum required")
        require(
            {rows[i]["control_seed"] for (i) in (np.flatnonzero(mask))} == set(SEEDS),
            "Complete compact controls required",
        )
        centered[mask] -= x[mask].mean(axis = 0)
        onehot[mask, column] = 1.0
        support.append(
            {
                "stratum": key,
                "mean_covariates": x[mask].mean(axis = 0).tolist(),
                "zero_support": hull_zero(x[mask]),
            }
        )
    scales = np.sqrt(np.mean(centered**2, axis = 0))
    nonzero = bool((scales > 0).all())
    safe = np.where(scales > 0, scales, 1.0)
    within = svd_rank(centered / safe)
    full = svd_rank(np.column_stack((onehot, x / safe)))
    identified = (
        nonzero
        and within["rank"] == 4
        and (full["rank"] == len(strata) + 4)
        and (within["condition_number"] is not None)
        and (within["condition_number"] <= 100000000.0)
        and (full["condition_number"] is not None)
        and (full["condition_number"] <= 100000000.0)
    )
    return {
        "status": "INPUT_ONLY_LINEAR_PROJECTION_IDENTIFIED"
        if (identified)
        else "NUMERICALLY_NOT_IDENTIFIED",
        "row_count": len(rows),
        "stratum_count": len(strata),
        "within_rms_scales": scales.tolist(),
        "within_svd": within,
        "full_svd": full,
        "strata": support,
        "outside_zero_hull_count": sum(
            (not item["zero_support"]["inside_convex_hull"] for (item) in (support))
        ),
        "distinct_covariate_vectors": len({tuple(row) for (row) in (x)}),
        "inference": "DESCRIPTIVE_NO_ADJUSTED_INFERENCE",
        "graph_independence_claim": False,
    }
