import itertools
from collections import Counter
from fractions import Fraction
import numpy as np

COHORTS = {"COAD_CCLE": 36, "LUAD_CCLE": 36, "COAD_TCGA": 396}
METHODS = ("DawnRank", "PersonaDrive", "PRODIGY")
TRANSITIONS = ("11_0_to_11_5", "11_5_to_12_0")
DIRECTIONS = ("insert", "retract")
CLASSES = ("PB", "PW", "DB", "DW", "PB_PW")
CONTROLS = (271828, 314159, 173205)
ENDPOINTS = tuple(
    itertools.product(
        ("common", "native"),
        (
            "NCG6_primary_all",
            "NCG6_known_subgroup",
            "NCG6_candidate_subgroup",
            "CancerMine2019_secondary",
        ),
        (10, 1, 5, 20),
    )
)
MASTERS = [104729, 130363, 155921, 196613, 228017]

def require(condition, message):
    if (not condition):
        raise ValueError(message)

def exact(value):
    if (value is None):
        return None
    require(
        isinstance(value, dict) and set(value) == {"numerator", "denominator"},
        "Exact inherited rational point required",
    )
    require(
        type(value["numerator"]) is int
        and type(value["denominator"]) is int
        and (value["denominator"] > 0),
        "Invalid rational point",
    )
    return Fraction(value["numerator"], value["denominator"])

def encode(value):
    if (isinstance(value, Fraction)):
        return {"numerator": value.numerator, "denominator": value.denominator}
    raise TypeError(type(value).__name__)

def identifier(cohort, method, endpoint):
    return "::".join((cohort, method, *map(str, endpoint)))

def compute(rows, admitted_design, originals, cohort, method, endpoint):
    require(len(rows) == len(originals) == 60, "Complete60 graph rows required")
    expected_strata = [
        "::".join((cohort, transition, direction, class_id))
        for (transition) in (TRANSITIONS)
        for (direction) in (DIRECTIONS)
        for (class_id) in (CLASSES)
    ]
    require(
        list(dict.fromkeys((row["stratum"] for (row) in (rows)))) == expected_strata,
        "Exact20 stratum order required",
    )
    require(
        [item["stratum"] for (item) in (admitted_design["strata"])] == expected_strata,
        "Admitted support correspondence differs",
    )
    x = np.asarray([row["covariates"] for (row) in (rows)], dtype = float)
    require(
        x.shape == (60, 4) and np.isfinite(x).all(),
        "Finite complete input covariates required",
    )
    points = [exact(row["original_signed_scalar"]["value"]) for (row) in (originals)]
    missing = [
        rows[index]["network_id"]
        for ((index, value)) in (enumerate(points))
        if (value is None)
    ]
    reason = (
        "INCOMPLETE_FULL_POPULATION_OUTCOMES"
        if (missing)
        else "NUMERICALLY_NOT_IDENTIFIED"
        if (admitted_design["status"] != "INPUT_ONLY_LINEAR_PROJECTION_IDENTIFIED")
        else None
    )
    (means_x, means_y, indices) = ([], [], [])
    for (stratum) in (expected_strata):
        positions = [
            index for ((index, row)) in (enumerate(rows)) if (row["stratum"] == stratum)
        ]
        require(
            len(positions) == 3
            and [rows[index]["control_seed"] for (index) in (positions)]
            == list(CONTROLS),
            "Exact three-control order required",
        )
        indices.append(positions)
        means_x.append(x[positions].mean(axis = 0))
        values = [points[index] for (index) in (positions)]
        means_y.append(
            sum(values, Fraction(0)) / 3
            if (all((value is not None for (value) in (values))))
            else None
        )
    coefficients = None
    residual_norm = None
    if (reason is None):
        scales = np.asarray(admitted_design["within_rms_scales"], dtype = float)
        require(
            scales.shape == (4,) and np.isfinite(scales).all() and (scales > 0).all(),
            "Admitted nonzero input scales required",
        )
        centered_x = x.copy()
        y = np.asarray([float(value) for (value) in (points)])
        centered_y = y.copy()
        for (positions, mean_x, mean_y) in (zip(indices, means_x, means_y)):
            centered_x[positions] -= mean_x
            centered_y[positions] -= float(mean_y)
        scaled_beta = np.linalg.lstsq(centered_x / scales, centered_y, rcond = None)[0]
        coefficients = scaled_beta / scales
        require(
            np.isfinite(coefficients).all(), "Finite projection coefficients required"
        )
        residual_norm = float(np.linalg.norm(centered_y - centered_x @ coefficients))
        coefficients = coefficients.tolist()
    strata = []
    for (index, stratum) in (enumerate(expected_strata)):
        support = admitted_design["strata"][index]["zero_support"]
        adjusted = (
            float(means_y[index])
            - float(np.asarray(means_x[index]) @ np.asarray(coefficients))
            if (coefficients is not None)
            else None
        )
        strata.append(
            {
                "stratum": stratum,
                "original_control_network_ids": [
                    rows[position]["network_id"] for (position) in (indices[index])
                ],
                "unadjusted_control_mean": means_y[index],
                "adjusted_at_zero": adjusted,
                "adjustment": adjusted - float(means_y[index])
                if (adjusted is not None)
                else None,
                "status": "DESCRIPTIVE_COMMON_SLOPE_EXTRAPOLATION"
                if (adjusted is not None and (not support["inside_convex_hull"]))
                else "DESCRIPTIVE_WITHIN_HULL_LINEAR_PROJECTION"
                if (adjusted is not None)
                else "UNAVAILABLE_" + reason,
                "zero_support": support,
                "adjusted_interval": None,
                "inference": "DESCRIPTIVE_NO_ADJUSTED_INFERENCE",
                "original_intervals_are_references_only": True,
                "common_slope_assumption": True,
            }
        )
    return {
        "projection_id": identifier(cohort, method, endpoint),
        "cohort": cohort,
        "method": method,
        "endpoint": list(endpoint),
        "status": "COMPLETE_DESCRIPTIVE_RESIDUAL_PROJECTION"
        if (reason is None)
        else "UNAVAILABLE_" + reason,
        "row_count": 60,
        "stratum_count": 20,
        "original_contrast_ids": [row["contrast_id"] for (row) in (originals)],
        "missing_outcome_network_ids": missing,
        "coefficients_raw_logdegree": coefficients,
        "centered_residual_l2": residual_norm,
        "strata": strata,
        "predictor_fit": False,
        "graph_independence_claim": False,
        "adjusted_inference": False,
    }

def project_all(assessment, scalars):
    require(
        len(scalars) == 17280
        and len({row["contrast_id"] for (row) in (scalars)}) == 17280,
        "Complete17280 unique original scalar references required",
    )
    lookup = {}
    for (row) in (scalars):
        key = (
            row["source_slice"].split("::")[0],
            row["method"],
            tuple(row["endpoint"]),
            row["network_id"],
        )
        require(key not in lookup, "Duplicate scalar identity")
        lookup[key] = row
    projections = []
    for (cohort) in (COHORTS):
        rows = [row for (row) in (assessment["rows"]) if (row["cohort"] == cohort)]
        for (method, endpoint) in (itertools.product(METHODS, ENDPOINTS)):
            originals = [
                lookup[cohort, method, endpoint, row["network_id"]] for (row) in (rows)
            ]
            projections.append(
                compute(
                    rows,
                    assessment["cohort_designs"][cohort],
                    originals,
                    cohort,
                    method,
                    endpoint,
                )
            )
    require(
        len(projections) == 288
        and sum((len(row["strata"]) for (row) in (projections))) == 5760,
        "Complete288/5760 accounting required",
    )
    return {
        "status": "COMPLETE_RESIDUAL_SENSITIVITY",
        "projections": projections,
        "projection_count": 288,
        "stratum_point_count": 5760,
        "original_scalar_count": 17280,
        "projection_status_counts": dict(
            Counter((row["status"] for (row) in (projections)))
        ),
        "stratum_status_counts": dict(
            Counter(
                (
                    point["status"]
                    for (row) in (projections)
                    for (point) in (row["strata"])
                )
            )
        ),
        "complete_populations": COHORTS,
        "algorithm_seeds": MASTERS,
        "compact_graph_control_seeds": list(CONTROLS),
        "ALL_has_no_matched_control_estimand": True,
        "legacy_control_seeds_separate": [161803, 141421],
        "no_bootstrap_or_predictor_fit": True,
        "original_intervals_reused_as_adjusted": False,
        "common_slope_linear_projection_not_causal_identification": True,
    }
