import math
import numpy as np
from . import mechanism as native

FEATURE_NAMES = (
    "affine_slope_signed_log1p",
    "affine_intercept_over_old_mean_scale",
    "residual_mean_absolute_over_old_mass",
    "residual_rms_over_old_mass",
    "residual_max_absolute_over_old_mass",
    "residual_concentration",
    "positive_residual_mass_over_old_mass",
    "negative_residual_mass_over_old_mass",
    "sigmoid_oldscore_weighted_residual_over_old_mass",
    "squared_standardised_oldscore_weighted_residual_over_old_mass",
    "candidate_entry_fraction",
    "candidate_exit_fraction",
    "shared_candidate_fraction",
    "degenerate_old_variance",
    "no_shared_candidates",
    "affine_reversal",
)
RTOL = 1e-10
ATOL = 1e-12
parse_pathways = native.parse_pathways
prepare_baseline = native.prepare_baseline

def require(condition, message):
    if (not condition):
        raise ValueError(message)

def affine_features(old_scores, new_scores, old_support, new_support):
    old_scores, new_scores = (
        np.asarray(old_scores, dtype = np.float64),
        np.asarray(new_scores, dtype = np.float64),
    )
    old_support, new_support = (np.asarray(old_support), np.asarray(new_support))
    require(
        old_scores.ndim == new_scores.ndim == old_support.ndim == new_support.ndim == 1,
        "One-dimensional aligned score and support vectors required",
    )
    require(
        old_scores.shape == new_scores.shape == old_support.shape == new_support.shape
        and old_support.dtype == new_support.dtype == np.dtype(bool),
        "Aligned Boolean candidate support required",
    )
    require(
        np.isfinite(old_scores).all()
        and np.isfinite(new_scores).all()
        and np.all(old_scores >= 0)
        and np.all(new_scores >= 0),
        "Finite nonnegative input-only scores required",
    )
    require(
        np.all(old_scores[~old_support] == 0) and np.all(new_scores[~new_support] == 0),
        "Scores outside native support must be zero",
    )
    shared, union = (old_support & new_support, old_support | new_support)
    old_n, new_n, shared_n, union_n = map(
        int, (old_support.sum(), new_support.sum(), shared.sum(), union.sum())
    )
    old_mass = math.fsum(map(float, old_scores[old_support]))
    mass_scale, mean_scale = (max(1.0, old_mass), max(1.0, old_mass / max(1, old_n)))
    require(math.isfinite(mass_scale), "Old input score scale overflow")
    slope = intercept = 0.0
    residual = np.empty(0, dtype = np.float64)
    standardised = np.empty(0, dtype = np.float64)
    degenerate = True
    if (shared_n):
        x, y = (old_scores[shared], new_scores[shared])
        x_scale, y_scale = (max(1.0, float(x.max())), max(1.0, float(y.max())))
        xn, yn = (x / x_scale, y / y_scale)
        x_mean, y_mean = (
            math.fsum(map(float, xn)) / shared_n,
            math.fsum(map(float, yn)) / shared_n,
        )
        xc, yc = (xn - x_mean, yn - y_mean)
        variance = math.fsum((float(value) * float(value) for (value) in (xc))) / shared_n
        degenerate = bool(np.all(x == x[0]))
        if (not degenerate):
            require(variance > 0, "Nonconstant old score variance underflow")
            covariance = (
                math.fsum((float(left) * float(right) for (left, right) in (zip(xc, yc))))
                / shared_n
            )
            coefficient = covariance / variance
            slope = coefficient * (y_scale / x_scale)
            intercept = (y_mean - coefficient * x_mean) * y_scale
            residual = (yc - coefficient * xc) * y_scale
            denominator = max(math.sqrt(variance), 1e-12 / x_scale)
            standardised = xc / denominator
        else:
            intercept = y_mean * y_scale
            residual = yc * y_scale
            standardised = np.zeros(shared_n, dtype = np.float64)
    require(
        math.isfinite(slope)
        and math.isfinite(intercept)
        and np.isfinite(residual).all()
        and np.isfinite(standardised).all(),
        "Nonfinite affine projection",
    )
    raw_residual_max = float(np.abs(residual).max()) if (shared_n) else 0.0
    residual_threshold = ATOL + RTOL * max(
        1.0, float(new_scores[shared].max()) if (shared_n) else 0.0
    )
    canonicalised = raw_residual_max <= residual_threshold
    if (canonicalised):
        residual = np.zeros(shared_n, dtype = np.float64)
    scaled = residual / mass_scale
    absolute = np.abs(scaled)
    largest = float(absolute.max()) if (shared_n) else 0.0
    norm_shape = scaled / largest if (largest) else np.zeros(shared_n)
    squares = math.fsum((float(value) ** 2 for (value) in (norm_shape)))
    concentration = (
        math.fsum((float(value) ** 4 for (value) in (norm_shape))) / squares**2
        if (squares)
        else 0.0
    )
    sigmoid_weights = np.empty(shared_n)
    nonnegative = standardised >= 0
    sigmoid_weights[nonnegative] = 1.0 / (1.0 + np.exp(-standardised[nonnegative]))
    negative_exp = np.exp(standardised[~nonnegative])
    sigmoid_weights[~nonnegative] = negative_exp / (1.0 + negative_exp)
    z_scale = max(1.0, float(np.abs(standardised).max())) if (shared_n) else 1.0
    square_weights = (standardised / z_scale) ** 2

    def weighted(weights):
        denominator = math.fsum(map(float, weights))
        return (
            math.fsum(
                (float(weight) * float(value) for (weight, value) in (zip(weights, scaled)))
            )
            / denominator
            if (denominator)
            else 0.0
        )

    values = (
        math.copysign(math.log1p(abs(slope)), slope) if (slope) else 0.0,
        intercept / mean_scale,
        math.fsum(map(float, absolute)) / shared_n if (shared_n) else 0.0,
        largest * math.sqrt(squares / shared_n) if (shared_n) else 0.0,
        largest,
        concentration,
        math.fsum((float(value) for (value) in (scaled) if (value > 0))),
        -math.fsum((float(value) for (value) in (scaled) if (value < 0))),
        weighted(sigmoid_weights),
        weighted(square_weights),
        int(np.count_nonzero(new_support & ~old_support)) / max(1, union_n),
        int(np.count_nonzero(old_support & ~new_support)) / max(1, union_n),
        shared_n / max(1, union_n),
        float(degenerate),
        float(shared_n == 0),
        float(slope < 0),
    )
    require(
        len(values) == len(FEATURE_NAMES) == 16
        and all((math.isfinite(value) for (value) in (values))),
        "Exactly sixteen finite invariant features required",
    )
    diagnostics = {
        "old_candidate_count": old_n,
        "new_candidate_count": new_n,
        "shared_candidate_count": shared_n,
        "union_candidate_count": union_n,
        "old_score_mass_scale": mass_scale,
        "old_score_mean_scale": mean_scale,
        "affine_slope": slope,
        "affine_intercept": intercept,
        "raw_residual_max": raw_residual_max,
        "residual_zero_threshold": residual_threshold,
        "residual_canonicalised_to_zero": canonicalised,
        "positive_affine_invariant_within_tolerance": bool(
            shared_n
            and (not degenerate)
            and (old_n == new_n == shared_n)
            and (slope > 0)
            and canonicalised
        ),
        "no_gene_order_or_identity_features": True,
    }
    return {
        "status": "AVAILABLE",
        "reason": None,
        "features": dict(zip(FEATURE_NAMES, values)),
        "diagnostics": diagnostics,
    }

def edited_state(baseline, changed_edges):
    require(
        isinstance(baseline, native.Baseline),
        "Bound input-only baseline arithmetic required",
    )
    changed = native.canonical_edges(changed_edges)
    require(
        all(
            (
                baseline.edges[edge] == changed[edge]
                for (edge) in (baseline.edges.keys() & changed.keys())
            )
        ),
        "Shared-edge confidence changes are outside the admitted geometry",
    )
    q, contacts, own_contacts = (
        baseline.state.q.copy(),
        baseline.state.contacts.copy(),
        baseline.state.own_contacts.copy(),
    )
    native.accumulate(
        baseline, baseline.edges.keys() - changed.keys(), -1, q, contacts, own_contacts
    )
    native.accumulate(
        baseline, changed.keys() - baseline.edges.keys(), 1, q, contacts, own_contacts
    )
    return native.state_from_counts(baseline, q, contacts, own_contacts)

def extract_features(baseline, changed_edges):
    after = edited_state(baseline, changed_edges)
    return {
        sample: affine_features(
            baseline.state.scores[index],
            after.scores[index],
            baseline.state.support[index],
            after.support[index],
        )
        for (index, sample) in (enumerate(baseline.population))
    }
