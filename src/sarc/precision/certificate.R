roundoff_gamma = function(count) {
    unit = .Machine$double.eps / 2
    stopifnot(all(is.finite(count)), all(count >= 0), all(count * unit < 0.5))
    return(count * unit / (1 - count * unit))
}

minimum_positive = function(values) {
    if (any(values > 0)) {
        return(min(values[values > 0]))
    }
    return(1)
}

certificate = function(transition, damping, q, values) {
    size = nrow(transition)
    stopifnot(inherits(transition, 'sparseMatrix'), size > 0L, ncol(transition) == size)
    stopifnot(length(damping) == size, length(q) == size, length(values) == size)
    stopifnot(all(is.finite(transition@x)), all(transition@x >= 0))
    stopifnot(all(is.finite(damping)), all(damping >= 0), all(damping < 1))
    stopifnot(all(is.finite(q)), all(q >= 0), is.finite(sum(q)), sum(q) > 0)
    stopifnot(all(is.finite(values)), all(values >= 0))
    rows = Matrix::rowSums(transition != 0)
    columns = Matrix::colSums(transition != 0)
    propagated = as.numeric(transition %*% values)
    constant = (1 - damping) * q
    residual = damping * propagated + constant - values
    stopifnot(all(is.finite(propagated)), all(is.finite(constant)), all(is.finite(residual)))
    allowance = roundoff_gamma(rows + 3) / (1 - roundoff_gamma(rows + 3)) * abs(damping * propagated)
    allowance = allowance + roundoff_gamma(2) / (1 - roundoff_gamma(2)) * abs(constant)
    allowance = allowance + roundoff_gamma(3) / (1 - roundoff_gamma(3)) * (abs(damping * propagated) + abs(constant) + abs(values))
    rounding = sum(allowance) / (1 - roundoff_gamma(size))
    upper = max(Matrix::colSums(Matrix::Diagonal(x = damping) %*% transition)) / (1 - roundoff_gamma(max(columns) + 3))
    minimum = min(minimum_positive(transition@x) * minimum_positive(values) * minimum_positive(damping),
        minimum_positive(1 - damping) * minimum_positive(q))
    stopifnot(is.finite(upper), upper < 1, minimum > .Machine$double.xmin)
    error = (sum(abs(residual)) / (1 - roundoff_gamma(size)) + rounding) / (1 - upper)
    stopifnot(is.finite(error), error >= 0, is.finite(rounding))
    return(list(error_bound = error, residual_l1 = sum(abs(residual)), operator_l1_upper = upper,
        roundoff_l1 = rounding, minimum_positive_product = minimum, maximum_column_nnz = max(columns)))
}
