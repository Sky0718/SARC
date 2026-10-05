solve_scores = function(transition, evidence, damping, epsilon = 1e-12, max_iterations = 10000L) {
    size = nrow(transition)
    stopifnot(inherits(transition, 'sparseMatrix'), size > 0L, ncol(transition) == size)
    stopifnot(length(evidence) == size, length(damping) == size)
    stopifnot(all(is.finite(transition@x)), all(transition@x >= 0))
    stopifnot(all(is.finite(evidence)), all(evidence >= 0), is.finite(sum(evidence)), sum(evidence) > 0)
    stopifnot(all(is.finite(damping)), all(damping >= 0), all(damping < 1))
    stopifnot(length(epsilon) == 1L, is.finite(epsilon), epsilon > 0)
    stopifnot(length(max_iterations) == 1L, max_iterations >= 1L, max_iterations == as.integer(max_iterations))
    q = evidence / sum(evidence)
    ranking = matrix(q, length(evidence), 1L)
    constant = q * (1 - damping)
    for (iteration in seq_len(max_iterations)) {
        updated = as.matrix((damping * (transition %*% ranking)) + constant)
        stopifnot(all(is.finite(updated)), all(updated >= 0))
        magnitude = sqrt(sum((updated - ranking) ^ 2))
        ranking = updated
        if (magnitude < epsilon) {
            break
        }
    }
    values = as.numeric(ranking)
    return(list(scores = values, iterations = iteration, terminal_update_l2 = magnitude,
        stopped_by_threshold = magnitude < epsilon, certificate = certificate(transition, damping, q, values)))
}

cutoff_certificate = function(scores, eligible, budget, error_bound) {
    stopifnot(length(scores) == length(eligible), is.logical(eligible), !anyNA(eligible))
    stopifnot(all(is.finite(scores)), length(budget) == 1L, budget >= 0L, budget == as.integer(budget))
    stopifnot(is.finite(error_bound), error_bound >= 0)
    ordered = which(eligible)
    ordered = ordered[order(-scores[ordered], method = 'radix')]
    complete = length(ordered) <= budget
    empty_budget = budget == 0L
    gap = if (complete || empty_budget) NA_real_ else scores[ordered[[budget]]] - scores[ordered[[budget + 1L]]]
    tied = if (complete || empty_budget) 0L else sum(scores[ordered] == scores[ordered[[budget]]])
    return(list(query_size = length(ordered), complete_menu = complete, empty_budget = empty_budget,
        cutoff_gap = gap, cutoff_tied = tied, certified = complete || empty_budget || gap > 2 * error_bound))
}
