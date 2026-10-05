dawn_transition = function(adjacency) {
    stopifnot(all(adjacency@x %in% c(0, 1)))
    degrees = Matrix::colSums(adjacency)
    transition = adjacency %*% Matrix::Diagonal(x = 1 / (degrees + 1e-16))
    dimnames(transition) = dimnames(adjacency)
    return(transition)
}

dawn_damping = function(adjacency, mu = 3) {
    degrees = Matrix::colSums(adjacency)
    return(degrees / (degrees + mu))
}

dawn_scores = function(transition, expression, mutation, damping, sample, epsilon, max_iterations) {
    solved = solve_scores(transition, expression, damping, epsilon, max_iterations)
    ranking = matrix(solved$scores, length(expression), 1L, dimnames = list(rownames(transition), NULL))
    output = data.frame(Rank = ranking, PercentRank = 100 * rank(ranking) / (length(ranking) + 1), isMutated = mutation)
    if (sum(output$isMutated) > 0L) {
        mutated = output[output$isMutated == 1, ]
        mutated = mutated[, -3L]
        mutated = data.frame(Gene = rownames(mutated), Patient = sample, mutated)
        rownames(mutated) = NULL
    }
    else {
        mutated = logical(0)
    }
    return(list(summaryOutput = output, mutatedRanks = mutated, convergenceIterations = solved$iterations,
        terminal_update_norm = solved$terminal_update_l2,
        termination_reason = if (solved$stopped_by_threshold) "UPDATE_THRESHOLD" else "MAXIMUM_ITERATIONS_EXHAUSTED",
        fixed_point_error_certified = is.finite(solved$certificate$error_bound), certificate = solved$certificate,
        cutoff = cutoff_certificate(solved$scores, as.logical(mutation), 10L, solved$certificate$error_bound)))
}

project_inherited = function(original, native_genes, genes, models) {
    stopifnot(is.list(original), all(c("values", "evaluable", "sample_status") %in% names(original)))
    stopifnot(is.matrix(original$values), is.numeric(original$values), identical(rownames(original$values),
        native_genes), identical(colnames(original$values), models))
    stopifnot(!anyDuplicated(native_genes), !anyDuplicated(genes), all(genes %in% native_genes), identical(names(original$evaluable),
        models))
    stopifnot(is.logical(original$evaluable), !anyNA(original$evaluable), length(original$sample_status) ==
        length(models))
    if (!is.null(names(original$sample_status))) {
        stopifnot(identical(names(original$sample_status), models))
    }
    projected = original$values[genes, models, drop = FALSE]
    consumed = setNames(vector("list", length(models)), models)
    ready = setNames(rep(FALSE, length(models)), models)
    reasons = setNames(rep("", length(models)), models)
    for (index in seq_along(models)) {
        available = original$evaluable[[index]]
        stopifnot(!available || all(is.finite(original$values[, index])))
        values = projected[, index]
        finite = all(is.finite(values))
        nonnegative = finite && all(values >= 0)
        mass = if (finite)
            sum(values)
        else NA_real_
        ready[[index]] = available && finite && nonnegative && is.finite(mass) && mass > 0
        reasons[[index]] = if (!available)
            "ORIGINAL_NORMALISATION_UNAVAILABLE"
        else if (!finite)
            "PROJECTED_NONFINITE"
        else if (!nonnegative)
            "PROJECTED_NEGATIVE"
        else if (!is.finite(mass))
            "NONFINITE_L1_MASS"
        else if (mass <= 0)
            "NONPOSITIVE_L1_MASS"
        else "INPUT_READY"
        if (ready[[index]]) {
            consumed[[index]] = values / sum(values)
            stopifnot(all(is.finite(consumed[[index]])))
        }
    }
    return(list(values = projected, consumed = consumed, evaluable = ready, sample_status = reasons))
}

prepare_native_dawn = function(specification) {
    tumour = read_native_dawn_matrix(specification$tumour)
    normal = read_native_dawn_matrix(specification$normal)
    mutation = read_native_dawn_matrix(specification$mutation)
    models = unlist(specification$models, use.names = FALSE)
    stopifnot(ncol(tumour) == specification$expected_models, ncol(normal) == specification$expected_normals,
        identical(colnames(tumour), models))
    stopifnot(identical(dimnames(tumour), dimnames(mutation)), identical(rownames(tumour), rownames(normal)),
        all(mutation %in% c(0, 1)))
    names = read.csv(specification$normalisation_names, check.names = FALSE, stringsAsFactors = FALSE)
    tumour_names = names[names$role == "tumour", , drop = FALSE]
    normal_names = names[names$role == "normal", , drop = FALSE]
    stopifnot(identical(tumour_names$canonical, colnames(tumour)), identical(normal_names$canonical,
        colnames(normal)))
    original = new.env(parent = .GlobalEnv)
    base::sys.source(specification$normalise_source, envir = original)
    normalised = normalise_dawn_input(tumour, normal, tumour_names$original, normal_names$original, original$DawnNormalize)
    saveRDS(normalised, file.path(specification$output, "normalised_input.rds"), version = 3, compress = FALSE)
    adjacency = aligned_dawn_graph(read_edges(specification$network), rownames(tumour))
    graph = list(adjacency = adjacency, transition = dawn_transition(adjacency), damping = dawn_damping(adjacency,
        specification$mu))
    return(list(genes = rownames(tumour), models = models, mutation = mutation, graphs = list(G0 = graph),
        evidence = list(H0 = normalised), cells = "G0_H0"))
}

prepare_factorial_dawn = function(specification) {
    genes = readLines(specification$common_genes, warn = FALSE)
    models = unlist(specification$models, use.names = FALSE)
    mutation = read_native_dawn_matrix(specification$mutation)
    stopifnot(length(genes) == specification$expected_genes, length(models) == specification$expected_models,
        !anyDuplicated(genes), !anyDuplicated(models))
    stopifnot(identical(rownames(mutation), genes), identical(colnames(mutation), models), all(mutation %in%
        c(0, 1)))
    graphs = lapply(c("G0", "G1"), function(key) {
        adjacency = aligned_dawn_graph(read_edges(specification[[paste0("edges_", key)]]), genes)
        return(list(adjacency = adjacency, transition = dawn_transition(adjacency), damping = dawn_damping(adjacency,
            specification$mu)))
    })
    names(graphs) = c("G0", "G1")
    evidence = lapply(c("H0", "H1"), function(key) {
        original = readRDS(specification[[paste0("inherited_", key)]])
        native_genes = readLines(specification[[paste0("native_genes_", key)]], warn = FALSE)
        return(project_inherited(original, native_genes, genes, models))
    })
    names(evidence) = c("H0", "H1")
    return(list(genes = genes, models = models, mutation = mutation, graphs = graphs, evidence = evidence,
        cells = c("G0_H0", "G1_H0", "G0_H1", "G1_H1")))
}

dawn_cell_identity = function(prepared, cell, sample) {
    pieces = strsplit(cell, "_", fixed = TRUE)[[1L]]
    graph = pieces[[1L]]
    evidence = pieces[[2L]]
    if (length(prepared$graphs) == 2L && identical(prepared$graphs$G0, prepared$graphs$G1)) {
        graph = "G0"
    }
    if (length(prepared$evidence) == 2L && prepared$evidence$H0$evaluable[[sample]] && prepared$evidence$H1$evaluable[[sample]] &&
        identical(prepared$evidence$H0$consumed[[sample]], prepared$evidence$H1$consumed[[sample]])) {
        evidence = "H0"
    }
    return(paste(graph, evidence, sep = "_"))
}

run_dawn = function(specification) {
    suppressPackageStartupMessages(library(Matrix))
    stopifnot(as.character(packageVersion("Matrix")) == "1.7.4")
    if (is.null(specification$mu)) {
        specification$mu = 3
    }
    stopifnot(specification$mu %in% c(1, 3, 10))
    stopifnot(!is.null(specification$epsilon), !is.null(specification$max_iterations))
    stopifnot(specification$epsilon %in% c(1e-4, 1e-6, 1e-8, 1e-12))
    stopifnot(specification$max_iterations == if (specification$epsilon == 1e-4) 100L else 10000L)
    parameter_id = paste0("mu_", specification$mu)
    prepared = if (specification$mode == "native")
        prepare_native_dawn(specification)
    else prepare_factorial_dawn(specification)
    saveRDS(prepared, file.path(specification$output, "prepared_input.rds"), version = 3, compress = FALSE)
    rank_fields = c("context", "sample_id", "method", "network_id", "evidence_id", "cell", "parameter_id",
        "gene", "rank_position", "score", "native_percent_rank")
    full_fields = c("cell", "sample_id", "gene", "score", "native_percent_rank", "is_mutated")
    ranked_path = file.path(specification$output, "ranked_scores.csv")
    full_path = file.path(specification$output, "all_gene_scores.csv")
    write.csv(empty_frame(rank_fields), ranked_path, row.names = FALSE)
    write.csv(empty_frame(full_fields), full_path, row.names = FALSE)
    native = file.path(specification$output, "native_results")
    stopifnot(dir.create(native))
    cache = new.env(parent = emptyenv())
    statuses = list()
    aliases = list()
    for (cell in prepared$cells) {
        pieces = strsplit(cell, "_", fixed = TRUE)[[1L]]
        graph = prepared$graphs[[pieces[[1L]]]]
        evidence = prepared$evidence[[pieces[[2L]]]]
        for (index in seq_along(prepared$models)) {
            sample = prepared$models[[index]]
            canonical = dawn_cell_identity(prepared, cell, sample)
            key = paste(canonical, sample, sep = "/")
            status = "UNAVAILABLE"
            reason = as.character(evidence$sample_status[[index]])
            count = 0L
            iteration = NA_integer_
            norm = NA_real_
            termination = "NOT_EXECUTED_UNAVAILABLE"
            attempted = FALSE
            certified = FALSE
            error_bound = NA_real_
            cutoff_certified = FALSE
            if (evidence$evaluable[[sample]]) {
                if (!exists(key, envir = cache, inherits = FALSE)) {
                  stopifnot(canonical == cell)
                  result = tryCatch(dawn_scores(graph$transition, evidence$values[, sample], prepared$mutation[,
                    sample], graph$damping, sample, specification$epsilon, specification$max_iterations), error = function(error) error)
                  assign(key, result, envir = cache)
                  attempted = TRUE
                  saveRDS(result, file.path(native, paste0(cell, "_", sprintf("%03d", index), ".rds")),
                    version = 3)
                }
                result = get(key, envir = cache, inherits = FALSE)
                if (inherits(result, "error")) {
                  status = "FAILED"
                  reason = conditionMessage(result)
                  termination = "NATIVE_ERROR"
                }
                else if (!all(is.finite(result$summaryOutput$Rank))) {
                  reason = "Native DawnRank returned non-finite scores"
                }
                else {
                  iteration = result$convergenceIterations
                  norm = result$terminal_update_norm
                  termination = result$termination_reason
                  certified = result$fixed_point_error_certified
                  error_bound = result$certificate$error_bound
                  cutoff_certified = result$cutoff$certified
                  complete = data.frame(cell = cell, sample_id = sample, gene = prepared$genes, score = sprintf("%.17g",
                    result$summaryOutput$Rank), native_percent_rank = sprintf("%.17g", result$summaryOutput$PercentRank),
                    is_mutated = prepared$mutation[, sample], stringsAsFactors = FALSE)
                  append_csv(complete, full_path)
                  candidates = result$mutatedRanks
                  count = if (is.data.frame(candidates))
                    nrow(candidates)
                  else 0L
                  status = if (count > 0L)
                    "SUCCESS"
                  else "SUCCESS_EMPTY"
                  reason = if (count > 0L)
                    ""
                  else "Original mutation-filtered output has no candidate"
                  if (count > 0L) {
                    candidates = candidates[order(-candidates$PercentRank, method = "radix"), , drop = FALSE]
                    network_id = if (specification$mode == "native")
                      specification$network_id
                    else pieces[[1L]]
                    rows = data.frame(context = specification$context, sample_id = sample, method = "DawnRank",
                      network_id = network_id, evidence_id = pieces[[2L]], cell = cell, parameter_id = parameter_id,
                      gene = candidates$Gene, rank_position = seq_len(count), score = sprintf("%.17g",
                        candidates$Rank), native_percent_rank = sprintf("%.17g", candidates$PercentRank),
                      stringsAsFactors = FALSE)
                    append_csv(rows, ranked_path)
                  }
                }
            }
            network_id = if (specification$mode == "native")
                specification$network_id
            else pieces[[1L]]
            statuses[[length(statuses) + 1L]] = data.frame(context = specification$context, sample_id = sample,
                method = "DawnRank", network_id = network_id, evidence_id = pieces[[2L]], cell = cell,
                parameter_id = parameter_id, status = status, reason = reason, candidate_count = count,
                canonical_cell = canonical, native_required = attempted, iterations = iteration, terminal_update_norm = norm,
                termination_reason = termination, fixed_point_error_certified = certified,
                fixed_point_error_bound = error_bound, cutoff_certified = cutoff_certified, stringsAsFactors = FALSE)
            aliases[[length(aliases) + 1L]] = data.frame(cell = cell, sample_id = sample, master_seed = unlist(specification$seeds),
                canonical_cell = canonical, status = status, seed_interpretation = "DETERMINISTIC_REFERENCE",
                stringsAsFactors = FALSE)
        }
    }
    statuses = do.call(rbind, statuses)
    write.csv(statuses, file.path(specification$output, "sample_status.csv"), row.names = FALSE)
    write.csv(do.call(rbind, aliases), file.path(specification$output, "seed_aliases.csv"), row.names = FALSE)
    write_json(list(method = "DawnRank", mode = specification$mode, model_states = nrow(statuses), physical_attempts = sum(statuses$native_required),
        successes = sum(statuses$status == "SUCCESS"), empty = sum(statuses$status == "SUCCESS_EMPTY"),
        unavailable = sum(statuses$status == "UNAVAILABLE"), failed = sum(statuses$status == "FAILED")),
        file.path(specification$output, "completion.json"))
    return(!any(statuses$status == "FAILED"))
}
