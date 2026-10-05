evidence_matrix = function(prepared, key) {
    value = prepared$evidence[[key]]
    if (is.list(value)) {
        if (!is.null(value$evaluable)) {
            stopifnot(all(value$evaluable))
        }
        value = value$values
    }
    stopifnot(is.matrix(value), all(is.finite(value)), all(value >= 0))
    return(value)
}

prepared_samples = function(prepared) {
    if (is.data.frame(prepared$models)) {
        return(as.character(prepared$models$sample_id))
    }
    return(as.character(prepared$models))
}

state_grid = function(mode) {
    states = list()
    for (a in 0:1) {
        for (b in 0:1) {
            for (h in 0:1) {
                if (mode == 'factorial' && a != b) {
                    next
                }
                name = if (mode == 'factorial') paste0('G', a, '_H', h) else paste0('P', a, '_D', b, '_H', h)
                states[[name]] = list(a = a, b = b, h = h)
            }
        }
    }
    return(states)
}

write_double = function(values, path) {
    stopifnot(!file.exists(path))
    writeBin(as.double(values), path, size = 8L, endian = 'little')
}

write_precise_csv = function(values, path) {
    for (name in names(values)) {
        if (is.double(values[[name]])) {
            values[[name]] = ifelse(is.na(values[[name]]), 'NA', sprintf('%.17g', values[[name]]))
        }
    }
    write.csv(values, path, row.names = FALSE, na = 'NA')
}

export_numerical_inputs = function(prepared, models, output) {
    writeLines(prepared$genes, file.path(output, 'genes.txt'))
    writeLines(models, file.path(output, 'models.txt'))
    write_double(prepared$mutation, file.path(output, 'mutation.bin'))
    for (key in c('G0', 'G1')) {
        graph = prepared$graphs[[key]]
        transition = as(graph$transition, 'dgCMatrix')
        write_double(transition@x, file.path(output, paste0(key, '_P_values.bin')))
        writeBin(as.integer(transition@i), file.path(output, paste0(key, '_P_rows.bin')), size = 4L, endian = 'little')
        writeBin(as.integer(transition@p), file.path(output, paste0(key, '_P_pointers.bin')), size = 4L, endian = 'little')
        write_double(graph$damping, file.path(output, paste0(key, '_damping.bin')))
        degrees = graph$degree
        if (!is.null(graph$adjacency)) {
            observed = Matrix::colSums(graph$adjacency)
            if (!is.null(degrees)) {
                stopifnot(identical(as.numeric(degrees), as.numeric(observed)))
            }
            degrees = observed
        }
        if (!is.null(degrees)) {
            stopifnot(length(degrees) == length(prepared$genes), all(is.finite(degrees)), all(degrees >= 0))
            write_double(degrees, file.path(output, paste0(key, '_degree.bin')))
        }
    }
    for (key in c('H0', 'H1')) {
        values = evidence_matrix(prepared, key)
        q = vapply(seq_len(ncol(values)), function(index) values[, index] / sum(values[, index]), numeric(nrow(values)))
        write_double(q, file.path(output, paste0(key, '_q.bin')))
    }
}

validate_prepared = function(prepared, specification) {
    models = prepared_samples(prepared)
    genes = as.character(prepared$genes)
    stopifnot(length(genes) == specification$expected_genes, length(models) == specification$expected_models)
    stopifnot(setequal(models, unlist(specification$expected_samples)))
    stopifnot(!anyDuplicated(genes), !anyDuplicated(models), identical(dim(prepared$mutation), c(length(genes), length(models))))
    stopifnot(identical(rownames(prepared$mutation), genes), identical(colnames(prepared$mutation), models))
    stopifnot(all(prepared$mutation %in% c(0, 1)))
    if (!is.null(prepared$seeds)) {
        stopifnot(identical(as.integer(prepared$seeds), c(104729L, 130363L, 155921L, 196613L, 228017L)))
    }
    for (key in c('G0', 'G1')) {
        graph = prepared$graphs[[key]]
        stopifnot(identical(dim(graph$transition), c(length(genes), length(genes))))
        stopifnot(length(graph$damping) == length(genes))
        if (!is.null(dimnames(graph$transition))) {
            stopifnot(identical(rownames(graph$transition), genes), identical(colnames(graph$transition), genes))
        }
    }
    for (key in c('H0', 'H1')) {
        values = evidence_matrix(prepared, key)
        stopifnot(identical(dimnames(values), list(genes, models)))
        stopifnot(all(vapply(seq_len(ncol(values)), function(index) sum(values[, index]) > 0, logical(1))))
    }
    return(models)
}

validate_reuse = function(prepared, original, records, specification) {
    stopifnot(identical(prepared$genes, original$genes), identical(prepared_samples(prepared), prepared_samples(original)))
    stopifnot(identical(prepared$mutation, original$mutation))
    for (key in c('G0', 'G1')) {
        stopifnot(identical(prepared$graphs[[key]]$transition, original$graphs[[key]]$transition))
        stopifnot(identical(prepared$graphs[[key]]$damping, original$graphs[[key]]$damping))
    }
    for (key in c('H0', 'H1')) {
        stopifnot(identical(evidence_matrix(prepared, key), evidence_matrix(original, key)))
    }
    stopifnot(records$epsilon == specification$epsilon, records$max_iterations == specification$max_iterations)
}

run_precision_batch = function(specification) {
    suppressPackageStartupMessages(library(Matrix))
    stopifnot(as.character(getRversion()) == specification$r_version)
    stopifnot(as.character(packageVersion('Matrix')) == specification$matrix_version)
    stopifnot(specification$mode %in% c('factorial', 'components'))
    prepared = readRDS(specification$prepared)
    models = validate_prepared(prepared, specification)
    output = specification$output
    states = state_grid(specification$mode)
    reuse = NULL
    if (specification$mode == 'components') {
        stopifnot(!is.null(specification$reuse_directory))
        reuse = readRDS(file.path(specification$reuse_directory, 'numerical.rds'))
        original = readRDS(file.path(specification$reuse_directory, 'prepared.rds'))
        validate_reuse(prepared, original, reuse, specification)
    }
    saveRDS(prepared, file.path(output, 'prepared.rds'), compress = FALSE, version = 3)
    export_numerical_inputs(prepared, models, output)
    all_records = list()
    metadata = list()
    for (name in names(states)) {
        state = states[[name]]
        transition = prepared$graphs[[paste0('G', state$a)]]$transition
        damping = prepared$graphs[[paste0('G', state$b)]]$damping
        evidence = evidence_matrix(prepared, paste0('H', state$h))
        reused = !is.null(reuse) && state$a == state$b
        source_name = paste0('G', state$a, '_H', state$h)
        vectors = matrix(NA_real_, length(prepared$genes), length(models), dimnames = list(prepared$genes, models))
        for (index in seq_along(models)) {
            sample = models[[index]]
            if (reused) {
                result = reuse$results[[source_name]][[sample]]
                stopifnot(!is.null(result), result$stopped_by_threshold)
            } else {
                result = solve_scores(transition, evidence[, index], damping, specification$epsilon, specification$max_iterations)
            }
            stopifnot(length(result$scores) == length(prepared$genes), all(is.finite(result$scores)))
            vectors[, index] = result$scores
            all_records[[name]][[sample]] = result
            for (budget in unlist(specification$budgets)) {
                cutoff = cutoff_certificate(result$scores, prepared$mutation[, index] == 1, budget, result$certificate$error_bound)
                metadata[[length(metadata) + 1L]] = data.frame(state = name, sample_id = sample, budget = budget,
                    iterations = result$iterations, terminal_update_l2 = result$terminal_update_l2,
                    stopped_by_threshold = result$stopped_by_threshold, result$certificate, cutoff, reused = reused)
            }
        }
        write_double(vectors, file.path(output, paste0(name, '_scores.bin')))
    }
    saveRDS(list(epsilon = specification$epsilon, max_iterations = specification$max_iterations, results = all_records),
        file.path(output, 'numerical.rds'), compress = FALSE, version = 3)
    write_precise_csv(do.call(rbind, metadata), file.path(output, 'cutoff_certificates.csv'))
    jsonlite::write_json(list(schema = 'dawn_precision_scores_v1', mode = specification$mode,
        epsilon = specification$epsilon, max_iterations = specification$max_iterations,
        states = states, models = length(models), genes = length(prepared$genes),
        r_version = as.character(getRversion()), matrix_version = as.character(packageVersion('Matrix'))),
        file.path(output, 'scores.json'), auto_unbox = TRUE, pretty = TRUE, digits = 17)
    return(invisible(NULL))
}
