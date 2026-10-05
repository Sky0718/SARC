read_explicit_evidence = function(path, common_genes) {
    value = jsonlite::fromJSON(path, simplifyVector = FALSE)
    stopifnot(identical(value$schema, "cdf_explicit_evidence_v1"), identical(value$status, "AVAILABLE"))
    stopifnot(isTRUE(value$non_null_empty), isTRUE(value$restriction_preserves_original_order_and_bits))
    genes = if (length(value$genes) == 0L)
        character()
    else unlist(value$genes, use.names = FALSE)
    stopifnot(is.character(genes), !anyNA(genes), !anyDuplicated(genes), all(nchar(genes) > 0L), all(genes %in%
        common_genes))
    bits = value$binary64_le
    stopifnot(is.character(bits), length(bits) == 1L, nchar(bits) == 16L * length(genes), value$size ==
        length(genes), grepl("^[0-9a-f]*$", bits))
    if (nchar(bits) == 0L) {
        differences = setNames(numeric(), character())
    }
    else {
        starts = seq.int(1L, nchar(bits), by = 2L)
        bytes = as.raw(strtoi(substring(bits, starts, starts + 1L), base = 16L))
        differences = setNames(readBin(bytes, what = "double", n = length(genes), size = 8L, endian = "little"),
            genes)
    }
    stopifnot(!is.null(differences), is.double(differences), all(is.finite(differences)), all(differences >=
        0), identical(names(differences), genes))
    return(differences)
}

binary64 = function(values) {
    return(paste(sprintf("%02x", as.integer(writeBin(as.double(values), raw(), size = 8L, endian = "little"))),
        collapse = ""))
}

prodigy_environment = function(specification, state) {
    execution = new.env(parent = .GlobalEnv)
    execution$profile_state = state
    stopifnot(length(specification$original_sources) == 8L)
    for (path in specification$original_sources) {
        base::sys.source(path, envir = execution)
    }
    original_enrichment = execution$get_enriched_pathways
    execution$get_enriched_pathways = function(...) {
        value = original_enrichment(...)
        state$enriched = names(value)
        state$enrichment_completed = TRUE
        return(value)
    }
    if (specification$mode == "factorial") {
        execution$get_DEGs = function(...) stop("Differential-expression refitting is not part of fixed inherited evidence")
        execution$get_diff_expressed_genes = function(...) stop("Differential-expression refitting is not part of fixed inherited evidence")
    }
    execution$normalmixEM = function(...) {
        entry = list(arguments = list(...), rng_before = .Random.seed)
        value = tryCatch(mixtools::normalmixEM(...), error = function(error) error)
        entry$rng_after = .Random.seed
        entry$status = if (inherits(value, "error"))
            "ERROR"
        else "RETURNED"
        entry$result = if (inherits(value, "error"))
            NULL
        else value
        entry$error = if (inherits(value, "error"))
            conditionMessage(value)
        else NULL
        state$mixtures[[length(state$mixtures) + 1L]] = entry
        if (inherits(value, "error")) {
            stop(value)
        }
        return(value)
    }
    invisible(trace("PRODIGY", , quote({
        profile_state$axes = list(mutation_genes = mutated_genes, evaluated_pathways = if (exists("Influence_matrix",
            inherits = FALSE)) rownames(Influence_matrix)[-1L] else NULL, influence_before_drop = if (exists("Influence_matrix",
            inherits = FALSE)) Influence_matrix else NULL, returned = returnValue())
    }), , FALSE, where = execution))
    invisible(trace("analyze_PRODIGY_results", , quote({
        profile_state$aggregation = if (exists("ranking", inherits = FALSE)) ranking else NULL
    }), , FALSE, where = execution))
    return(execution)
}

prodigy_patient = function(specification, output, state) {
    execution = prodigy_environment(specification, state)
    counts = read_numeric_matrix(specification$counts)
    mutations = read_numeric_matrix(specification$mutations)
    network = original_prodigy_network(read_edges(specification$network, parser = "fread"))
    origins = unlist(specification$sample_origins, use.names = FALSE)
    models = unlist(specification$models, use.names = FALSE)
    validate_prodigy_expression(counts, origins, models, network)
    stopifnot(identical(colnames(mutations), models), identical(rownames(counts), rownames(mutations)),
        all(mutations %in% c(0, 1)), specification$sample %in% models)
    pathways = readRDS(specification$pathways)
    selected_mutations = rownames(mutations)[mutations[, specification$sample] == 1]
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(as.integer(specification$patient_seed))
    state$rng$pipeline_entry = .Random.seed
    if (specification$mode == "factorial") {
        differences = read_explicit_evidence(specification$explicit_evidence, rownames(counts))
    }
    else {
        differences = execution$get_DEGs(expression_matrix = counts, samples = specification$sample,
            sample_origins = origins, beta = 2, gamma = 0.05)[[specification$sample]]
    }
    state$rng$after_DEG = .Random.seed
    saveRDS(list(differences = differences, rng = state$rng), file.path(output, "DEG_checkpoint.rds"),
        compress = FALSE)
    if (!is.null(differences) && is.double(differences) && !is.null(names(differences))) {
        write_json(list(schema = "cdf_explicit_evidence_v1", status = "AVAILABLE", genes = as.list(names(differences)),
            binary64_le = binary64(differences), size = length(differences), original_size = length(differences),
            non_null_empty = TRUE, restriction_preserves_original_order_and_bits = TRUE), file.path(output,
            "evidence.json"))
    }
    raw = execution$PRODIGY(mutated_genes = selected_mutations, expression_matrix = counts, network = network,
        sample = specification$sample, diff_genes = differences, alpha = specification$alpha, pathway_list = pathways,
        num_of_cores = 1, sample_origins = origins, write_results = FALSE, beta = 2, gamma = 0.05, delta = 0.05)
    state$rng$after_influence = .Random.seed
    stopifnot(identical(raw, state$axes$returned))
    if (is.null(raw)) {
        cause = if (length(state$axes$mutation_genes) == 0L)
            "NO_MUTATIONS_IN_NETWORK"
        else if (state$enrichment_completed && length(state$enriched) == 0L)
            "NO_ENRICHED_PATHWAYS"
        else "UNEXPLAINED_NULL"
        stopifnot(cause != "UNEXPLAINED_NULL")
        restored = list(status = "SUCCESSFUL_EMPTY", values = NULL)
    }
    else {
        cause = NULL
        restored = restore_prodigy_axes(raw, state$axes$mutation_genes, state$axes$evaluated_pathways)
        stopifnot(restored$status != "UNAVAILABLE_UNVERIFIED_AXES")
    }
    identity = specification[c("context", "sample", "network_id", "normal_pool", "master_seed", "patient_seed",
        "alpha")]
    saveRDS(list(specification = identity, differences = differences, raw_influence = raw, axes = state$axes,
        restored = restored, rng = state$rng), file.path(output, "influence_checkpoint.rds"), compress = FALSE)
    ranking = execution$analyze_PRODIGY_results(restored$values)
    state$rng$pipeline_exit = .Random.seed
    genes = if (length(ranking) == 0L || is.null(ranking[[1L]]))
        character()
    else ranking[[1L]]
    stopifnot(is.character(genes), !anyDuplicated(genes), all(genes %in% selected_mutations))
    scores = as.numeric(state$aggregation[genes])
    stopifnot(length(scores) == length(genes), all(is.finite(scores)), all(scores > 0))
    result = list(status = if (length(genes) > 0L) "SUCCESS" else "SUCCESS_EMPTY", context = specification$context,
        sample = specification$sample, sample_id = specification$sample, method = "PRODIGY", network_id = specification$network_id,
        cell = specification$cell, master_seed = specification$master_seed, patient_seed = specification$patient_seed,
        alpha = specification$alpha, parameter_id = paste0("alpha_", specification$alpha), candidate_count = length(genes),
        output_count = length(genes), genes = genes, scores = scores, scores_decimal17 = as.list(sprintf("%.17g",
            scores)), scores_little_endian_float64_hex = binary64(scores), reason = cause, empty_cause = cause,
        axis_status = restored$status, score_semantics = "Original positive aggregate influence before the native mixture cutoff",
        pre_method_mutation_count = length(selected_mutations), selected_mutations = length(selected_mutations),
        enriched_pathways = length(state$enriched))
    saveRDS(list(result = result, ranking = ranking, aggregate_influence = state$aggregation, rng = state$rng,
        mixtures = state$mixtures), file.path(output, "result.rds"), compress = FALSE)
    if (length(genes) == 0L && is.null(result$reason)) {
        result$reason = "Original postprocessor returned no candidates"
    }
    result$genes = as.list(result$genes)
    result$scores = as.list(result$scores)
    write_json(result, file.path(output, "result.json"))
    return(result)
}

run_prodigy = function(specification) {
    for (package in c("DESeq2", "igraph", "PCSF", "plyr", "MASS", "mixtools", "data.table")) {
        suppressPackageStartupMessages(library(package, character.only = TRUE))
    }
    outcomes = vector("list", length(specification$tasks))
    for (index in seq_along(specification$tasks)) {
        task = modifyList(specification, specification$tasks[[index]])
        task$tasks = NULL
        output = file.path(specification$output, sprintf("task_%06d", index))
        stopifnot(dir.create(output))
        state = new.env(parent = emptyenv())
        state$axes = NULL
        state$aggregation = NULL
        state$enriched = NULL
        state$enrichment_completed = FALSE
        state$mixtures = list()
        state$rng = list()
        result = tryCatch({
            if (identical(task$input_status, "UNAVAILABLE")) {
                value = list(status = "UNAVAILABLE", context = task$context, sample_id = task$sample,
                  method = "PRODIGY", cell = task$cell, network_id = task$network_id, master_seed = task$master_seed,
                  patient_seed = task$patient_seed, reason = task$input_reason)
                write_json(value, file.path(output, "result.json"))
                value
            }
            else {
                prodigy_patient(task, output, state)
            }
        }, error = function(error) {
            value = list(status = "FAILED", context = task$context, sample_id = task$sample, method = "PRODIGY",
                cell = task$cell, network_id = task$network_id, master_seed = task$master_seed, patient_seed = task$patient_seed,
                reason = conditionMessage(error))
            saveRDS(list(result = value, rng = state$rng, axes = state$axes, mixtures = state$mixtures),
                file.path(output, "failure.rds"), compress = FALSE)
            write_json(value, file.path(output, "result.json"))
            return(value)
        })
        outcomes[[index]] = list(task = index, cell = task$cell, sample_id = task$sample, master_seed = task$master_seed,
            status = result$status, result = file.path(basename(output), "result.json"))
    }
    write_json(outcomes, file.path(specification$output, "task_results.json"))
    return(!any(vapply(outcomes, function(value) value$status == "FAILED", logical(1))))
}
