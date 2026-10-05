normalise_counts = function(counts) {
    stopifnot(all(is.finite(counts)), all(counts >= 0), all(colSums(counts) > 0))
    value = edgeR::DGEList(counts = counts, lib.size = colSums(counts), remove.zeros = FALSE)
    value = edgeR::calcNormFactors(value, method = "TMM", refColumn = NULL, logratioTrim = 0.3,
        sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e10)
    result = edgeR::cpm(value, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2)
    stopifnot(all(is.finite(result)), identical(dimnames(result), dimnames(counts)))
    return(list(values = result, samples = value$samples))
}

read_count_matrix = function(path) {
    handle = if (grepl("[.]gz$", path)) gzfile(path, "rt") else file(path, "rt")
    on.exit(close(handle))
    result = as.matrix(read.csv(handle, row.names = 1L, check.names = FALSE,
        stringsAsFactors = FALSE))
    storage.mode(result) = "double"
    stopifnot(!anyDuplicated(rownames(result)), !anyDuplicated(colnames(result)),
        all(is.finite(result)))
    return(result)
}

prepare_oesophageal = function(specification) {
    stopifnot(as.character(getRversion()) == "4.5.2",
        as.character(packageVersion("Matrix")) == "1.7.4",
        as.character(packageVersion("edgeR")) == "4.8.2")
    work = specification$output
    inputs = jsonlite::read_json(specification$inputs, simplifyVector = TRUE)
    models = inputs$models
    normals = inputs$normal_samples
    stopifnot(nrow(models) == 59L, length(unique(models$donor_id)) == 58L,
        length(normals) == 403L, length(unique(normals)) == 403L)
    tumour = read_count_matrix(specification$tumour_counts)
    normal = read_count_matrix(specification$normal_counts)
    mutation = read_count_matrix(specification$mutation_binary)
    stopifnot(identical(colnames(tumour), models$sample_id),
        identical(colnames(mutation), models$sample_id), identical(colnames(normal), normals),
        identical(rownames(tumour), rownames(normal)), nrow(tumour) == 39102L)
    stopifnot(all(tumour >= 0), all(normal >= 0), all(tumour == floor(tumour)),
        all(normal == floor(normal)), all(mutation %in% c(0, 1)),
        length(intersect(colnames(tumour), colnames(normal))) == 0L)
    colnames(normal) = paste0("NORMAL_", normals)
    stopifnot(length(intersect(colnames(tumour), colnames(normal))) == 0L)
    normalised = normalise_counts(cbind(normal, tumour))
    write.csv(data.frame(sample = rownames(normalised$samples), normalised$samples,
        check.names = FALSE), file.path(work, "normalisation_factors.csv"), row.names = FALSE)
    saveRDS(normalised, file.path(work, "normalised_counts.rds"), compress = FALSE, version = 3)
    original = readLines(specification$common_genes, warn = FALSE)
    support = intersect(rownames(tumour), rownames(mutation))
    genes = original[original %in% support]
    stopifnot(length(original) == 7399L, length(genes) == 7388L,
        anyDuplicated(genes) == 0L, identical(genes, inputs$shared_original_query))
    writeLines(genes, file.path(work, "genes.txt"))
    writeLines(models$sample_id, file.path(work, "samples.txt"))
    write.csv(models, file.path(work, "models.csv"), row.names = FALSE)
    writeLines(normals, file.path(work, "normal_ids.txt"))
    write.csv(data.frame(gene = genes, mutation[genes, , drop = FALSE], check.names = FALSE),
        file.path(work, "mutation.csv"), row.names = FALSE)
    original_method = new.env(parent = .GlobalEnv)
    base::sys.source(specification$normalise_source, envir = original_method)
    evidence = list()
    graphs = list()
    diagnostics = list()
    for (version in 0:1) {
        native = readLines(specification[[paste0("native_genes_H", version)]], warn = FALSE)
        native = native[native %in% support]
        values = normalised$values[native, , drop = FALSE]
        expression = original_method$DawnNormalize(
            values[, models$sample_id, drop = FALSE],
            values[, paste0("NORMAL_", normals), drop = FALSE])
        stopifnot(all(is.finite(expression)), all(expression >= 0), all(genes %in% native))
        selected = expression[genes, , drop = FALSE]
        stopifnot(all(colSums(selected) > 0))
        evidence[[paste0("H", version)]] = selected
        writeLines(native, file.path(work, paste0("native_genes_H", version, ".txt")))
        edges = read_edges(specification[[paste0("edges_G", version)]])
        adjacency = aligned_dawn_graph(edges, genes)
        degree = Matrix::colSums(adjacency)
        transition = adjacency %*% Matrix::Diagonal(x = 1 / (degree + 1e-16))
        dimnames(transition) = dimnames(adjacency)
        damping = degree / (degree + 3)
        graphs[[paste0("G", version)]] = list(transition = transition, damping = damping,
            degree = degree)
        diagnostics[[version + 1L]] = list(native_genes = length(native),
            edges = length(adjacency@x) / 2, isolates = sum(degree == 0))
    }
    stopifnot(identical(vapply(diagnostics, function(value) value$native_genes, integer(1)),
        c(8994L, 9538L)))
    prepared = list(genes = genes, models = models, normals = normals,
        mutation = mutation[genes, , drop = FALSE], graphs = graphs, evidence = evidence)
    saveRDS(prepared, file.path(work, "prepared.rds"), compress = FALSE, version = 3)
    jsonlite::write_json(list(models = 59L, donors = 58L, normals = 403L,
        genes = length(genes), graph_states = diagnostics, status = "PREPARED"),
        file.path(work, "preparation.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
    return(invisible(prepared))
}
