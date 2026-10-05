write_json = function(value, path) {
    jsonlite::write_json(value, path, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null", na = "null")
}

read_numeric_matrix = function(path) {
    value = data.table::fread(path, data.table = FALSE, check.names = FALSE)
    genes = as.character(value[[1L]])
    value = as.matrix(value[, -1L, drop = FALSE])
    storage.mode(value) = "double"
    rownames(value) = genes
    stopifnot(!anyNA(genes), !anyDuplicated(genes), all(nchar(genes) > 0L), !anyDuplicated(colnames(value)))
    return(value)
}

read_native_dawn_matrix = function(path) {
    return(as.matrix(read.csv(path, row.names = 1L, check.names = FALSE, stringsAsFactors = FALSE)))
}

read_edges = function(path, parser = "base") {
    edges = if (parser == "fread")
        data.table::fread(path, data.table = FALSE, check.names = FALSE)
    else read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    stopifnot(identical(colnames(edges), c("gene1", "gene2", "combined_score")))
    stopifnot(all(is.finite(edges$combined_score)), all(edges$combined_score > 700), all(edges$combined_score <=
        1000), all(edges$combined_score == round(edges$combined_score)))
    stopifnot(all(edges$gene1 != edges$gene2), !anyDuplicated(paste(pmin(edges$gene1, edges$gene2), pmax(edges$gene1,
        edges$gene2), sep = "\t")))
    return(edges)
}

aligned_dawn_graph = function(edges, genes) {
    stopifnot(length(genes) > 0L, !anyDuplicated(genes))
    selected = edges[(edges$gene1 %in% genes) & (edges$gene2 %in% genes), , drop = FALSE]
    first = match(selected$gene1, genes)
    second = match(selected$gene2, genes)
    return(Matrix::sparseMatrix(i = c(first, second), j = c(second, first), x = rep(1, 2 * nrow(selected)),
        dims = c(length(genes), length(genes)), dimnames = list(genes, genes)))
}

normalise_dawn_input = function(tumour, normal, original_tumour_names, original_normal_names, original_normalise) {
    stopifnot(is.matrix(tumour), is.matrix(normal), identical(rownames(tumour), rownames(normal)))
    stopifnot(ncol(tumour) == length(original_tumour_names), ncol(normal) == length(original_normal_names))
    stopifnot(!anyDuplicated(original_tumour_names), !anyDuplicated(original_normal_names), all(is.finite(tumour)),
        all(is.finite(normal)))
    canonical = colnames(tumour)
    stopifnot(!is.null(canonical), length(canonical) == ncol(tumour), !anyDuplicated(canonical), !anyDuplicated(rownames(tumour)))
    if (length(intersect(original_tumour_names, original_normal_names)) == 1L) {
        stop("The original single-matched-column drop requires a separately declared compatibility decision")
    }
    colnames(tumour) = original_tumour_names
    colnames(normal) = original_normal_names
    normalised = original_normalise(tumorMat = tumour, normalMat = normal)
    stopifnot(identical(dim(normalised), dim(tumour)), identical(rownames(normalised), rownames(tumour)))
    stopifnot(!anyDuplicated(colnames(normalised)), setequal(colnames(normalised), original_tumour_names))
    normalised = normalised[, original_tumour_names, drop = FALSE]
    colnames(normalised) = canonical
    evaluable = vapply(seq_len(ncol(normalised)), function(column) all(is.finite(normalised[, column])),
        logical(1))
    names(evaluable) = canonical
    status = ifelse(evaluable, "INPUT_READY", "UNAVAILABLE_NONFINITE_ORIGINAL_NORMALISATION")
    return(list(values = normalised, sample_status = status, evaluable = evaluable))
}

original_prodigy_network = function(edges) {
    network = as.matrix(data.frame(src = edges$gene1, dest = edges$gene2, score = edges$combined_score / 1000,
        stringsAsFactors = FALSE))
    stopifnot(identical(colnames(network), c("src", "dest", "score")))
    return(network)
}

validate_prodigy_expression = function(counts, origins, tumours, network) {
    stopifnot(is.matrix(counts), length(origins) == ncol(counts), all(origins %in% c("normal", "tumor")))
    stopifnot(all(is.finite(counts)), all(counts >= 0), !anyDuplicated(rownames(counts)), !anyDuplicated(colnames(counts)))
    stopifnot(sum(origins == "normal") > 0L, !anyDuplicated(tumours), setequal(colnames(counts)[origins ==
        "tumor"], tumours))
    stopifnot(all(rownames(counts) %in% unique(c(network[, "src"], network[, "dest"]))))
    return(invisible(NULL))
}

empty_frame = function(fields) {
    return(as.data.frame(setNames(replicate(length(fields), character(), simplify = FALSE), fields),
        stringsAsFactors = FALSE))
}

append_csv = function(value, path) {
    write.table(value, path, sep = ",", row.names = FALSE, col.names = FALSE, quote = TRUE, append = TRUE)
}
