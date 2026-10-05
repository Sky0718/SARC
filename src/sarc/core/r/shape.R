restore_prodigy_axes = function(raw, mutation_genes, evaluated_pathways) {
    unavailable = function(reason) list(status = "UNAVAILABLE_UNVERIFIED_AXES", reason = reason, values = NULL)
    valid_axis = function(axis) is.character(axis) && length(axis) > 0L && !anyNA(axis) && all(nchar(axis) >
        0L) && anyDuplicated(axis) == 0L
    if (!valid_axis(mutation_genes) || !valid_axis(evaluated_pathways)) {
        return(unavailable("Complete unique authoritative mutation and evaluated-pathway axes are required"))
    }
    if (is.null(raw) || !is.numeric(raw)) {
        return(unavailable("No numeric original influence result"))
    }
    expected_dimensions = c(length(evaluated_pathways), length(mutation_genes))
    if (is.matrix(raw)) {
        if (!identical(dim(raw), expected_dimensions) || !identical(rownames(raw), evaluated_pathways) ||
            !identical(colnames(raw), mutation_genes)) {
            return(unavailable("Original matrix does not match the authoritative axes"))
        }
        return(list(status = "AXES_VERIFIED", reason = NULL, values = raw))
    }
    if (!is.null(dim(raw)) || length(raw) != prod(expected_dimensions)) {
        return(unavailable("Vector length cannot establish the declared matrix"))
    }
    if (length(mutation_genes) == 1L && length(evaluated_pathways) > 1L && identical(names(raw), evaluated_pathways)) {
        restored = matrix(unname(raw), ncol = 1L, dimnames = list(evaluated_pathways, mutation_genes))
    }
    else if (length(evaluated_pathways) == 1L && length(mutation_genes) > 1L && identical(names(raw),
        mutation_genes)) {
        return(list(status = "VALID_ORIGINAL_VECTOR_BRANCH_PRESERVED", reason = NULL, values = raw))
    }
    else if (all(expected_dimensions == 1L) && (is.null(names(raw)) || identical(names(raw), evaluated_pathways) ||
        identical(names(raw), mutation_genes))) {
        return(list(status = "SCALAR_GENE_IDENTITY_VERIFIED", reason = NULL, values = setNames(as.vector(raw),
            mutation_genes)))
    }
    else {
        return(unavailable("Dropped vector lacks an exact authoritative singleton-axis interpretation"))
    }
    stopifnot(identical(as.vector(restored), as.vector(raw)))
    return(list(status = "AXES_RESTORED_EXACT_NUMERICAL_VALUES", reason = NULL, values = restored))
}
