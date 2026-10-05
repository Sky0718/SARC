arguments = commandArgs(trailingOnly = TRUE)

stopifnot(length(arguments) == 1L)

specification = jsonlite::fromJSON(arguments[[1L]], simplifyVector = FALSE)

if (!is.null(specification$library_paths)) {
    .libPaths(unlist(specification$library_paths, use.names = FALSE))
}

stopifnot(getRversion() == "4.5.2")

Sys.setlocale("LC_COLLATE", "C")

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

suppressPackageStartupMessages(library(data.table))

data.table::setDTthreads(1L)

entry = grep("^--file=", commandArgs(), value = TRUE)

directory = dirname(normalizePath(sub("^--file=", "", entry[[1L]]), mustWork = TRUE))

source(file.path(directory, "common.R"))

source(file.path(directory, "shape.R"))

source(file.path(directory, "..", "..", "precision", "certificate.R"))

source(file.path(directory, "..", "..", "precision", "solver.R"))

source(file.path(directory, "dawn.R"))

source(file.path(directory, "prodigy.R"))

for (package in names(specification$versions)) {
    stopifnot(as.character(packageVersion(package)) == specification$versions[[package]])
}

status = tryCatch({
    if (specification$method == "PRODIGY") {
        run_prodigy(specification)
    }
    else if (specification$method == "DawnRank") {
        run_dawn(specification)
    }
    else {
        stop("Unsupported R method")
    }
}, error = function(error) {
    write_json(list(status = "ERROR", error = conditionMessage(error)), file.path(specification$output,
        "failure.json"))
    return(FALSE)
})

writeLines(capture.output(sessionInfo()), file.path(specification$output, "R_session.txt"))

quit(status = if (isTRUE(status)) 0L else 1L)
