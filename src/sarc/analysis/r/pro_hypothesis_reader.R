arguments = commandArgs(trailingOnly = TRUE)

stopifnot(length(arguments) == 1L, as.character(getRversion()) == "4.5.2")

invisible(Sys.setlocale("LC_COLLATE", "C"))

hex = function(value) paste(sprintf("%02x", as.integer(value)), collapse = "")

strings = function(value) {
    if (is.null(value))
        return("null")
    paste0("[", paste(vapply(value, function(x) if (is.na(x))
        "null"
    else paste0("\"", hex(charToRaw(enc2utf8(x))), "\""), character(1)), collapse = ","), "]")
}

node = function(value) {
    if (is.null(value))
        return("{\"t\":\"null\"}")
    kind = typeof(value)
    stopifnot(kind %in% c("double", "integer", "logical", "character", "list"))
    prefix = paste0("{\"t\":\"", kind, "\",\"names\":", strings(names(value)))
    dimensions = dim(value)
    prefix = paste0(prefix, ",\"dim\":", if (is.null(dimensions))
        "null"
    else paste0("[", paste(dimensions, collapse = ","), "]"))
    labels = dimnames(value)
    prefix = paste0(prefix, ",\"dimnames\":", if (is.null(labels))
        "null"
    else paste0("[", paste(vapply(labels, strings, character(1)), collapse = ","), "]"))
    if (kind == "list")
        return(paste0(prefix, ",\"values\":[", paste(lapply(value, node), collapse = ","), "]}"))
    if (kind == "character")
        return(paste0(prefix, ",\"values\":", strings(as.vector(value)), "}"))
    if (kind == "double") {
        states = rep("F", length(value))
        states[is.na(value)] = "A"
        states[is.nan(value)] = "N"
        states[is.infinite(value) & value > 0] = "P"
        states[is.infinite(value) & value < 0] = "M"
        return(paste0(prefix, ",\"bits\":\"", hex(writeBin(as.vector(value), raw(), size = 8L, endian = "little")),
            "\",\"states\":\"", paste(states, collapse = ""), "\"}"))
    }
    paste0(prefix, ",\"bits\":\"", hex(writeBin(as.integer(value), raw(), size = 4L, endian = "little")),
        "\"}")
}

select_value = function(value, mode) {
    if (mode == "checkpoint") {
        stopifnot(is.list(value), all(c("specification", "differences", "raw_influence", "axes", "restored") %in%
            names(value)))
        identity = value$specification[c("context", "sample", "network_id", "normal_pool", "master_seed",
            "patient_seed", "alpha")]
        return(list(identity = identity, differences = value$differences, raw_influence = value$raw_influence,
            axes = value$axes, restored = value$restored))
    }
    if (mode == "result") {
        stopifnot(is.list(value), "result" %in% names(value))
        return(list(result = value$result, ranking = value$ranking, aggregate = value$aggregate_influence,
            provenance_mode = value$provenance_mode, provenance = value$provenance, differences = value$differences))
    }
    if (mode == "pathways") {
        stopifnot(is.list(value), !is.null(names(value)), !anyDuplicated(names(value)), !anyNA(names(value)))
        return(lapply(value, function(pathway) {
            if (!(is.matrix(pathway) || identical(class(pathway), "data.frame"))) return(list(status = "UNAVAILABLE_UNSUPPORTED_PATHWAY_CLASS",
                class = class(pathway)))
            stopifnot(ncol(pathway) >= 2L)
            members = unique(c(pathway[, 1], pathway[, 2]))
            if (!is.character(members) || anyNA(members)) return(list(status = "UNAVAILABLE_NONCHARACTER_PATHWAY_MEMBERS",
                class = class(pathway)))
            list(status = "EXACT_SOURCE_MEMBERSHIP", members = members)
        }))
    }
    stop("Unsupported scientific RDS projection")
}

request = readLines(arguments[[1]], encoding = "UTF-8", warn = FALSE)

for (line in request) {
    fields = strsplit(line, "\t", fixed = TRUE)[[1]]
    stopifnot(length(fields) == 3L, !file.exists(fields[[3]]))
    value = readRDS(fields[[2]])
    selected = select_value(value, fields[[1]])
    connection = file(fields[[3]], open = "wb")
    tryCatch(writeBin(charToRaw(paste0(node(selected), "\n")), connection), finally = close(connection))
}
