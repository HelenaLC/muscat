#' @rdname aggregateData
#' @title Aggregation of single-cell to pseudo-bulk data
#' 
#' @description ...
#' 
#' @param x a \code{[SingleCellExperiment]{SingleCellExperiment}}.
#' @param data a character string. 
#'   Specifies the assay slot to use as input data.
#' @param fun a character string.
#'   Specifies the function to use as summary statistic.
#' @param scale logical.
#' 
#' @return a list of sample-wise pseudo-bulk data for each cluster.
#' 
#' @author Helena L. Crowel \email{helena.crowell@uzh.ch} and Mark D. Robinson.
#' 
#' @importFrom Matrix colSums rowMeans rowSums
#' @importFrom matrixStats rowMedians
#' @importFrom methods is
#' @importFrom SummarizedExperiment assayNames assays colData
#' 
#' @export

aggregateData <- function(x, data, 
    fun = c("sum", "mean", "median"), 
    scale = FALSE) {
    
    # validity checks for input arguments
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot(c("cluster_id", "sample_id") %in% colnames(colData(x)))
    stopifnot(is.character(data), length(data) == 1, data %in% assayNames(x))
    stopifnot(is.logical(scale), length(scale) == 1)

    # get aggregation function
    fun <- match.arg(fun)
    fun <- switch(fun,
        sum = "rowSums",
        mean = "rowMeans",
        median = "rowMedians")
    
    # split cells by cluster-sample
    dt <- data.table(
        cell = colnames(x),
        cluster_id = colData(x)$cluster_id,
        sample_id = colData(x)$sample_id)
    idx <- split(dt, 
        by = c("cluster_id", "sample_id"), 
        sorted = TRUE, flatten = FALSE)
    idx <- lapply(idx, lapply, "[[", "cell")
    
    # compute pseudo-bulks
    pb <- lapply(idx, vapply, function(i)
        get(fun)(assays(x)[[data]][, i, drop = FALSE]),
        numeric(nrow(x)))
    
    # scale
    if (scale) {
        if (data == "counts" & fun == "rowSums") {
            pb_counts <- pb
        } else {
            pb_counts <- lapply(idx, vapply, function(i)
                rowSums(assays(x)[[data]][, i, drop = FALSE]),
                numeric(nrow(x)))
        }
        n_samples <- nlevels(dt$sample_id)
        lib_sizes <- vapply(pb, colSums, numeric(n_samples))
        for (k in names(pb))
            pb[[k]] <- pb[[k]] * lib_sizes[, k] / 1e6
    }
    return(pb)
}
