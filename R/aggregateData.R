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
#' @import SingleCellExperiment
#' @importFrom Matrix rowMeans rowSums
#' @importFrom matrixStats rowMedians
#' @importFrom methods is
#' @export

aggregateData <- function(x, data, 
    fun = c("sum", "mean", "median"), 
    scale = FALSE) {
    
    # validity checks for input arguments
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot(all(c("cluster_id", "sample_id") %in% colnames(colData(x))))
    stopifnot(is.character(data), length(data) == 1, data %in% assayNames(x))
    stopifnot(is.logical(scale), length(scale) == 1)

    # get aggregation function
    fun <- match.arg(fun)
    fun <- switch(fun,
        sum = "rowSums",
        mean = "rowMeans",
        median = "rowMedians")
    fun <- get(fun)
    
    # split cells by cluster-sample
    cluster_ids <- factor(colData(x)$cluster_id)
    sample_ids <- factor(colData(x)$sample_id)
    idx <- split(seq_len(ncol(x)), list(cluster_ids, sample_ids))
    
    # compute pseudo-bulks
    pb <- sapply(idx, function(i) 
        fun(assays(x)[[data]][, i, drop = FALSE]))
    
    # scale
    if (scale) {
        if (data == "counts" & fun == "sum") {
            pb_counts <- pb
        } else {
            pb_counts <- sapply(idx, function(i) 
                rowSums(assays(x)[[data]][, i, drop = FALSE]))
        }
        lib_sizes <- colSums(pb_counts)
        pb <- pb * lib_sizes / 1e6
    }
    
    # unwrap & split by cluster
    df <- data.frame(
        index = seq_len(ncol(pb)),
        cluster_id = rep(levels(cluster_ids), nlevels(sample_ids)),
        sample_id = rep(levels(sample_ids), each = nlevels(cluster_ids)))
    dfs <- split(df, df$cluster_id)
    
    lapply(dfs, function(df) {
        x <- pb[, df$index, drop = FALSE]
        colnames(x) <- df$sample_id
        return(x)
    })
}
