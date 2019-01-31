#' @rdname aggregateData
#' @title Aggregation of single-cell to pseudo-bulk data
#' 
#' @description ...
#' 
#' @param x a \code{[SingleCellExperiment]{SingleCellExperiment}}.
#' @param assay a character string. 
#'   Specifies the assay slot to use as input data.
#' @param fun a character string.
#'   Specifies the function to use as summary statistic.
#' @param scale logical.
#' 
#' @return a list of sample-wise pseudo-bulk data for each cluster.
#' 
#' @examples 
#' data(kang)
#' 
#' pb <- aggregateData(kang, assay = "counts", fun = "sum")
#' names(pb)
#' head(pb[[1]])
#' 
#' counts <- assay(kang)
#' assays(kang)$cpm <- edgeR::cpm(counts)
#' pb <- aggregateData(kang, assay = "cpm", fun = "sum", scale = TRUE)
#' head(pb[[1]])
#' 
#' @author Helena L. Crowel \email{helena.crowell@uzh.ch} and Mark D. Robinson.
#' 
#' @importFrom Matrix colSums rowMeans rowSums
#' @importFrom matrixStats rowMedians
#' @importFrom methods is
#' @importFrom SummarizedExperiment assayNames assays colData
#' 
#' @export

aggregateData <- function(x, assay, 
    fun = c("sum", "mean", "median"), 
    scale = FALSE) {
    
    # validity checks for input arguments
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot(c("cluster_id", "sample_id") %in% colnames(colData(x)))
    stopifnot(is.character(assay), length(assay) == 1, assay %in% assayNames(x))
    stopifnot(is.logical(scale), length(scale) == 1)

    # get aggregation function
    fun <- match.arg(fun)
    fun <- switch(fun,
        sum = "rowSums",
        mean = "rowMeans",
        median = "rowMedians")
    
    # split cells by cluster-sample
    cells_by_cluster_sample <- .split_cells(x)
    
    # compute pseudo-bulks
    pb <- lapply(cells_by_cluster_sample, vapply, function(i)
        get(fun)(assays(x)[[assay]][, i, drop = FALSE]),
        numeric(nrow(x)))
    
    # scale
    if (scale) {
        if (assay == "counts" & fun == "rowSums") {
            pb_counts <- pb
        } else {
            pb_counts <- lapply(cells_by_cluster_sample, vapply, function(i)
                rowSums(assays(x)[[assay]][, i, drop = FALSE]),
                numeric(nrow(x)))
        }
        n_samples <- nlevels(colData(x)$sample_id)
        lib_sizes <- vapply(pb_counts, colSums, numeric(n_samples))
        cluster_ids <- levels(colData(x)$cluster_id)
        names(cluster_ids) <- cluster_ids
        pb <- lapply(cluster_ids, function(k) 
          pb[[k]] / lib_sizes[, k] * 1e6)
    }
    return(pb)
}
