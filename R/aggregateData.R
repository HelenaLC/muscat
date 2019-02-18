#' @rdname aggregateData
#' @title Aggregation of single-cell to pseudo-bulk data
#' 
#' @description ...
#' 
#' @param x a \code{[SingleCellExperiment]{SingleCellExperiment}}.
#' @param assay character string specifying the assay slot to use as 
#'   input data. Defaults to the 1st available (\code{assayNames(x)[1]}).
#' @param by character vector specifying which 
#'   \code{colData(x)} columns to summarize by.
#' @param fun a character string.
#'   Specifies the function to use as summary statistic.
#' @param scale logical. Should pseudo-bulks be scaled
#'   with the effective library size & multiplied by 1M?
#' 
#' @return a \code{[SingleCellExperiment]{SingleCellExperiment}}.
#' 
#' @examples 
#' library(SummarizedExperiment)
#' data(kang)
#' 
#' # pseudo-bulk counts by cluster-sample
#' pb <- aggregateData(kang)
#' assayNames(pb)  # one sheet per cluster
#' head(assay(pb)) # n_genes x n_samples
#' 
#' # scaled CPM
#' assays(kang)$cpm <- edgeR::cpm(assay(kang))
#' pb <- aggregateData(kang, assay = "cpm", scale = TRUE)
#' head(assay(pb)) 
#' 
#' # aggregate by cluster
#' pb <- aggregateData(kang, by = "cluster_id")
#' length(assays(pb)) # single assay
#' head(assay(pb))    # n_genes x n_clusters
#' 
#' @author Helena L. Crowell \email{helena.crowell@uzh.ch} and Mark D. Robinson.
#' 
#' @importFrom Matrix colSums rowMeans rowSums
#' @importFrom matrixStats rowMedians
#' @importFrom methods is
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData
#' 
#' @export

aggregateData <- function(x, assay,
    by = c("cluster_id", "sample_id"),
    fun = c("sum", "mean", "median"), 
    scale = FALSE) {
    
    if (missing("assay"))
        assay <- assayNames(x)[1]

    # validity checks for input arguments
    .check_sce(x, req_group = FALSE)
    .check_arg_assay(x, assay)
    stopifnot(is.character(by), by %in% colnames(colData(x)))
    stopifnot(is.logical(scale), length(scale) == 1)
    
    # get aggregation function
    fun <- match.arg(fun)
    fun <- switch(fun,
        sum = "rowSums",
        mean = "rowMeans",
        median = "rowMedians")
    
    # split cells & compute pseudo-bulks
    cs <- .split_cells(x, by)
    pb <- .pb(x, cs, assay, fun)
    
    # scale
    if (scale) {
        if (assay == "counts" & fun == "rowSums") {
            pb_counts <- pb
        } else {
            pb_counts <- .pb(x, cs, assay, "rowSums")
        }
        pb <- map_depth(pb_counts, -2, function(u)
            u / colSums(u)[col(u)] * 1e6)
    }
    pb <- map_depth(pb, -2, as.matrix)
    
    # return SCE
    md <- metadata(x)
    md$agg_pars <- list(assay = assay, fun = fun)
    
    SingleCellExperiment(
        assays = pb,
        metadata = md)
}
