#' @rdname aggregateData
#' @title Aggregation of single-cell to pseudobulk data
#' 
#' @description ...
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param assay character string specifying the assay slot to use as 
#'   input data. Defaults to the 1st available (\code{assayNames(x)[1]}).
#' @param by character vector specifying which 
#'   \code{colData(x)} columns to summarize by (at most 2!).
#' @param fun a character string.
#'   Specifies the function to use as summary statistic.
#' @param scale logical. Should pseudo-bulks be scaled
#'   with the effective library size & multiplied by 1M?
#' 
#' @return a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' \itemize{
#' \item{If \code{length(by) == 2}, each sheet (\code{assay}) contains 
#'   pseudobulks for each of \code{by[1]}, e.g., for each cluster when 
#'   \code{by = "cluster_id"}. Rows correspond to genes, columns to 
#'   \code{by[2]}, e.g., samples when \code{by = "sample_id"}}.
#' \item{If \code{length(by) == 1}, the returned SCE will contain only 
#'   a single \code{assay} with rows = genes and colums = \code{by}.}}
#'   
#'   Aggregation parameters (\code{assay, by, fun, scaled}) are stored in 
#'   \code{metadata()$agg_pars}, and the number of cells that were aggregated 
#'   are accessible in \code{metadata()$n_cells}.
#' 
#' @examples 
#' data(sce)
#' library(SingleCellExperiment)
#' 
#' # pseudobulk counts by cluster-sample
#' pb <- aggregateData(sce)
#' 
#' assayNames(sce)  # one sheet per cluster
#' head(assay(sce)) # n_genes x n_samples
#' 
#' # scaled CPM
#' assays(sce)$cpm <- edgeR::cpm(assay(sce))
#' pb <- aggregateData(sce, assay = "cpm", scale = TRUE)
#' head(assay(pb)) 
#' 
#' # aggregate by cluster only
#' pb <- aggregateData(sce, by = "cluster_id")
#' length(assays(pb)) # single assay
#' head(assay(pb))    # n_genes x n_clusters
#' 
#' @author Helena L Crowell & Mark D Robinson
#' 
#' @references 
#' Crowell, HL, Soneson, C, Germain, P-L, Calini, D, 
#' Collin, L, Raposo, C, Malhotra, D & Robinson, MD: 
#' On the discovery of population-specific state transitions from 
#' multi-sample multi-condition single-cell RNA sequencing data. 
#' \emph{bioRxiv} \strong{713412} (2018). 
#' doi: \url{https://doi.org/10.1101/713412}
#' 
#' @importFrom Matrix colSums rowSums rowMeans
#' @importFrom matrixStats rowMedians
#' @importFrom S4Vectors DataFrame metadata
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData colData<-
#' @export

aggregateData <- function(x,
    assay = NULL, by = c("cluster_id", "sample_id"),
    fun = c("sum", "mean", "median"), scale = FALSE) {
    
    # check validity of input arguments
    fun <- match.arg(fun)
    if (is.null(assay)) 
        assay <- assayNames(x)[1] 
    .check_arg_assay(x, assay)
    .check_args_aggData(as.list(environment()))
    
    # store aggregation parameters &
    # nb. of cells that went into aggregation
    md <- metadata(x)
    md$agg_pars <- list(assay = assay, by = by, fun = fun, scale = scale)
    md$n_cells <- table(as.data.frame(colData(x)[, by]))
    
    # get aggregation function
    fun <- switch(fun,
        sum = "rowSums",
        mean = "rowMeans",
        median = "rowMedians")
    
    # split cells & compute pseudo-bulks
    cs <- .split_cells(x, by)
    pb <- .pb(x, cs, assay, fun)
    if (scale & length(by) == 2) {
        ls <- lapply(.pb(x, cs, "counts", "rowSums"), colSums)
        pb <- lapply(seq_along(pb), function(i) pb[[i]] / 1e6 * ls[[i]])
        names(pb) <- names(ls)
    }
    
    # construct SCE
    pb <- SingleCellExperiment(pb, metadata = md)
    
    # propagate 'colData' columns that are unique across 2nd 'by'
    if (length(by) == 2) {
        cd <- colData(x)
        ids <- colnames(pb)
        counts <- vapply(ids, function(u) {
            m <- as.logical(match(cd[, by[2]], u, nomatch = 0))
            vapply(cd[m, ], function(u) length(unique(u)), numeric(1))
        }, numeric(ncol(colData(x))))
        cd_keep <- apply(counts, 1, function(u) all(u == 1))
        cd_keep <- setdiff(names(which(cd_keep)), by)
        if (length(cd_keep) != 0) {
            m <- match(ids, cd[, by[2]], nomatch = 0)
            cd <- cd[m, cd_keep, drop = FALSE]
            rownames(cd) <- ids
            colData(pb) <- cd
        }
    }
    return(pb)
}
