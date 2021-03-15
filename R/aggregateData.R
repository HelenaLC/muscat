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
#' @param BPPARAM a \code{\link[BiocParallel]{BiocParallelParam}}
#'   object specifying how aggregation should be parallelized.
#' @param verbose logical. Should information on progress be reported?
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
#'   are accessible in \code{int_colData()$n_cells}.
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
#' @importFrom Matrix colSums
#' @importFrom purrr map
#' @importFrom S4Vectors DataFrame metadata
#' @importFrom SingleCellExperiment SingleCellExperiment int_colData<-
#' @importFrom SummarizedExperiment colData colData<-
#' @export

aggregateData <- function(x,
    assay = NULL, by = c("cluster_id", "sample_id"),
    fun = c("sum", "mean", "median"), scale = FALSE, 
    verbose = TRUE, BPPARAM = SerialParam(progressbar = verbose)) {
    
    # check validity of input arguments
    fun <- match.arg(fun)
    if (is.null(assay)) 
        assay <- assayNames(x)[1] 
    .check_arg_assay(x, assay)
    .check_args_aggData(as.list(environment()))
    stopifnot(is(BPPARAM, "BiocParallelParam"))
    
    # assure 'by' colData columns are factors
    # so that missing combinations aren't dropped
    for (i in by) 
        if (!is.factor(x[[i]])) 
            x[[i]] <- factor(x[[i]])
    
    # compute pseudo-bulks
    pb <- .pb(x, by, assay, fun, BPPARAM)
    if (scale & length(by) == 2) {
        # compute library sizes
        cs <- if (assay == "counts" && fun == "sum")
            pb else .pb(x, by, "counts", "sum", BPPARAM)
        ls <- lapply(cs, colSums)
        # scale pseudobulks by CPM
        pb <- lapply(seq_along(pb), function(i) pb[[i]] / 1e6 * ls[[i]])
        names(pb) <- names(ls)
    }
    
    # construct SCE
    md <- metadata(x)
    md$agg_pars <- list(assay = assay, by = by, fun = fun, scale = scale)
    pb <- SingleCellExperiment(pb, metadata = md)
    
    # tabulate number of cells
    cd <- data.frame(colData(x)[, by])
    for (i in names(cd))
        if (is.factor(cd[[i]]))
            cd[[i]] <- droplevels(cd[[i]])
    ns <- table(cd)
    if (length(by) == 2) {
        ns <- asplit(ns, 2)
        ns <- map(ns, ~c(unclass(.)))
    } else ns <- c(unclass(ns))
    int_colData(pb)$n_cells <- ns
    
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
