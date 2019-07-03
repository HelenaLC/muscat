#' @rdname aggregateData
#' @title Aggregation of single-cell to pseudo-bulk data
#' 
#' @description ...
#' 
#' @param x a \code{[SingleCellExperiment]{SingleCellExperiment}}.
#' @param assay character string specifying the assay slot to use as 
#'   input data. Defaults to the 1st available (\code{assayNames(x)[1]}).
#' @param by character vector specifying which 
#'   \code{colData(x)} columns to summarize by (at most 2!).
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
#' @importFrom dplyr last
#' @importFrom Matrix colSums
#' @importFrom purrr map_depth
#' @importFrom S4Vectors metadata
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData colData<-
#' 
#' @export

aggregateData <- function(x, 
    assay = NULL,
    by = c("cluster_id", "sample_id"),
    fun = c("sum", "mean", "median"), 
    scale = FALSE) {
   
    # validity checks for input arguments
    if (is.null(assay)) {
        assay <- assayNames(x)[1] 
    } else {
        .check_arg_assay(x, assay)
    }
    stopifnot(is.character(by), by %in% colnames(colData(x)), length(by) <= 2)
    stopifnot(is.logical(scale), length(scale) == 1)

    # store aggregation parameters &
    # nb. of cells that went into aggregation
    md <- metadata(x)
    md$agg_pars <- list(assay = assay, fun = fun, scale = scale)
    md$n_cells <- table(as.data.frame(colData(x)[, by]))
    
    # get aggregation function
    fun <- switch(match.arg(fun),
        sum = "rowSums",
        mean = "rowMeans",
        median = "rowMedians")
    
    # split cells & compute pseudo-bulks
    cs <- .split_cells(x, by)
    pb <- .pb(x, cs, assay, fun)
    if (scale) {
        if (assay = "countts" & fun == "rowSums") {
            pb_sum <- pb
        } else {
            pb_sum <- .pb(x, cs, "counts", "rowSums")
        }
        pb <- map_depth(pb_sum, -2, function(u)
            u / colSums(u)[col(u)] * 1e6)
    }

    # construct SCE
    pb <- SingleCellExperiment(pb, metadata = md)
    
    # propagate colData columns that are unique across 2nd 'by'
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
    } else {
        colData(pb) <- DataFrame()
    }
    return(pb)
}
