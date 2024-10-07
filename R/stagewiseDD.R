#' @importFrom dplyr bind_rows
.res_DX <- function(res_DS, res_DD) {
    # for each contrast...
    names(cts) <- cts <- names(res_DS$table) 
    lapply(cts, function(ct) { 
        # for each cluster...
        names(kids) <- kids <- names(res_DS$table[[ct]]) 
        lapply(kids, function(kid) { 
            # get DS/DD results
            DS <- res_DS$table[[ct]][[kid]]
            DD <- res_DD$table[[ct]][[kid]]
            # add missing genes
            DS <- bind_rows(DS, data.frame(gene=setdiff(DD$gene, DS$gene)))
            DD <- bind_rows(DD, data.frame(gene=setdiff(DS$gene, DD$gene)))
            # reorder & return both
            DD <- DD[match(DS$gene, DD$gene), ]
            return(list(DS=DS, DD=DD))
        })
    })
}

#' @rdname stagewise_DS_DD
#' @title Perform two-stage testing on DS and DD analysis results
#'
#' @param res_DS a list of DS testing results as returned 
#'   by \code{\link{pbDS}} or \code{\link{mmDS}}.
#' @param res_DD a list of DD testing results as returned 
#'   by \code{\link{pbDD}} (or \code{\link{pbDS}} with \code{method="DD"}).
#' @param sce (optional) \code{SingleCellExperiment} object containing the data 
#'   that underlies testing, prior to summarization with \code{\link{aggregateData}}. 
#'   Used for validation of inputs in order to prevent unexpected failure/results.
#' @param verbose logical. Should information on progress be reported?
#'
#' @return 
#' A list of \code{DFrame}s containing results for each contrast and cluster.
#' Each table contains DS and DD results for genes shared between analyses,
#' as well as results from stagewise testing analysis, namely:
#' \itemize{
#' \item{\code{p_adj}: FDR adjusted p-values for the
#'   screening hypothesis that a gene is neither DS nor DD
#'   (see \code{?stageR::getAdjustedPValues} for details)}
#' \item{\code{p_val.DS/D}: confirmation stage p-values for DS/D}}
#'
#' @examples
#' data(example_sce)
#' 
#' pbs_sum <- aggregateData(example_sce, assay="counts", fun="sum")
#' pbs_det <- aggregateData(example_sce, assay="counts", fun="num.detected")
#'   
#' res_DS <- pbDS(pbs_sum, min_cells=0, filter="none", verbose=FALSE)
#' res_DD <- pbDD(pbs_det, min_cells=0, filter="none", verbose=FALSE)
#' 
#' res <- stagewise_DS_DD(res_DS, res_DD)
#' head(res[[1]][[1]]) # results for 1st cluster
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom purrr map_depth
#' @export

stagewise_DS_DD <- function(res_DS, res_DD, sce=NULL, verbose=FALSE) {
    if (!requireNamespace("stageR", quietly=TRUE))
        stop("Install 'stageR' to use this function.")
        
    # validity checks
    # TODO: helper to check validity of 'res_DS/D'
    # against each other and, optionally, 'sce'
    stopifnot(
        # same coefs/constrasts
        names(x <- res_DS$table) == 
        names(y <- res_DD$table), 
        # any shared clusters
        sum(mapply(\(i, j) 
            length(intersect(i, j)), 
            i=lapply(x, names), 
            j=lapply(y, names))) > 0)
    if (!is.null(sce)) {
        .check_sce(sce)
        . <- map_depth(list(x, y), 3, \(df) df$gene %in% rownames(sce))
        stopifnot("gene(s) present in 'res_DS/D' not found in 'sce'"=unlist(.))
        . <- map_depth(list(x, y), 3, \(df) df$cluster_id %in% sce$cluster_id)
        stopifnot("cluster(s) present in 'res_DS/D' not found in 'sce'"=unlist(.))
    }
    
    # assure that results contain same set of genes, in the same order
    # (indepedent of different filtering criteria for the two analyses)
    res_DX <- .res_DX(res_DS=res_DS, res_DD=res_DD)
    
    # perform harmonic mean p-value aggregation according to
    # (https://www.pnas.org/doi/full/10.1073/pnas.1814092116)
    .mu <- \(x) 1/mean(1/x, na.rm=TRUE)
    
    # perform stagewise testing
    res <- map_depth(res_DX, 2, \(x) {
        ps <- data.frame(
            p_val.DS=x$DS$p_val,
            p_val.DD=x$DD$p_val,
            row.names=x$DS$gene)
        qs <- apply(ps, 1, .mu); names(qs) <- x$DS$gene
        obj <- stageR::stageR(qs, as.matrix(ps), FALSE)
        eva <- expression({
            obj <- stageR::stageWiseAdjustment(obj, 
                method="none", alpha=0.05, allowNA=TRUE)
            res <- stageR::getAdjustedPValues(obj, 
                onlySignificantGenes=FALSE, order=FALSE)
        })
        res <- if (verbose) eval(eva) else suppressMessages({eval(eva)})
        # TODO: communicate this better with the user?
        if (is.null(res)) res <- NA
        colnames(res)[1] <- "p_adj"
        return(res)
    })
    names(cs) <- cs <- names(res)
    lapply(cs, \(c) {
        names(ks) <- ks <- names(res[[c]])
        lapply(ks, \(k) {
            gs <- rownames(df <- res[[c]][[k]])
            res_DS <- I((. <- res_DS$table[[c]][[k]])[match(gs, .$gene), ])
            res_DD <- I((. <- res_DD$table[[c]][[k]])[match(gs, .$gene), ])
            DataFrame(gene=gs, df, cluster_id=k, contrast=c, res_DS, res_DD)
        })
    })
}
