#' @rdname plotDiffGenes
#' @title Heatmap of mean-marker expression by cluster-sample
#' 
#' @description ...
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param y a list of differential testing results 
#'   as returned by \code{\link{runDS}}.
#' @param c a character string specifying the comparison 
#'   (contrast of coefficient) for which results should be plotted.
#'   Defaults to the first one available.
#' @param top_n single numeric specifying the number of genes to include.
#' @param sort_by a character string specifying a \code{y} column to sort by.
#' @param clusters a character string specifying which cluster(s) results
#'   should be included for. If NULL (the default), the \code{top_n} hits
#'   for each cluster will be included.
#' @param fdr single numeric specifying the threshold on adjusted p-values
#'   below which results should be considered significant.
#' @param lfc single numeric specifying the threshold on absolute logFCs
#'   above which results should be included.
#' 
#' @return a \code{\link{HeatmapList-class}} object.
#' 
#' @author Helena L. Crowell \email{helena@crowells.eu}
#' 
#' @import ComplexHeatmap
#' @importFrom dplyr %>% bind_rows filter
#' @importFrom grid gpar
#' @importFrom methods is
#' @importFrom purrr modify_depth
#' @importFrom SummarizedExperiment assayNames assays colData
#' @importFrom viridis viridis
#' @export

plotDiffGenes <- function(x, y, c = NULL, g = NULL, k = NULL, 
    top_n = 20, sort_by = c("FDR", "logFC"), fdr = 0.05, lfc = 1) {
    
    # validity checks of input arguments
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot("logcounts" %in% assayNames(x))
    stopifnot(is.null(c) || c %in% names(y$table))
    stopifnot(is.numeric(top_n), length(top_n) == 1, top_n > 0)
    sort_by <- match.arg(sort_by)
    
    # default to 1st contrast/coef
    y <- y$table
    if (is.null(c)) c <- names(y)[1]
    y <- y[[c]]
    y <- y[!sapply(y, is.null)]
    
    # filter results
    if (!is.null(g)) {
        y <- lapply(y, filter, gene %in% g)
    }
    if (!is.null(k)) {
        if (k == "all") {
            cluster_ids <- colData(x)$cluster_id
            k <- levels(cluster_ids)
        }
        y <- y[k]
    }
    y <- lapply(y, filter, FDR < fdr, abs(logFC) > lfc)
    
    # order & subset top_n results
    y <- switch(sort_by,
        FDR = bind_rows(lapply(y, function(u) {
            o <- order(u$FDR)
            o <- o[seq_len(top_n)]
            u[o[!is.na(o)], ]
        }))
        ,
        logFC = bind_rows(lapply(y, function(u) {
            o <- order(u$logFC, decreasing = TRUE)
            o <- o[seq_len(top_n)]
            u[o[!is.na(o)], ]
        }))
    )
    
    es <- as.matrix(logcounts(x)[unlist(y$gene), ])
    es0 <- t(CATALYST:::scale_exprs(t(es)))
    
    # split cells by cluster-sample
    cells_by_cluster_sample <- .split_cells(x)
    
    ms <- t(apply(y, 1, function(u) {
        g <- u[["gene"]]
        k <- u[["cluster_id"]]
        vapply(cells_by_cluster_sample[[k]], 
            function(i) mean(es0[g, i]), numeric(1))
    }))
    rownames(ms) <- sprintf("%s(%s)", y$gene, y$cluster_id)
    colnames(ms) <- levels(colData(x)$sample_id)
    
    Heatmap(ms,
        col = viridis(10),
        name = "avg. scaled\nexpression",
        column_title = sprintf("%s (top %s %s)", c, top_n, sort_by),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 6))
}
