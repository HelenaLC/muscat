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
#' @param g,k character strings specifying the gene and cluster, 
#'   respectively, to include. Defaults to NULL (all genes/clusters).
#' @param top_n single numeric specifying the number of genes to include.
#' @param sort_by a character string specifying a \code{y} column to sort by.
#' @param fdr single numeric specifying the threshold on adjusted p-values
#'   below which results should be considered significant.
#' @param lfc single numeric specifying the threshold on absolute logFCs
#'   above which results should be included.
#' 
#' @return a \code{\link{HeatmapList-class}} object.
#' 
#' @examples 
#' data(sce)
#' 
#' # compute pseudo-bulk counts
#' pb <- aggregateData(sce)
#' 
#' # specify design & contrast matrix
#' ei <- metadata(sce)$experiment_info
#' design <- model.matrix(~ 0 + ei$group_id)
#' dimnames(design) <- list(ei$sample_id, levels(ei$group_id))
#' contrast <- limma::makeContrasts("stim-ctrl", levels = design)
#' 
#' # test for cluster-specific DE 
#' res <- runDS(sce, pb, design, contrast, method = "edgeR")
#' 
#' @author Helena L. Crowell \email{helena@crowells.eu}
#' 
#' @import ComplexHeatmap
#' @importFrom dplyr %>% bind_rows filter_
#' @importFrom grid gpar
#' @importFrom methods is
#' @importFrom purrr modify_depth
#' @importFrom scales hue_pal
#' @importFrom SummarizedExperiment assayNames assays colData
#' @importFrom viridis viridis
#' @export

plotDiffGenes <- function(x, y, c = NULL, g = NULL, k = NULL, 
    top_n = 20, sort_by = c("p_adj", "logFC"), fdr = 0.05, lfc = 1) {
    
    # validity checks of input arguments
    .check_sce(x)
    .check_res(x, y)
    .check_arg_assay(x, "logcounts")
    stopifnot(is.null(c) || c %in% names(y$table))
    stopifnot(is.null(g) || is.character(g) & all(g %in% rownames(x)))
    stopifnot(is.null(k) || is.character(k) & all(k %in% levels(x$cluster_id)))
    stopifnot(is.numeric(top_n), length(top_n) == 1, top_n > 0)
    sort_by <- match.arg(sort_by)
    
    # default to 1st contrast/coef
    y <- y$table
    if (is.null(c)) c <- names(y)[1]
    y <- y[[c]]
    y <- y[!vapply(y, is.null, logical(1))]
    
    # filter results
    if (!is.null(g)) {
        y <- lapply(y, filter_, ~gene %in% g)
    }
    if (!is.null(k)) {
        if (k == "all") {
            cluster_ids <- colData(x)$cluster_id
            k <- levels(cluster_ids)
        }
        y <- y[k]
    }
    y <- lapply(y, filter_, ~p_adj < fdr, ~abs(logFC) > lfc)
    
    # order & subset top_n results
    y <- switch(sort_by,
        p_adj = bind_rows(lapply(y, function(u) {
            o <- order(u$p_adj)
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
    
    es <- assays(x)$logcounts
    es <- es[unlist(y$gene), ]
    
    # split cells by cluster-sample
    cells_by_ks <- .split_cells(x)
    
    # compute cluster-sample means
    ms <- t(vapply(seq_len(nrow(y)), function(i) {
        g <- y$gene[i]
        k <- y$cluster_id[i]
        vapply(cells_by_ks[[k]], function(j)
            mean(es[g, j]), numeric(1))
    }, numeric(nlevels(x$sample_id))))
    ms <- .scale(ms)
    rownames(ms) <- sprintf("%s(%s)", y$gene, y$cluster_id)
    
    # column annoation
    ei <- metadata(x)$experiment_info
    m <- match(levels(x$sample_id), ei$sample_id)
    cols <- setNames(hue_pal()(nlevels(x$group_id)), levels(x$group_id))
    col_anno <- data.frame(group_id = ei$group_id[m])
    col_anno <- columnAnnotation(col_anno,
        col = list(group_id = cols),
        gp = gpar(col = "white"))
    
    Heatmap(ms,
        col = viridis(10),
        name = "avg. scaled\nexpression",
        column_title = sprintf("%s (top %s %s)", c, top_n, sort_by),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 6),
        top_annotation = col_anno)
}
