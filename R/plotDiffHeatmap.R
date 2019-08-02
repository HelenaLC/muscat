#' @rdname plotDiffHeatmap
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
#' @param assay a character string specifying which assay 
#'   in \code{assays(x)}to obtain expression values from.
#' @param normalize logical specifying whether
#'   mean-expression values be z-normazlied.
#' @param colors character vector of colors to use for plotting.
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
#' @author Helena L Crowell & Mark D Robinson
#' 
#' @import ComplexHeatmap
#' @importFrom dplyr %>% bind_rows filter_
#' @importFrom grDevices colorRampPalette
#' @importFrom grid gpar
#' @importFrom magrittr set_rownames
#' @importFrom methods is
#' @importFrom purrr modify_depth
#' @importFrom scales hue_pal
#' @importFrom SummarizedExperiment assayNames assays colData
#' @importFrom viridis viridis
#' @export

plotDiffHeatmap <- function(x, y, c = NULL, g = NULL, k = NULL, 
    top_n = 20, sort_by = "p_adj.loc", decreasing = FALSE, fdr = 0.05, lfc = 1,
    assay = "logcounts", normalize = TRUE, colors = viridis(10)) {
    
    # validity checks of input arguments
    .check_sce(x)
    .check_res(x, y)
    .check_arg_assay(x, assay)
    stopifnot(is.null(c) || c %in% names(y$table))
    stopifnot(is.null(g) || is.character(g) & all(g %in% rownames(x)))
    stopifnot(is.null(k) || is.character(k) & all(k %in% levels(x$cluster_id)))
    stopifnot(is.numeric(top_n), length(top_n) == 1, top_n > 0)
    
    # default to 1st contrast/coef
    y <- y$table
    if (is.null(c)) 
        c <- names(y)[1]
    y <- y[[c]]
    y <- y[!vapply(y, is.null, logical(1))]
    
    stopifnot(is.character(sort_by), 
        sort_by %in% names(y[[1]]),
        is.numeric(y[[1]][[sort_by]]))

    # filter results
    if (!is.null(g))
        y <- lapply(y, filter_, ~gene %in% g)
    if (is.null(k))
        k <- levels(sce$cluster_id)
    y <- y[k]
    y <- lapply(y, filter_, ~p_adj.loc < fdr, ~abs(logFC) > lfc)
    
    # get cluster IDs & nb. of clusters
    kids <- names(y)
    names(kids) <- kids
    nk <- length(kids)
    
    # order & subset top_n results
    if (is.null(top_n)) {
        ns <- vapply(y, nrow, numeric(1))
    } else {
        ns <- rep(top_n, length(y))
    }
    names(ns) <- kids
    y <- lapply(kids, function(k) {
        u <- y[[k]]
        o <- order(u[[sort_by]], decreasing = decreasing)
        o <- o[seq_len(ns[k])]
        u[o[!is.na(o)], ]
    }) %>% bind_rows
    
    es <- assays(x)[[assay]]
    es <- es[unlist(y$gene), , drop = FALSE]
    
    # split cells by cluster-sample
    cells_by_ks <- .split_cells(x)
    
    # compute cluster-sample means
    ms <- t(vapply(seq_len(nrow(y)), function(i) {
        g <- y$gene[i]
        k <- y$cluster_id[i]
        vapply(cells_by_ks[[k]], function(j)
            mean(es[g, j]), numeric(1))
    }, numeric(nlevels(x$sample_id)))) %>% 
        set_rownames(y$gene)
    if (normalize) ms <- .z_norm(ms)
    
    # row annotation
    if (length(kids) > 1) {
        if (nk > length(cluster_colors))
            cluster_colors <- colorRampPalette(cluster_colors)(nk)
        cols <- setNames(cluster_colors[seq_len(nk)], kids)
        row_anno <- data.frame(cluster_id = y$cluster_id)
        row_anno <- rowAnnotation(row_anno,
            col = list(cluster_id = cols),
            gp = gpar(col = "white"))
    } else {
        row_anno <- NULL
    }
    
    # column annoation
    ei <- metadata(x)$experiment_info
    m <- match(levels(x$sample_id), ei$sample_id)
    cols <- setNames(hue_pal()(nlevels(x$group_id)), levels(x$group_id))
    col_anno <- data.frame(group_id = ei$group_id[m])
    col_anno <- columnAnnotation(col_anno,
        col = list(group_id = cols),
        gp = gpar(col = "white"))
    
    main <- sprintf("%s %s [%s%s]", c, 
        ifelse(length(k) == 1, sprintf("(%s)", k), ""), 
        ifelse(is.null(top_n), "", sprintf("top_n = %s, ", top_n)), 
        paste("sort_by =", dQuote(sort_by)))
    
    row_anno + Heatmap(ms,
        col = colors,
        name = paste0("z-normalized\n"[normalize], "mean expr."),
        column_title = main,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 6),
        top_annotation = col_anno,
        split = y$cluster_id,
        combined_name_fun = NULL)
}
