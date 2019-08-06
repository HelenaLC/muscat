#' @rdname pbHeatmap
#' @title Heatmap of cluster-sample pseudobulks
#' 
#' @description ...
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param y a list of DS analysis results as returned by 
#'   \code{\link{pbDS}} or \code{\link{mmDS}}.
#' @param k character vector; specifies which cluster ID(s) to retain.
#'   Defaults to \code{levels(x$cluster_id)}.
#' @param g character vector; specifies which genes to retain.
#'   Defaults to considering all genes.
#' @param c character string; specifies which contrast/coefficient to retain.
#'   Defaults to \code{names(y$table)[1]}.
#' @param top_n single numeric; number of genes to retain per cluster.
#' @param fdr,lfc single numeric; FDR and logFC cutoffs to filter results by.
#'   The specified FDR threshold is applied to \code{p_adj.loc} values.
#' @param sort_by character string specifying 
#'   a numeric results table column to sort by.
#' @param decreasing logical; whether to sort 
#'   in decreasing order of \code{sort_by}.
#' @param assay character string; specifies which assay to use;
#'   should be one of \code{assayNames(x)}.
#' @param fun function to use as summary statistic, 
#'   e.g., mean, median, sum (depending on the input assay).
#' @param normalize logical; whether to apply a z-normalization 
#'   to each row (gene) of the cluster-sample pseudobulk data.
#' @param col character vector of colors or color mapping function
#'   generated with \code{\link[circlize]{colorRamp2}}. Passed to 
#'   argument \code{col} in \code{\link[ComplexHeatmap]{Heatmap}}
#'   (see \code{?ComplexHeatmap::Heatmap} for details).
#' @param row_anno,col_anno logical; whether to render
#'   annotations of cluster and group IDs, respectively.
#' 
#' @return a \code{\link{HeatmapList-class}} object.
#' 
#' @examples 
#' data(sce)
#' 
#' # compute pseudobulks & run DS analysis
#' pb <- aggregateData(sce)
#' res <- pbDS(pb)
#' 
#' # cluster-sample expression means
#' pbHeatmap(sce, res)
#' 
#' # include only a single cluster
#' pbHeatmap(sce, res, k = levels(sce$cluster_id)[1])
#' 
#' # plot specific gene across all clusters
#' pbHeatmap(sce, res, g = rownames(sce)[1])
#' 
#' # show sum of counts instead
#' pbHeatmap(sce, res, assay = "counts", fun = sum)
#' 
#' @author Helena L Crowell
#' 
#' @importFrom ComplexHeatmap Heatmap columnAnnotation rowAnnotation
#' @importFrom dplyr bind_rows filter_
#' @importFrom grid gpar
#' @importFrom purrr map
#' @importFrom scales hue_pal
#' @importFrom viridis viridis
#' @export

pbHeatmap <- function(x, y, 
    k = NULL, g = NULL, c = NULL, 
    top_n = 20, fdr = 0.05, lfc = 1, 
    sort_by = "p_adj.loc", decreasing = FALSE,
    assay = "logcounts", fun = mean, normalize = TRUE,
    col = viridis(10), row_anno = TRUE, col_anno = TRUE) {
    
    # check validity of input arguments
    .check_sce(x, req_group = TRUE)
    .check_arg_assay(x, assay)
    .check_args_pbHeatmap(as.list(environment()))
    #.check_res(x, y)
    
    # defaults for NULL arguments
    if (is.null(k)) k <- levels(x$cluster_id)
    if (is.null(c)) c <- names(y$table)[1]
    
    # subset specified contrast/coef & cluster(s)
    y <- y$table[[c]][k]
    y <- y[!vapply(y, is.null, logical(1))]
    
    # filter results
    if (!is.null(g)) y <- lapply(y, filter_, ~gene %in% g)
    y <- lapply(y, filter_, ~p_adj.loc < fdr, ~abs(logFC) > lfc)
    
    # get cluster IDs & nb. of clusters
    nk <- length(names(kids) <- kids <- names(y))
    
    # subset 'top_n' results
    if (is.null(top_n)) {
        ns <- vapply(y, nrow, numeric(1))
    } else {
        ns <- vapply(y, function(u) min(nrow(u), top_n), numeric(1))
    }
    
    # re-order results
    vs <- map(y, sort_by)
    os <- lapply(vs, order, decreasing = decreasing)
    y <- lapply(kids, function(k) {
        o <- os[[k]][seq_len(ns[k])]
        y[[k]][o, ]
    })
    y <- bind_rows(y)
    
    # subset 'assay' data
    es <- assays(x)[[assay]]
    es <- es[unlist(y$gene), , drop = FALSE]
    
    # compute cluster-sample summary values (e.g., means)
    cells_by_ks <- .split_cells(x)
    xs <- t(mapply(function(g, k)
        vapply(cells_by_ks[[k]], function(cs)
            fun(es[g, cs]), numeric(1)), 
        g = y$gene, k = y$cluster_id))
    if (normalize) xs <- .z_norm(xs)
    
    # plotting -----------------------------------------------------------------
    # row & column annotation
    lgd_aes <- list(labels_gp = gpar(fontsize = 6),
        title_gp = gpar(fontface = "bold", fontsize = 8))
    
    if (row_anno & nk > 1) {
        cols <- .cluster_colors
        if (nk > length(cols))
            cols <- colorRampPalette(cols)(nk)
        cols <- setNames(cols[seq_len(nk)], kids)
        row_anno <- rowAnnotation(
            df = data.frame(cluster_id = y$cluster_id),
            col = list(cluster_id = cols),
            gp = gpar(col = "white"),
            show_annotation_name = FALSE,
            annotation_legend_param = lgd_aes)
    } else {
        row_anno <- NULL
    }
    
    ei <- metadata(x)$experiment_info
    m <- match(levels(x$sample_id), ei$sample_id)
    if (col_anno) {
        cols <- setNames(hue_pal()(nlevels(x$group_id)), levels(x$group_id))
        col_anno <- columnAnnotation(
            df = data.frame(group_id = ei$group_id[m]),
            col = list(group_id = cols),
            gp = gpar(col = "white"),
            show_annotation_name = FALSE,
            annotation_legend_param = lgd_aes)
    } else {
        col_anno <- NULL
    }
    
    Heatmap(matrix = xs, col = col, 
        name = sprintf("%s%s %s", 
            c("", "z-normalized\n")[normalize + 1], 
            deparse(substitute(fun)), assay),
        row_title = NULL, column_title = NULL,
        cluster_rows = FALSE, cluster_columns = FALSE,
        left_annotation = row_anno, top_annotation = col_anno,
        split = y$cluster_id, column_split = ei$group_id[m],
        heatmap_legend_param = lgd_aes,
        row_names_gp = gpar(fontsize = 6), 
        column_names_gp = gpar(fontsize = 8),
        column_title_gp = gpar(fontface = "bold", fontsize = 10))
}