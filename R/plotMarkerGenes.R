#' @rdname plotMarkerGenes
#' @title Heatmap of marker-gene expression by cluster
#' 
#' @description ...
#' 
#' @param x a \code{[SingleCellExperiment]{SingleCellExperiment}}.
#' @param marker_genes a names list of marker genes for each cluster.
#' @param scale logical. 
#'   Specifies wether scaled expression values should be scaled 
#'   between 0 and 1 using 1% and 99% quantiles as boundaries.
#' @param cluster_columns logical. 
#'   Should columns be hierarchically clustered using euclidean distances.
#' 
#' @return a \code{\link{HeatmapList-class}} object.
#' 
#' @author Helena L. Crowell \email{helena@crowells.eu}
#' 
#' @import ComplexHeatmap
#' @import SingleCellExperiment
#' @importFrom dplyr %>% group_by summarise_all
#' @importFrom grid gpar
#' @importFrom viridis viridis
#' @export

plotMarkerGenes <- function(x, marker_genes, scale = TRUE, cluster_columns = TRUE) {
    cluster_ids <- colData(x)$cluster_id
    es <- logcounts(x)[unlist(marker_genes), ]
    es <- t(as.matrix(es))
    if (scale)
        es <- CATALYST:::scale_exprs(es)
    df <- data.frame(es, cluster_id = cluster_ids)
    means_by_cluster <- df %>% 
        group_by(.data$cluster_id) %>% 
        summarise_all(mean) %>% 
        data.frame(row.names = 1)
    mat <- as.matrix(means_by_cluster)
    mat <- mat[, !is.na(mat[1, ])]
    row_anno <- levels(cluster_ids)
    if (cluster_columns) {
        hm <- Heatmap(mat)
        row_o <- row_order(hm)[[1]]
        col_o <- unlist(marker_genes[row_o])
        col_o <- make.names(col_o)
        col_o <- col_o[!duplicated(col_o)]
        col_o <- col_o[col_o %in% colnames(mat)]
        mat <- mat[row_o, col_o]
        row_anno <- row_anno[row_o]
    }
    row_anno <- Heatmap(
        mat = row_anno,
        col = setNames(CATALYST:::cluster_cols, levels(cluster_ids)),
        name = "cluster_id",
        cluster_rows = FALSE,
        rect_gp = gpar(col = "white"))
    hm <- Heatmap(
        matrix = mat, 
        col = viridis(10),
        name = sprintf("mean%s\nexpression", ifelse(scale, " scaled", "")),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 8))
    row_anno + hm
}
