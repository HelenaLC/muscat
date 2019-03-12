#' @rdname pbMDS
#' @title Pseudo-bulk level MDS plot
#' 
#' @description ...
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}} 
#'   containing pseudo-bulk data as returned by \code{\link{aggregateData}}.
#' 
#' @return a \code{ggplot} object.
#' 
#' @examples 
#' sce <- toyData()
#' sce <- prepData(sce, "cluster_id", "sample_id", "group_id")
#' pb <- aggregateData(sce)
#' pbMDS(pb)
#' 
#' @author Helena L. Crowell \email{helena.crowells@uzh.ch} and Mark D. Robinson.
#' 
#' @import ggplot2
#' @importFrom edgeR calcNormFactors cpm DGEList plotMDS.DGEList
#' @importFrom dplyr bind_cols
#' @importFrom grDevices colorRampPalette
#' @importFrom Matrix rowSums
#' @importFrom methods as
#' @importFrom SummarizedExperiment assays
#' @importFrom S4Vectors metadata
#' @export

pbMDS <- function(x) {
    y <- as(assays(x), "list")
    y <- do.call("cbind", y)
    d <- suppressMessages(DGEList(y))
    d <- calcNormFactors(d)
    
    mds <- plotMDS.DGEList(d, plot = FALSE)
    ei <- metadata(x)$experiment_info
    m <- match(colnames(x), ei$sample_id)
    nk <- length(assays(x))
    df <- data.frame(
        MDS1 = mds$x, 
        MDS2 = mds$y, 
        cluster_id = rep(assayNames(x), each = ncol(x)),
        group_id = rep(ei$group_id[m], nk))
    
    cols <- cluster_colors
    if (nk > length(cols)) 
        cols <- colorRampPalette(cols)(nk)
    
    ggplot(df, aes_string(x = "MDS1", y = "MDS2", 
        col = "cluster_id", shape = "group_id")) +
        scale_color_manual(values = cols) +
        geom_point(size = 5, alpha = .8) + 
        guides(color = guide_legend(override.aes = list(alpha = 1))) +
        theme_bw() + theme(aspect.ratio = 1,
            axis.text = element_text(color = "black"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(size = 0.2, color = "lightgrey"))
}
