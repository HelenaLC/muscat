#' @rdname pbMDS
#' @title Pseudobulk-level MDS plot
#' 
#' @description Renders a multidimensional scaling (MDS) 
#'   where each point represents a cluster-sample instance; 
#'   with points colored by cluster ID and shaped by group ID. 
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}} 
#'   containing cluster-sample pseudobulks as returned by 
#'   \code{\link{aggregateData}} with argument 
#'   \code{by = c("cluster_id", "sample_id")}.
#' 
#' @return a \code{ggplot} object.
#' 
#' @examples 
#' data(example_sce)
#' pb <- aggregateData(example_sce)
#' pbMDS(pb)
#' 
#' @author Helena L Crowell & Mark D Robinson
#' 
#' @importFrom edgeR calcNormFactors cpm DGEList plotMDS.DGEList
#' @importFrom SummarizedExperiment assays
#' @importFrom grDevices colorRampPalette
#' @importFrom S4Vectors metadata
#' @importFrom Matrix rowSums
#' @importFrom rlang .data
#' @import ggplot2
#' @export

pbMDS <- function(x) {
    # check validity of input pseudobulk-SCE
    # (cells should have been aggregated by cluster-sample)
    .check_pbs(x, check_by = TRUE)
    
    y <- as.list(assays(x))
    y <- do.call("cbind", y)
    y <- y[, (j <- c(t(.n_cells(x))) != 0)]
    d <- DGEList(unname(y), remove.zeros = TRUE)
    d <- calcNormFactors(d)
    
    mds <- plotMDS.DGEList(d, plot = FALSE)
    nk <- length(kids <- assayNames(x))

    ss <- rep(colnames(x), nk)
    ks <- rep(kids, each = ncol(x))
    
    if (any(!j)) {
        txt <- paste(sQuote(ks[!j]), sQuote(ss[!j]), sep = "-")
        message("Removing cluster-sample instance(s) ", 
            paste(txt, collapse = ", "))
    }

    df <- data.frame(
        MDS1 = mds$x, MDS2 = mds$y, 
        cluster_id = factor(ks[j], levels = kids), 
        group_id = rep(x$group_id, nk)[j])
    
    cols <- .cluster_colors
    if (nk > length(cols)) 
        cols <- colorRampPalette(cols)(nk)
    
    ggplot(df, aes(.data$MDS1, .data$MDS2, 
        col=.data$cluster_id, shape=.data$group_id)) +
        scale_color_manual(values = cols) +
        geom_point(size = 3, alpha = 0.8) + 
        guides(color = guide_legend(override.aes = list(alpha = 1))) +
        theme_bw() + theme(aspect.ratio = 1,
            axis.text = element_text(color = "black"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(linewidth = 0.2, color = "lightgrey"))
}
