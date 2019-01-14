#' calcExprFreq
#' 
#' Calculates gene expression frequencies
#' 
#' \code{calcExprFreq} computes for each cluster-sample, the fraction of cells
#' that express a given gene. Here, a gene is considered to be expressed when
#' the specified measurement value lies above the specified threshold value.
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param assay a character string specifying which assay to use.
#' @param th numeric threshold value above which
#'   a gene should be considered to be expressed.
#' @param melt logical. Should results be metled into a single data.frame.
#' 
#' @return returns a list of length #clusters containing, for each cluster, 
#'   a data.frame of dimensions #genes x #samples giving the fraction of cells
#'   that express each gene in each sample. 
#'   If \code{melt = TRUE}, a data.frame with columns \code{"cluster_id"}, 
#'   \code{"sample_id"}, and \code{"freq"} will be returned.
#'   
#' @examples
#' data(kang)
#' freqs <- calcExprFreq(kang)
#' names(freqs)
#' head(freqs[[1]])
#' 
#' freqs <- calcExprFreq(kang, melt = TRUE)
#' head(freqs)
#' 
#' @author Helena L. Crowell \email{helena.crowell@uzh.ch} and Mark D. Robinson.
#' 
#' @importFrom data.table data.table
#' @importFrom Matrix rowMeans
#' @importFrom methods is
#' @importFrom reshape2 melt
#' @importFrom SummarizedExperiment assays colData
#' 
#' @export
    
calcExprFreq <- function(x, assay = "counts", th = 0, melt = FALSE) {
    
    # check validity of input arguments
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot(all(c("cluster_id", "sample_id") %in% colnames(colData(x))))
    stopifnot(is.character(assay), length(assay) == 1, assay %in% names(assays(x)))
    stopifnot(is.numeric(th), length(th) == 1)
    stopifnot(is.logical(melt), length(melt) == 1)
    
    # split cells by cluster-sample
    dt <- data.table(
        cell = colnames(x),
        cluster_id = colData(x)$cluster_id,
        sample_id = colData(x)$sample_id)
    dt_split <- split(dt, 
        by = c("cluster_id", "sample_id"), 
        sorted = TRUE, flatten = FALSE)
    cells_by_cluster_sample <- lapply(dt_split, lapply, "[[", "cell")
    
    # for each gene, compute fraction of cells 
    # w/ assay value above threshold in each cluster-sample
    freqs <- lapply(cells_by_cluster_sample, sapply, function(i)
        rowMeans(assays(x)[[assay]][, i, drop = FALSE] > th))
    
    if (melt) {
        df <- data.frame(
            check.names = FALSE,
            do.call(rbind, freqs),
            cluster_id = rep(names(freqs), nrow(x)))
        freqs <- melt(df, 
            id.var = "cluster_id", 
            variable.name = "sample_id",
            value.name = "freq")
    }
    return(freqs)
}
