#' calcExprFreqs
#' 
#' Calculates gene expression frequencies
#' 
#' \code{calcExprFreq} computes, for each e.g. cluster, cluster-sample, 
#' or cluster-group, the fraction of cells that express a given gene. 
#' Here, a gene is considered to be expressed when the specified 
#' measurement value lies above the specified threshold value.
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param assay a character string specifying which assay to use.
#' @param th numeric threshold value above which
#'   a gene should be considered to be expressed.
#' 
#' @return returns a list of length #clusters containing, for each cluster, 
#'   a data.frame of dimensions #genes x #samples giving the fraction of cells
#'   that express each gene in each sample. 
#'   If \code{melt = TRUE}, a data.frame with columns \code{"cluster_id"}, 
#'   \code{"sample_id"}, and \code{"freq"} will be returned.
#'   
#' @examples
#' data(kang)
#' frqs <- calcExprFreqs(kang, by = "cluster_id")
#' head(frqs)
#' 
#' frqs <- calcExprFreqs(kang, by = c("cluster_id", "sample_id"))
#' names(frqs)
#' head(frqs[[1]])
#' 
#' @author Helena L. Crowell \email{helena.crowell@uzh.ch} and Mark D. Robinson.
#' 
#' @importFrom data.table data.table
#' @importFrom dplyr bind_cols bind_rows
#' @importFrom Matrix rowMeans
#' @importFrom methods is
#' @importFrom SummarizedExperiment assayNames assays colData
#' 
#' @export
    
calcExprFreqs <- function(x, assay = "counts", th = 0) {
    
    # check validity of input arguments
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot(is.character(assay), length(assay) == 1, assay %in% assayNames(x))
    stopifnot(is.numeric(th), length(th) == 1)
    
    # split cells according to 'by'
    dt <- data.table(data.frame(colData(x)), cell = colnames(x))
    dt_split <- split(dt, 
        by = c("cluster_id", "sample_id"), 
        sorted = TRUE, flatten = FALSE)
    cells <- lapply(dt_split, lapply, "[[", "cell")

    # for each gene, compute fraction of cells 
    # w/ assay value above threshold in each group of cells
    data <- as.matrix(assays(x)[[assay]])
    fq_by_sample <- lapply(cells, vapply, function(i)
        rowMeans(data[, i, drop = FALSE] > th),
        numeric(nrow(x)))

    cluster_ids <- colData(x)$cluster_id
    n_cells_by_sample <- table(cluster_ids, colData(x)$sample_id)
    n_cells_by_group <- table(cluster_ids, colData(x)$group_id)
    
    ei <- metadata(x)$experiment_info
    samples_by_group <- split(levels(colData(x)$sample_id), ei$group_id) 
    groups <- names(samples_by_group)
    
    cluster_ids <- levels(cluster_ids)
    names(cluster_ids) <- cluster_ids
    
    lapply(cluster_ids, function(k) {
        ns_by_sample <- t(t(fq_by_sample[[k]]) * n_cells_by_sample[k, ])
        fq_by_group <- vapply(groups, function(g)
            rowSums(ns_by_sample[, samples_by_group[[g]]]) / n_cells_by_group[k, g], 
            numeric(nrow(x)))
        data.frame(fq_by_group, fq_by_sample[[k]], check.names = FALSE)
    })
}
