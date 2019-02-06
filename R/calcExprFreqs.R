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
#' @return a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   containing, for each cluster, an assay of dimensions #genes x #samples 
#'   giving the fraction of cells that express each gene in each sample. 
#'   If \code{colData(x)} contains a \code{"group_id"} column, the fraction 
#'   of expressing cells in each each group will be included as well.
#'   
#' @examples
#' library(SummarizedExperiment)
#' 
#' colnames(colData(kang)) # contains sample_id only
#' frq <- calcExprFreqs(kang)
#' assayNames(frq)  # one assay per cluster
#' head(assay(frq)) # expr. freqs. by sample
#' 
#' sim <- simData(kang, n_genes = 5, n_cells = 1e3) 
#' frq <- calcExprFreqs(sim)
#' colnames(colData(sim)) # contains group_id!
#' assay(frq)             # freqs. by group included
#' 
#' @author Helena L. Crowell \email{helena.crowell@uzh.ch} and Mark D. Robinson.
#' 
#' @importFrom Matrix rowMeans
#' @importFrom methods is
#' @importFrom SummarizedExperiment assays colData
#' 
#' @export
    
calcExprFreqs <- function(x, assay = "counts", th = 0) {
    
    # check validity of input arguments
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot(c("cluster_id", "sample_id") %in% colnames(colData(x)))
    stopifnot(is.numeric(th), length(th) == 1)
    .check_arg_assay(x, assay)
    
    # split cells by cluster-sample
    cells_by_cluster_sample <- .split_cells(x)

    # for each gene, compute fraction of cells 
    # w/ assay value above threshold in each sample
    data <- as.matrix(assays(x)[[assay]])
    fq <- lapply(cells_by_cluster_sample, vapply, 
        function(i) rowMeans(data[, i, drop = FALSE] > th),
        numeric(nrow(x)))

    if ("group_id" %in% colnames(colData(x))) {
        cluster_ids <- colData(x)$cluster_id
        n_cells_by_s <- table(cluster_ids, colData(x)$sample_id)
        n_cells_by_g <- table(cluster_ids, colData(x)$group_id)
        
        ei <- metadata(x)$experiment_info
        sample_by_g <- split(levels(colData(x)$sample_id), ei$group_id) 
        groups <- names(sample_by_g)
        
        cluster_ids <- levels(cluster_ids)
        names(cluster_ids) <- cluster_ids
        
        fq <- lapply(cluster_ids, function(k) {
            n_cells_s <- fq[[k]] *  n_cells_by_s[k, ][col(fq[[k]])]
            fq_by_g <- vapply(groups, function(g) {
                n_cells_g <- n_cells_s[, sample_by_g[[g]]]
                rowSums(n_cells_g) / n_cells_by_g[k, g]
            }, numeric(nrow(x)))
            cbind(fq[[k]], fq_by_g)
        })
    }
    # return SCE
    SummarizedExperiment(assays = fq)
}
