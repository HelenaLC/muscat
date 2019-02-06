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
#' @importFrom SummarizedExperiment assays colData SummarizedExperiment
#' @export
    
calcExprFreqs <- function(x, assay = "counts", th = 0) {
    
    # check validity of input arguments
    .check_sce(x)
    .check_arg_assay(x, assay)
    stopifnot(is.numeric(th), length(th) == 1)
    
    # split cells by cluster-sample
    cells_by_cluster_sample <- .split_cells(x)

    # for each gene, compute fraction of cells 
    # w/ assay value above threshold in each sample
    fun <- getFromNamespace("rowMeans", "Matrix")
    fq <- lapply(cells_by_cluster_sample, vapply, function(i) 
        fun(assays(x)[[assay]][, i, drop = FALSE] > th),
        numeric(nrow(x)))

    if ("group_id" %in% colnames(colData(x))) {
        kids <- colData(x)$cluster_id
        sids <- colData(x)$sample_id
        gids <- colData(x)$group_id
        
        n_cells_by_ks <- table(kids, sids)
        n_cells_by_kg <- table(kids, gids)
        
        ei <- metadata(x)$experiment_info
        samples_by_g <- split(levels(sids), ei$group_id) 
        gids <- names(samples_by_g)
        
        kids <- levels(kids)
        names(kids) <- kids
        
        fq <- lapply(kids, function(k) {
            n_cells_s <- fq[[k]] * n_cells_by_ks[k, ][col(fq[[k]])]
            fq_by_g <- vapply(gids, function(g) {
                n_cells_g <- n_cells_s[, samples_by_g[[g]], drop = FALSE]
                rowSums(n_cells_g) / n_cells_by_kg[k, g]
            }, numeric(nrow(x)))
            cbind(fq[[k]], fq_by_g)
        })
    }
    # return SCE
    SummarizedExperiment(assays = fq)
}
