#' calcExprFreqs
#' 
#' Calculates gene expression frequencies
#' 
#' \code{calcExprFreq} computes, for each sample and group (in each cluster),
#' the fraction of cells that express a given gene. Here, a gene is considered 
#' to be expressed when the specified measurement value (\code{assay}) 
#' lies above the specified threshold value (\code{th}).
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param assay a character string specifying which assay to use.
#' @param th numeric threshold value above which
#'   a gene should be considered to be expressed.
#' 
#' @return a \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   containing, for each cluster, an assay of dimensions #genes x #samples 
#'   giving the fraction of cells that express each gene in each sample. 
#'   If \code{colData(x)} contains a \code{"group_id"} column, the fraction 
#'   of expressing cells in each each group will be included as well.
#'   
#' @examples
#' data(example_sce)
#' library(SingleCellExperiment)
#' 
#' frq <- calcExprFreqs(example_sce)
#' 
#' # one assay per cluster
#' assayNames(frq)  
#' 
#' # expression frequencies by
#' # sample & group; 1st cluster:
#' head(assay(frq))
#' 
#' @author Helena L Crowell & Mark D Robinson
#' 
#' @importFrom Matrix rowMeans
#' @importFrom methods is
#' @importFrom purrr set_names
#' @importFrom SummarizedExperiment assays colData SummarizedExperiment
#' @export
    
calcExprFreqs <- function(x, assay = "counts", th = 0) {
    # check validity of input arguments
    .check_sce(x, req_group = FALSE)
    .check_arg_assay(x, assay)
    stopifnot(is.numeric(th), length(th) == 1)
    
    # split cells by cluster-sample
    cs_by_ks <- .split_cells(x)

    # for each gene, compute fraction of cells 
    # w/ assay value above threshold in each sample
    y <- assays(x)[[assay]]
    fq <- lapply(cs_by_ks, vapply, function(i) {
        if (is.null(i)) return(rep(0, nrow(x)))
        Matrix::rowMeans(y[, i, drop = FALSE] > th)
    }, numeric(nrow(x)))

    # same for ea. group (if colData column "group_id" exists)
    if ("group_id" %in% colnames(colData(x))) {
        kids <- x$cluster_id
        nc_by_ks <- table(kids, x$sample_id)
        nc_by_kg <- table(kids, x$group_id)
        
        ei <- metadata(x)$experiment_info
        s_by_g <- split(ei$sample_id, ei$group_id) 
        gids <- names(s_by_g)
        kids <- set_names(levels(kids))
        
        fq <- lapply(kids, function(k) {
            nc_by_s <- fq[[k]] * nc_by_ks[k, ][col(fq[[k]])]
            fq_by_g <- vapply(gids, function(g) {
                nc_g <- nc_by_s[, s_by_g[[g]], drop = FALSE]
                rowSums(nc_g) / nc_by_kg[k, g]
            }, numeric(nrow(x)))
            cbind(fq[[k]], fq_by_g)
        })
    }
    # return SCE
    SingleCellExperiment(assays = fq)
}
