#' pbFlatten
#' Flatten pseudobulk SCE
#'
#' Flattens a pseudobulk \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' as returned by \code{\link{aggregateData}} such that all cell subpopulations 
#' are represented as a single assay.
#'
#' @param pb a pseudobulk \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   as returned by \code{\link{aggregateData}}, with different subpopulations as assays.
#' @param normalize logical specifying whether to compute a \code{logcpm} assay.
#'
#' @return a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' 
#' @examples 
#' data(example_sce)
#' library(SingleCellExperiment)
#' pb_stack <- aggregateData(example_sce)
#' (pb_flat <- pbFlatten(pb_stack))
#' ncol(pb_flat) == ncol(pb_stack)*length(assays(pb_stack))
#' 
#' @importFrom methods is
#' @importFrom edgeR cpm calcNormFactors DGEList
#' @importFrom S4Vectors DataFrame metadata as.list
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment assay assay<- assays colData rowData
#' @export

pbFlatten <- function(pb, normalize = TRUE){
    # check validity of input arguments
    stopifnot(is.logical(normalize), length(normalize) == 1,
        is(pb, "SingleCellExperiment"), length(assays(pb)) > 1)
    # concatenate assay data
    as <- assays(pb)
    a <- do.call(cbind, as.list(as))
    sids <- rep(colnames(pb), length(as))
    kids <- rep(names(as), each = ncol(pb))
    colnames(a) <- paste(sids, kids, sep = ".")
    pb$sample_id <- colnames(pb)
    # construct cell metadata
    cd <- lapply(seq_along(as), function(u) colData(pb))
    cd <- do.call(rbind, cd)
    rownames(cd) <- colnames(a)
    cd$cluster_id <- kids
    # construct single-assay SCE
    sce <- SingleCellExperiment(
        assays = list(counts = a),
        colData = cd, 
        rowData = rowData(pb),
        metadata = metadata(pb))
    # (optionally) add number of cells per cluster-sample
    if (!is.null(n_cells <- .n_cells(pb))) {
        n_cells <- tryCatch(mapply(
            function(k, s) n_cells[k, s],
            k = as.character(kids), s = as.character(sids)),
            error = function(e) {warning(e); NULL})
        if (!is.null(n_cells)) 
            sce$n_cells <- as.numeric(n_cells)
    }
    # (optionally) do logCPM normalization
    if (normalize){
	# remove empty columns (samples that lack a cluster)
	sce <- sce[,colSums(a!=0)>0]
        assay(sce, "logcpm") <-
          log1p(cpm(calcNormFactors(DGEList(assay(sce)))))
    }
    return(sce)
}
