#' @rdname prepData
#' @title Prepare SCE for DS analysis
#' 
#' @description ...
#' 
#' @param x a \code{[SingleCellExperiment]{SingleCellExperiment}}.
#' @param cluster_id,sample_id,group_id character strings specifying
#'   the \code{colData(x)} columns containing cluster assignments,
#'   unique sample identifiers, and group IDs (e.g., treatment).
#' 
#' @examples
#' # generate random counts
#' n_genes <- 50
#' n_cells <- 100
#' counts <- matrix(
#'     sample(n_genes * n_cells),
#'     nrow = n_genes, ncol = n_cells)
#'     
#' # generate some cell metadata
#' sids <- sample(paste0("sample", 1:5), n_cells, TRUE)  # sample IDs
#' cids <- sample(paste0("cluster", 1:5), n_cells, TRUE) # cluster IDs
#' condition <- sample(c("ctrl", "stim"), n_cells, TRUE) # conditions
#' batch <- sample(1:3, n_cells, TRUE)                   # batch nb.
#' 
#' # construct SCE
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(
#'     assays = list(counts = counts),
#'     colData = data.frame(id = sids, cluster = cids, batch, condition))
#'     
#' # prep. for workflow
#' sce <- prepData(sce, 
#'     cluster_id = "cluster", 
#'     sample_id = "id", 
#'     factors = c("batch", "condition"))
#'     
#' head(colData(sce))
#' metadata(sce)$experiment_info
#' 
#' @return a \code{[SingleCellExperiment]{SingleCellExperiment}}.
#' 
#' @importFrom SingleCellExperiment reducedDims SingleCellExperiment
#' @importFrom SummarizedExperiment assays colData rowData
#' @export

prepData <- function(x, cluster_id, sample_id, group_id) {
  
  # construct colData
  ids <- as.list(match.call()[-(1:2)])
  cd <- data.frame(
    row.names = colnames(x),
    colData(x)[unlist(ids)])
  colnames(cd) <- names(ids)
  
  # construct experimental design table
  m <- match(levels(cd$sample_id), cd$sample_id)
  ei <- data.frame(
    row.names = NULL,
    sample_id = levels(cd$sample_id,)
    group_id = cd$group_id[m],
    n_cells = as.numeric(table(cd$sample_id)))
  md <- list(experiment_info = ei)
  
  # construct SingleCellExperiment
  SingleCellExperiment(
    assays = assays(x),
    colData = cd,
    rowData = rowData(x),
    reducedDims = reducedDims(x),
    metadata = md)
}
