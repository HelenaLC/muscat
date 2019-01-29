#' @rdname prepData
#' @title Prepare SCE for DS analysis
#' 
#' @description ...
#' 
#' @param x a \code{[SingleCellExperiment]{SingleCellExperiment}}.
#' @param cluster_id a character string
#'   specifying the column in \code{colData(x)}
#'   that contains cluster assignments.
#' @param sample_id a character string
#'   specifying the column in \code{colData(x)}
#'   that contains unique sample identifies.
#' @param factors a character string
#'   specifying columns  in \code{colData(x)}
#'   that containfactors of interest.
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

prepData <- function(x, cluster_id, sample_id, factors) {
    
    # construct column data
    sample_ids <- factor(colData(x)[[sample_id]])
    col_data <- data.frame(
        row.names = colnames(x),
        sample_id = sample_ids,
        cluster_id = colData(x)[[cluster_id]],
        colData(x)[factors])
    
    # construct experimental design table
    m <- match(levels(sample_ids), sample_ids)
    ei <- data.frame(
        row.names = NULL,
        sample_id = sample_ids[m],
        colData(x)[m, ][factors],
        n_cells = as.numeric(table(sample_ids)))
    
    # construct SingleCellExperiment
    SingleCellExperiment(
        assays = assays(x),
        colData = col_data,
        rowData = rowData(x),
        metadata = list(experiment_info = ei),
        reducedDims = reducedDims(x))
}
