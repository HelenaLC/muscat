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
#' @importFrom SingleCellExperiment reducedDims SingleCellExperiment
#' @importFrom SummarizedExperiment assays colData rowData
#' @export

prepData <- function(x, cluster_id, sample_id, factors) {
    # construct experimental design table

  sample_ids <- factor(colData(x)[[sample_id]])
    col_data <- data.frame(
        row.names = colnames(x),
        sample_id = sample_ids,
        cluster_id = colData(x)[[cluster_id]],
        sapply(factors, function(i) colData(x)[[i]]))

    
    m <- match(levels(sample_ids), sample_ids)
    ei <- data.frame(
      row.names = NULL,
      sample_id = sample_ids[m],
      sapply(factors, function(i) colData(x)[[i]][m]),
      n_cells = as.numeric(table(sample_ids)))
    
    SingleCellExperiment(
        assays = assays(x),
        colData = col_data,
        rowData = rowData(x),
        metadata = list(experiment_info = ei),
        reducedDims = reducedDims(x))
}
    