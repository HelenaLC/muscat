#' @rdname aggregateData
#' @title Aggregation of single-cell to pseudo-bulk data
#' 
#' @description ...
#' 
#' @param x a \code{[SingleCellExperiment]{SingleCellExperiment}}.
#' @param data a character string specifying the assay slot to use as input data.
#' @param fun a character string specifying the function to use as summary statistic.
#' @param scale logical.
#' 
#' @import SingleCellExperiment
#' @importFrom DelayedArray t
#' @importFrom dplyr group_by_ summarise_at ungroup select
#' @importFrom scater calculateCPM
#' @importFrom tidyr complete
#' @export

aggregateData <- function(x, data, fun, scale = FALSE) {
    
    # validity checks for input arguments
    stopifnot(class(x) == "SingleCellExperiment")
    stopifnot(all(c("cluster_id", "sample_id") %in% colnames(colData(x))))
    stopifnot(is.character(data), length(data) == 1, data %in% assayNames(x))
    stopifnot(is.character(fun), length(fun) == 1, exists(fun, mode = "function"))
    stopifnot(is.logical(scale), length(scale) == 1)

    df <- data.frame(
        t(assays(x)[[data]]),
        cluster_id = colData(x)$cluster_id,
        sample_id = colData(x)$sample_id,
        check.names = FALSE)
    
    # compute pseudo-bulks
    zeros <- as.list(numeric(nrow(sce)))
    zeros <- setNames(zeros, rownames(x))
    pb <- df %>% 
        group_by_(~cluster_id, ~sample_id) %>% 
        summarise_at(rownames(x), get(fun)) %>% 
        complete(sample_id, fill = zeros)
    
    # scale
    if (scale) {
        if (data == "counts" & fun == "sum") {
            lib_sizes <- rowSums(pb %>% select(rownames(x)))
        } else {
            lib_sizes <- df %>% 
                group_by_(~cluster_id, ~sample_id) %>%
                summarise_at(rownames(x), sum) %>%
                complete(sample_id, fill = zeros) %>%
                ungroup() %>% select(rownames(x))
            lib_sizes <- rowSums(lib_sizes)
        }
        pb <- pb * lib_sizes / 1e6
    }
    
    # split by cluster
    pb <- split(pb, pb$cluster_id)
    lapply(pb, function(x) {
        x <- x %>% ungroup() %>% select(-"cluster_id")
        t(data.frame(x, row.names = 1))
    })
}