#' @rdname aggregateData
#' @title Aggregation of single-cell to pseudo-bulk data
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
#' @param method a character string
#'   specifying the method to use for aggregation.
#' 
#' @import SingleCellExperiment
#' @importFrom DelayedArray t
#' @importFrom dplyr group_by_ summarise_at ungroup select
#' @importFrom scater calculateCPM
#' @importFrom tidyr complete
#' @export

aggregateData <- function(x, method = c("rawCounts", "normedCounts", "scaledCPM")) {
    
    # validity checks for input arguments
    stopifnot(class(x) == "SingleCellExperiment")
    method <- match.arg(method)
    
    # get input data according to 'method'
    data <- switch(method,
        rawCounts = assays(x)$counts,
        normedCounts = {
            if (!"normcounts" %in% assayNames(x))
                x <- normalize(x, return_log = FALSE)
            assays(x)$normcounts
        },
        scaledCPM = calculateCPM(assays(x)$counts))
    
    df <- data.frame(t(data),
        cluster_id = colData(x)$cluster_id,
        sample_id = colData(x)$sample_id,
        check.names = FALSE)
    
    # compute pseudo-bulks
    pb <- df %>% 
        group_by_(~cluster_id, ~sample_id) %>% 
        complete(sample_id) %>% 
        summarise_at(rownames(x), sum)
    
    if (method == "scaledCPM") {
        df <- data.frame(
            t(assays(x)$counts), 
            cluster_id = colData(x)$cluster_id,
            sample_id = colData(x)$sample_id,
            check.names = FALSE)
        lib_sizes <- rowSums(df %>%
            group_by_(~cluster_id, ~sample_id) %>% 
            complete(sample_id) %>% 
            summarise_at(rownames(x), sum) %>%
            ungroup() %>% select(rownames(x)))
        pb <- pb * lib_sizes / 1e6
    }
    
    # split by cluster
    pb <- split(pb, pb$cluster_id)
    pb <- lapply(pb, function(x) {
        x <- x %>% ungroup() %>% select(-"cluster_id")
        t(data.frame(x, row.names = 1))
    })
    return(pb)
}

