#' @rdname mmDS
#' @title DS analysis using mixed-models (MM)
#'
#' @description Performs cluster-wise DE analysis by fitting cell-level models.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param coef character specifying the coefficient to test. 
#'   If NULL (default), will test the last level of \code{"group_id"}.
#' @param covs character vector of \code{colData(x)} 
#'   column names to use as covariates.
#' @param method a character string. 
#'   Either \code{"dream"} (default, lme4 with voom-weights), 
#'   \code{"vst"} (variance-stabilizing transformation), 
#'   or a function accepting the same arguments.
#' @param n_cells number of cells per cluster-sample 
#'   required to consider a sample for testing.
#' @param n_samples number of samples per group 
#'   required to consider a cluster for testing.
#' @param min_count numeric. For a gene to be tested in a given cluster, 
#'   at least \code{min_cells} must have a count >= \code{min_count}.
#' @param min_cells number (or fraction, if < 1) of cells with a count >
#'   \code{min_count} required for a gene to be tested in a given cluster.
#' @param n_threads number of threads to use.
#' @param verbose logical specifying whether messages 
#'   on progress and a progress bar should be displayed.
#' @param ... additional parameters passed to the method function.
#'
#' @return a data.frame
#' 
#' @author Pierre-Luc Germain and Helena L. Crowell.
#' 
#' @importFrom dplyr %>% mutate bind_rows
#' @importFrom progress progress_bar
#' @importFrom purrr map_depth
#' @importFrom tibble add_column
#' @importFrom SingleCellExperiment colData counts
#' @importFrom stats p.adjust
#' @export
mmDS <- function(x, 
    coef = NULL, covs = NULL, method = c("dream", "vst"),
    n_cells = 10, n_samples = 2, min_count = 1, min_cells = 20,
    n_threads = 8, verbose = TRUE, ...) {
    
    .check_sce(x, req_group = TRUE)
    .check_arg_assay(x, "counts")

    if (!is.null(covs) && !all(covs %in% names(colData(x))))
        stop(paste("Some of the specified covariates couldn't be found:",
            paste(setdiff(covs, names(colData(x))), collapse=", ")))
    
    # get method function
    method <- match.arg(method)
    fun <- ifelse(is.function(method), method, get(paste0(".mm_", method)))
    
    # counts cells per cluster-sample
    n_cells_by_ks <- table(x$cluster_id, x$sample_id)
    
    # filter clusters w/ >= n_cells in >= n_samples
    ei <- metadata(x)$experiment_info
    m <- match(levels(x$sample_id), ei$sample_id)
    gids <- ei$group_id[m]
    ks_keep <- apply(n_cells_by_ks > n_cells, 1, 
        function(u) all(table(gids[u]) >= n_samples))
    if (sum(ks_keep) == 0) 
        stop(paste("No cluster has at least", n_samples, 
            "samples with at least", n_cells, "cells."))
    
    kids <- levels(x$cluster_id)
    if (sum(ks_keep) < length(kids))
        message(paste("Skipping cluster(s)", 
            paste(dQuote(kids[!ks_keep]), collapse = ", "), 
            "\ndue to an insufficient number of samples",
            "with a sufficient number of cells."))
    kids <- kids[ks_keep]
    names(kids) <- kids
    
    # split cells by cluster
    cells_by_k <- split(colnames(x), x$cluster_id)
    
    if (verbose) pb <- progress_bar$new(total = length(kids))
    res <- lapply(kids, function(k) {
        y <- x[, cells_by_k[[k]]]
        if (min_count < 1) 
            min_count <- floor(min_count * min(n_cells_by_ks[k, ]))
        y <- y[rowSums(counts(y) >= min_count) > min_cells, ]
        if (verbose) 
            message(sprintf(
                "Testing %s genes across %s cells in cluster %s...", 
                nrow(y), ncol(y), dQuote(k)))
        
        # call to .mm_dream/.mm_vst
        fun(y, coef, covs, n_threads, verbose, ...) %>% 
            add_column(.before = 1, gene = rownames(y), cluster_id = k) %>% 
            set_rownames(NULL)
    }) 
    if (verbose) pb$terminate()
    
    # global p-value adjustment
    bind_rows(res) %>% 
        add_column(.after = "p_adj.loc", p_adj.glb = p.adjust(.$p_val)) %>% 
        split(.$cluster_id)
}
