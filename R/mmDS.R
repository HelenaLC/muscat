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
#'   \code{"poisson"} (poisson GLM-MM), 
#'   \code{"nbinom"} (negative binomial GLM-MM), 
#'   \code{"hybrid"} (combination of pseudobulk and poisson methods)
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
#' @examples 
#' data(sce)
#' # subset "B cells" cluster
#' sce <- sce[, sce$cluster_id == "B cells"]
#' sce$cluster_id <- droplevels(sce$cluster_id)
#' 
#' # downsample to 100 genes
#' cs_by_s <- split(colnames(sce), sce$sample_id)
#' gs <- sample(nrow(sce), 100)
#' sce <- sce[gs, ]
#' 
#' res <- mmDS(sce, method = "dream", 
#'     n_threads = 1, verbose = FALSE)
#' head(res$`B cells`)
#' 
#' @author Pierre-Luc Germain & Helena L Crowell
#' 
#' @references 
#' Crowell, HL, Soneson, C, Germain, P-L, Calini, D, 
#' Collin, L, Raposo, C, Malhotra, D & Robinson, MD: 
#' On the discovery of population-specific state transitions from 
#' multi-sample multi-condition single-cell RNA sequencing data. 
#' \emph{bioRxiv} \strong{713412} (2018). 
#' doi: \url{https://doi.org/10.1101/713412}
#' 
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateDispersions
#'   sizeFactors varianceStabilizingTransformation
#' @importFrom dplyr %>% mutate bind_rows
#' @importFrom matrixStats rowMins
#' @importFrom progress progress_bar
#' @importFrom purrr map_depth
#' @importFrom tibble add_column
#' @importFrom sctransform vst
#' @importFrom SingleCellExperiment counts counts<- 
#'   colData sizeFactors sizeFactors<-
#' @importFrom stats p.adjust
#' @export

mmDS <- function(x, coef = NULL, covs = NULL, 
    method = c("dream", "vst", "poisson", "nbinom", "hybrid"),
    n_cells = 10, n_samples = 2, min_count = 1, min_cells = 20,
    n_threads = 8, verbose = TRUE, 
    dup_corr = FALSE, trended = FALSE,
    vst = c("sctransform", "DESeq2"), 
    bayesian = FALSE, blind = TRUE, REML = TRUE,
    ddf = c("Satterthwaite", "Kenward-Roger", "lme4")) {
    
    .check_sce(x, req_group = TRUE)
    .check_arg_assay(x, "counts")
    
    args <- as.list(environment())
    args$method <- match.arg(method)
    args$vst <- match.arg(vst)
    args$ddf <- match.arg(ddf)

    if (!is.null(covs) && !all(covs %in% names(colData(x)))) {
        txt <- paste(dQuote(setdiff(covs, names(colData(x)))), collapse = ", ")
        stop("Some of the specified covariates couldn't be found: ", txt)
            
    }
    
    # counts cells per cluster-sample
    n_cells_by_ks <- table(x$cluster_id, x$sample_id)
    
    # filter clusters w/ >= n_cells in >= n_samples
    ei <- metadata(x)$experiment_info
    m <- match(levels(x$sample_id), ei$sample_id)
    gids <- ei$group_id[m]
    ks_keep <- apply(n_cells_by_ks > n_cells, 1, 
        function(u) all(tabulate(gids[u]) >= n_samples))
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

    if (min_count < 1) {
        min_count <- floor(min_count * rowMins(n_cells_by_ks))
    } else {
        min_count <- rep(min_count, length(kids))
    }
    names(min_count) <- kids
    
    # variance-stabilizing transformation
    if (args$method == "vst") {
        vst_call <- switch(args$vst,
            sctransform = expression(.vst_sctransform(x, verbose)),
            DESeq2 = expression(.vst_DESeq2(x, covs, blind)))
        if (verbose) {
            counts(x) <- eval(vst_call)
        } else {
            counts(x) <- suppressMessages(eval(vst_call))
        }
    }
    
    # get method function & construct correct call
    fun <- ifelse(is.function(args$method), 
        args$method, get(paste0(".mm_", args$method)))
    args_use <- names(formals(fun))
    args <- args[names(args) %in% args_use]
    
    if (verbose) pb <- progress_bar$new(total = length(kids))
    res <- lapply(kids, function(k) {
        y <- x[, cells_by_k[[k]]]
        y <- y[rowSums(counts(y) >= min_count[k]) > min_cells, ]
        if (verbose) 
            message("Testing ", nrow(y), " genes across ", 
                ncol(y), " cells in cluster ", dQuote(k), "...")
        
        # call to .mm_dream/.mm_vst
        args$x <- y
        do.call(fun, args) %>% 
            add_column(.before = 1, gene = rownames(y), cluster_id = k) %>% 
            set_rownames(NULL)
    }) 
    if (verbose) pb$terminate()
    
    # assemble results from all cluster
    res <- bind_rows(res)
    # global p-value adjustment
    p_adj.glb <- p.adjust(res$p_val)
    res <- add_column(res, p_adj.glb, .after = "p_adj.loc")
    # re-split by cluster
    split(res, res$cluster_id)
}
