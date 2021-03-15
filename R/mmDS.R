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
#'   Either \code{"dream2"} (default, lme4 with voom-weights),
#'   \code{"dream"} (previous implementation of the dream method),
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
#' @param BPPARAM a \code{\link[BiocParallel]{BiocParallelParam}}
#'   object specifying how differential testing should be parallelized.
#' @param verbose logical specifying whether messages
#'   on progress and a progress bar should be displayed.
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
#' gs <- sample(nrow(sce), 100)
#' sce <- sce[gs, ]
#'
#' res <- mmDS(sce, method = "dream", verbose = FALSE)
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
#' @importFrom sctransform vst
#' @importFrom SingleCellExperiment counts counts<-
#'   colData sizeFactors sizeFactors<-
#' @importFrom stats p.adjust
#' @export

mmDS <- function(x, coef = NULL, covs = NULL,
    method = c("dream2", "dream", "vst", "poisson", "nbinom", "hybrid"),
    n_cells = 10, n_samples = 2, min_count = 1, min_cells = 20,
    verbose = TRUE, BPPARAM = SerialParam(progressbar = verbose), 
    vst = c("sctransform", "DESeq2"),
    ddf = c("Satterthwaite", "Kenward-Roger", "lme4"),
    dup_corr = FALSE, trended = FALSE, bayesian = FALSE, 
    blind = TRUE, REML = TRUE, moderate = FALSE) {
    
    # check validity of input arguments
    .check_sce(x, req_group = TRUE)
    .check_arg_assay(x, "counts")
    .check_args_mmDS(as.list(environment()))
    stopifnot(is(BPPARAM, "BiocParallelParam"))
    
    args <- as.list(environment())
    args$method <- match.arg(method)
    args$vst <- match.arg(vst)
    args$ddf <- match.arg(ddf)
    
    # counts cells per cluster-sample
    n_cells_by_ks <- table(x$cluster_id, x$sample_id)
    
    # filter clusters w/ >= n_cells in >= n_samples
    if (!is.null(metadata(x)$experiment_info$group_id)) {
        ei <- metadata(x)$experiment_info
        m <- match(levels(x$sample_id), ei$sample_id)
        gids <- ei$group_id[m]
    } else {
        gids <- x$group_id
    }
    ks_keep <- apply(n_cells_by_ks > n_cells, 1,
        function(u) all(tabulate(gids[u]) >= n_samples))
    if (sum(ks_keep) == 0)
        stop(paste("No cluster has at least", n_samples,
            "samples with at least", n_cells, "cells."))
    
    kids <- levels(x$cluster_id)
    if (verbose && sum(ks_keep) < length(kids))
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
        z <- do.call(fun, args)
        z <- cbind(stringsAsFactors = FALSE,
            gene = rownames(y), cluster_id = k, z)
        rownames(z) <- NULL; z
    })
    if (verbose) pb$terminate()
    
    # assemble results from all cluster
    res <- bind_rows(res)
    # global p-value adjustment
    p_adj.glb <- p.adjust(res$p_val)
    i <- which(colnames(res) == "p_adj.loc")
    res[["p_adj.glb"]] <- p_adj.glb
    res <- res[, c(seq_len(i), ncol(res), seq_len(ncol(res)-1)[-seq_len(i)])]
    # re-split by cluster
    split(res, res$cluster_id)
}
