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
#'   \code{"poisson"} (poisson GLM-MM), \code{"nbinom"} (negative binomial GLM-MM), 
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
#' # downsample to 100 genes & 50 cells per sample
#' cells_by_sample <- split(colnames(sce), sce$sample_id)
#' genes_keep <- sample(nrow(sce), 100)
#' cells_keep <- sapply(cells_by_sample, sample, 50)
#' sce <- sce[genes_keep, cells_keep]
#' 
#' res <- mmDS(sce, method = "dream", verbose = FALSE)
#' head(res$`B cells`)
#' 
#' @author Pierre-Luc Germain and Helena L. Crowell.
#' 
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateDispersions
#'   sizeFactors varianceStabilizingTransformation
#' @importFrom dplyr %>% mutate bind_rows
#' @importFrom matrixStats rowMins
#' @importFrom progress progress_bar
#' @importFrom purrr map_depth
#' @importFrom tibble add_column
#' @importFrom scran computeSumFactors
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

    if (!is.null(covs) && !all(covs %in% names(colData(x))))
        stop(paste("Some of the specified covariates couldn't be found:",
            paste(setdiff(covs, names(colData(x))), collapse=", ")))
    
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

    if (min_count < 1) {
        min_count <- floor(min_count * rowMins(n_cells_by_ks))
    } else {
        min_count <- rep(min_count, length(kids))
    }
    names(min_count) <- kids
    
    # variance-stabilizing transformation
    if (args$method == "vst") {
        vst_call <- switch(args$vst,
            sctransform = {
                # assure correct vst() function is used
                fun <- getFromNamespace("vst", "sctransform")
                expression(fun(assay(x), show_progress = verbose)$y)
            },
            DESeq2 = {
                if (is.null(sizeFactors(x)))
                    x <- computeSumFactors(x)
                formula <- as.formula(paste("~", 
                    paste(c(covs, "sample_id"), collapse = "+")))
                y <- suppressMessages(DESeqDataSetFromMatrix(
                    as.matrix(counts(x)), colData(x), formula))
                sizeFactors(y) <- sizeFactors(x)
                if (!blind) y <- estimateDispersions(y)
                expression(assay(varianceStabilizingTransformation(y, blind)))
            })
        if (!verbose) vst_call <- parse(text = 
                sprintf("suppressMessages(%s)", paste(vst_call)))
        counts(x) <- eval(vst_call)
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
            message(sprintf(
                "Testing %s genes across %s cells in cluster %s...", 
                nrow(y), ncol(y), dQuote(k)))
        
        # call to .mm_dream/.mm_vst
        args$x <- y
        do.call(fun, args) %>% 
            add_column(.before = 1, gene = rownames(y), cluster_id = k) %>% 
            set_rownames(NULL)
    }) 
    if (verbose) pb$terminate()
    
    # global p-value adjustment
    bind_rows(res) %>% 
        add_column(.after = "p_adj.loc", p_adj.glb = p.adjust(.$p_val)) %>% 
        split(.$cluster_id)
}
