#' @rdname pbDS
#' @title Cluster-specific DE analysis using
#'
#' @description \code{run_edgeR} tests for cluster-specific 
#'   differential expression by aggregating single-cell 
#'   measurements and using \code{edgeR} for testing.
#'
#' @param x 
#'   a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param pb a named list of pseudo-bulk data for each cluster
#'   computed with \code{aggregateData}.
#' @param design 
#'   a design matrix with row and column names
#'   created with \code{\link[stats]{model.matrix}}.
#' @param contrast 
#'   a matrix of contrasts created with \code{\link[limma]{makeContrasts}}.
#' @param coef 
#'   passed to \code{\link[edgeR]{glmQLFTest}}.
#'   Ignored if \code{contrast} is not NULL.
#' @param method 
#'   a character string.
#' @param min_cells a numeric. 
#'   Specifies the minimum number of cells in a given cluster-sample 
#'   required to consider the sample for differential testing.
#' @param verbose 
#'   logical. Should information on progress be reported?
#'
#' @details \code{run_edgeR} tests for cluster-specific 
#'   differential expression by aggregating single-cell 
#'   measurements. Depending on the selected \code{method},
#'   differential testing is performed on pseudo-bulk data
#'   obtained via...
#'   \describe{
#'   \item{\code{raw_counts}}{
#'     summing, for every gene, raw counts for each cluster-sample.}
#'   \item{\code{normed_counts}}{
#'     summing, for every gene, normalized counts for each cluster-sample.}
#'   \item{\code{scaled_cpm}}{
#'     summing, for every gene, scaled CPM for each cluster-sample.
#'     Scaled CPM are obtained by multiplying pseudo-bulk raw counts
#'     by effective library sizes and dividing by 1M.}}
#'
#' @return a list containing 
#' \itemize{
#' \item a data.frame with differential testing results,
#' \item a \code{\link[edgeR]{DGEList}} object of length nb.-clusters, and
#' \item the \code{design} matrix, and \code{contrast} or \code{coef} used.
#' }
#'
#' @examples
#' # simulate 5 clusters, 20% of DE genes
#' data(sce)
#'     
#' # compute pseudo-bulk counts
#' pb <- aggregateData(sce)
#' 
#' # specify design & contrast matrix
#' ei <- metadata(x)$experiment_info
#' design <- model.matrix(~ 0 + ei$group_id)
#' dimnames(design) <- list(ei$sample_id, levels(ei$group_id))
#' contrast <- limma::makeContrasts("B-A", levels = design)
#' 
#' # test for cluster-specific DE 
#' res <- pbDS(sce, pb, design, contrast, method = "edgeR")
#'
#' names(res)
#' names(res[[1]])
#' lapply(res[[1]]$`stim-ctrl`, head)
#' 
#' # count nb. of DE genes by cluster
#' vapply(res[[1]]$`stim-ctrl`, function(x) 
#'   sum(x$p_adj < 0.05), numeric(1))
#' 
#' # get top 5 hits for ea. cluster w/ abs(logFC) > 1
#' library(dplyr)
#' lapply(res[[1]]$`stim-ctrl`, function(u)
#'   filter(u, abs(logFC) > 1) %>% 
#'     arrange(p_adj) %>% 
#'     slice(seq_len(5)))
#'
#' @author Helena L. Crowell \email{helena.crowell@uzh.ch} and Mark D. Robinson.
#'
#' @importFrom edgeR calcNormFactors DGEList 
#'   estimateDisp glmQLFit glmQLFTest topTags
#' @importFrom dplyr rename
#' @importFrom limma contrasts.fit eBayes lmFit topTable voom
#' @importFrom SummarizedExperiment colData
#' @importFrom tibble add_column
#' @export

pbDS <- function(x, pb, 
    design, contrast = NULL, coef = NULL, 
    method = c("edgeR", "limma-trend", "limma-voom"),
    min_cells = 10, verbose = TRUE) {
    
    # check validty of input arguments
    .check_sce(x, req_group = TRUE)
    .check_pb(x, pb)
    stopifnot(is.null(design) | is.matrix(design))
    stopifnot(is.null(contrast) | is.matrix(contrast))
    stopifnot(is.null(coef) | is.numeric(coef))
    stopifnot(is.numeric(min_cells), length(min_cells) == 1)
    stopifnot(is.logical(verbose), length(verbose) == 1)
    method <- match.arg(method)
    
    if (!is.null(contrast)) {
        ctype <- "contrast"
        cs <- colnames(contrast)
        names(cs) <- cs
    } else {
        ctype <- "coef"
        cs <- vapply(coef, function(i) 
            paste(colnames(design)[i], collapse = "-"),
            character(1))
        names(cs) <- names(coef) <- cs
    }
    
    # compute cluster-sample counts
    n_cells <- table(x$cluster_id, x$sample_id)
    kids <- levels(x$cluster_id)
    names(kids) <- kids
    
    # for ea. cluster, run DEA
    res <- lapply(kids, function (k) {
        if (verbose) cat(k, "..", sep = "")
        y <- assays(pb)[[k]]
        y <- y[, n_cells[k, ] >= min_cells]
        d <- design[colnames(y), ]
        if (any(colSums(d) < 2)) 
            return(NULL)
        if (method == "edgeR") {
            y <- suppressMessages(DGEList(y, remove.zeros = TRUE))
            y <- calcNormFactors(y)
            y <- estimateDisp(y, d)
            fit <- glmQLFit(y, d)
            tt <- lapply(cs, function(c) {
                qlf <- glmQLFTest(fit, coef[[c]], contrast[, c])
                tt <- topTags(qlf, n = Inf, sort.by = "none")
                res_df(k, tt, ctype, c) %>% 
                    rename(p_val = "PValue", p_adj.loc = "FDR")
            })
        } else {
            if (method == "limma-trend") {
                trend <- robust <- TRUE
                y <- switch(metadata(pb)$agg_pars$assay,
                    counts = { # raw counts > compute logCPM
                        y <- suppressMessages(DGEList(y, remove.zeros = TRUE))
                        y <- calcNormFactors(y)
                        cpm(y, log = TRUE, prior.count = 3)
                    },
                    logcounts = y, # log-normcounts > do nothing
                    log2(y + 1))   # CPM, scaledCPM, normcounts > take log
            } else if (method == "limma-voom") {
                trend <- robust <- FALSE
                y <- suppressMessages(DGEList(y, remove.zeros = TRUE))
                y <- calcNormFactors(y)
                y <- voom(y, d)
            }
            w <- n_cells[k, colnames(y)]
            fit <- lmFit(y, d, weights = w)
            tt <- lapply(cs, function(c) {
                cfit <- contrasts.fit(fit, contrast[, c], coef[[c]])
                efit <- eBayes(cfit, trend = trend, robust = robust)
                tt <- topTable(efit, number = Inf, sort.by = "none")  
                res_df(k, tt, ctype, c) %>% 
                    rename(p_val = "P.Value", p_adj.loc = "adj.P.Val")
            })
        }
        return(list(tt = tt, data = y))
    })
    # remove empty clusters
    skipped <- vapply(res, is.null, logical(1))
    if (any(skipped))
        message(paste("Cluster(s)", dQuote(kids[skipped]), "skipped due to an",
            "insufficient number of cells in at least 2 samples per group."))
    res <- res[!skipped]
    kids <- kids[names(res)]

    # re-organize by contrast & 
    # do global p-value adjustment
    tt <- lapply(res, "[[", "tt")
    tt <- lapply(cs, function(c) map(tt, c))
    tt <- .p_adj_global(tt)

    # return results
    data <- lapply(res, "[[", "data")
    list(table = tt, 
        data = data, 
        design = design, 
        contrast = contrast, 
        coef = coef)
}
