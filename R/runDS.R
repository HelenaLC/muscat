#' @rdname runDS
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
#'   a matrix of contrasts created with \code{\link[edgeR]{makeContrasts}}.
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
#' data(kang)
#' sim <- simData(kang, n_genes = 100, n_cells = 200, 
#'     p_dd = c(0.8, 0, 0.2, 0, 0, 0), fc = 4)
#'     
#' # compute pseudo-bulk counts
#' pb <- aggregateData(sim, assay = "counts", fun = "sum")
#' 
#' # specify design & contrast matrix
#' ei <- metadata(sim)$experiment_info
#' design <- model.matrix(~ 0 + ei$group)
#' dimnames(design) <- list(ei$sample_id, levels(ei$group))
#' contrast <- limma::makeContrasts("B-A", levels = design)
#' 
#' # test for cluster-specific DE 
#' res <- runDS(sim, pb, design, contrast, method = "edgeR")
#'
#' names(res)
#' names(res[[1]])
#' lapply(res[[1]]$`B-A`, head)
#' 
#' # count nb. of DE genes by cluster
#' n_de <- sapply(res[[1]]$`B-A`, function(x) sum(x$FDR < 0.05))
#' 
#' # get top_n hits for ea. cluster w/ abs(logFC) > 1
#' top_n <- 5
#' lapply(res[[1]]$`B-A`, function(x) {
#'   x <- x[abs(x$logFC) > 1, ]
#'   x[order(x$FDR)[seq_len(top_n)], ]
#' })
#'
#' @author Helena L. Crowell \email{helena.crowell@uzh.ch} and Mark D. Robinson.
#'
#' @importFrom edgeR calcNormFactors DGEList estimateDisp glmQLFit glmQLFTest topTags
#' @importFrom dplyr rename
#' @importFrom limma contrasts.fit eBayes lmFit topTable
#' @importFrom SummarizedExperiment colData
#'
#' @export

runDS <- function(x, pb, 
    design, contrast = NULL, coef = NULL, 
    method = c("edgeR", "limma"),
    min_cells = 10, verbose = TRUE) {
    
    # check validty of input arguments
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot(all(c("sample_id", "cluster_id") %in% colnames(colData(x))))
    stopifnot(all(names(pb) %in% levels(colData(x)$cluster_id)))
    stopifnot(is.matrix(design))
    stopifnot(!is.null(contrast) | !is.null(coef))
    stopifnot(is.null(contrast) | is.matrix(contrast))
    stopifnot(is.null(coef) | is.numeric(coef))
    method <- match.arg(method)
    
    # compute cluster-sample counts
    n_cells <- table(colData(x)$cluster_id, colData(x)$sample_id)

    if (!is.null(contrast)) {
        ctype <- "contrast"
        cs <- colnames(contrast)
        names(cs) <- cs
    } else {
        ctype <- "coef"
        cs <- vapply(coef, function(i) 
            paste(colnames(design)[i], collapse = "--"),
            character(1))
        names(cs) <- names(coef) <- cs
    }
    cluster_ids <- levels(colData(x)$cluster_id)
    names(cluster_ids) <- cluster_ids
    
    # wrapper to create output tables
    res_df <- function(k, tt, ctype, c) {
        df <- data.frame(gene = rownames(tt), cluster_id = k,
            tt, row.names = NULL, stringsAsFactors = FALSE)
        df[[ctype]] <- c
        return(df)
    }
    # for ea. cluster, run DEA
    res <- lapply(cluster_ids, function (k) {
        if (verbose) cat(k, "..", sep = "")
        y <- pb[[k]][, n_cells[k, ] >= min_cells]
        d <- design[colnames(y), ]
        if (any(colSums(d) < 2)) return(NULL)
        tt <- switch(method,
            limma = {
                w <- n_cells[k, colnames(y)]
                fit <- lmFit(y, d, weight = w)
                lapply(cs, function(c) {
                    cfit <- contrasts.fit(fit, contrast[, c], coef[[c]])
                    efit <- eBayes(cfit, trend = TRUE)
                    tt <- topTable(efit, number = Inf, sort.by = "none")
                    res_df(k, tt, ctype, c) %>% 
                        rename(p_val = "P.Value", p_adj = "adj.P.Val")
                })
            },
            edgeR = {
                y <- suppressMessages(DGEList(y, remove.zeros = TRUE))
                y <- calcNormFactors(y)
                y <- estimateDisp(y, d)
                fit <- glmQLFit(y, d)
                lapply(cs, function(c) {
                    qlf <- glmQLFTest(fit, coef = coef[[c]], contrast = contrast[, c])
                    tt <- topTags(qlf, n = Inf, sort.by = "none")
                    res_df(k, tt, ctype, c) %>% 
                        rename(p_val = "PValue", p_adj = "FDR")
                })
            })
        return(list(tt = tt, data = y))
    })
    # remove empty clusters
    res <- res[!vapply(res, is.null, logical(1))]
    data <- lapply(res, "[[", "data")
    tt <- lapply(res, "[[", "tt")
    
    # re-organize results by comparison
    cluster_ids <- cluster_ids[names(res)]
    tt <- lapply(cs, function(c) 
        lapply(cluster_ids, function(k) tt[[k]][[c]]))
    
    # return results
    list(table = tt, 
        data = data, 
        design = design, 
        contrast = contrast, 
        coef = coef)
}

    
