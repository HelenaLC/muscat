#' @rdname run_edgeR
#' @title Cluster-specific DE analysis using \code{edgeR} 
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
#'   a character string (see details).
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
#' sim <- simData(kang, n_genes = 10, n_cells = 100, 
#'     p_dd = c(0.8, 0, 0.2, 0, 0, 0))
#'     
#' # compute pseudo-bulk counts
#' pb <- aggregateData(sim, data = "counts", fun = "sum")
#' 
#' # specify design & contrast matrix
#' ei <- metadata(sim)$experiment_info
#' design <- model.matrix(~ 0 + ei$group)
#' dimnames(design) <- list(ei$sample_id, levels(ei$group))
#' contrast <- makeContrasts("B-A", levels = design)
#' 
#' # test for cluster-specific DE 
#' res <- run_edgeR(sim, pb, design, contrast)
#' head(res[[1]])
#' 
#' # count nb. of DE genes by cluster
#' sapply(split(res[[1]], res[[1]]$cluster_id), 
#'     function(x) table(x$FDR < 0.05))
#'
#' @author Helena L. Crowell \email{helena.crowell@uzh.ch} and Mark D. Robinson.
#'
#' @import SingleCellExperiment
#' @importFrom edgeR calcNormFactors DGEList estimateDisp glmQLFit glmQLFTest topTags
#' @importFrom dplyr group_by_ select summarise_at ungroup %>%
#' @importFrom reshape2 dcast
#' @importFrom scater calculateCPM normalize
#' @importFrom tidyr complete
#'
#' @export

run_edgeR <- function(x, pb, 
    design, contrast = NULL, coef = NULL, 
    min_cells = 10, verbose = TRUE) {
    
    # check validty of input arguments
    stopifnot(class(x) == "SingleCellExperiment")
    stopifnot(all(names(pb) %in% levels(colData(x)$cluster_id)))
    stopifnot(is.matrix(design))
    stopifnot(!is.null(contrast) | !is.null(coef))
    stopifnot(is.null(contrast) | is.matrix(contrast))
    stopifnot(is.null(coef) | is.numeric(coef))
    
    # compute cluster-sample counts
    df <- data.frame(
        t(assays(x)$counts), colData(x),
        check.names = FALSE)
    n_cells <- df %>% 
        count_(c("cluster_id", "sample_id")) %>% 
        acast(cluster_id ~ sample_id, value.var = "n", fill = 0)
    
    # for each gene, compute percentage of cells 
    # w/ non-zero counts in each cluster-sample & group
    cluster_ids <- factor(colData(x)$cluster_id)
    sample_ids <- factor(colData(x)$sample_id)
    
    idx <- split(seq_len(ncol(x)), list(cluster_ids, sample_ids))
    p_cells <- sapply(idx, function(i) 
        rowMeans(assays(x)$counts[, i, drop = FALSE] > 0))
    df <- data.frame(
        index = seq_len(ncol(p_cells)),
        cluster_id = rep(levels(cluster_ids), nlevels(sample_ids)),
        sample_id = rep(levels(sample_ids), each = nlevels(cluster_ids)))
    
    dfs <- split(df, df$cluster_id)
    p_cells <- lapply(dfs, function(df) {
        x <- p_cells[, df$index, drop = FALSE]
        colnames(x) <- df$sample_id
        return(x)
    })

    # for ea. cluster, run DEA w/ edgeR
    res <- lapply(levels(cluster_ids), function(k) {
        if (verbose) message(k, "..", appendLF = FALSE)
        y <- pb[[k]]
        # remove samples w/ less than min_cells
        y <- y[, n_cells[k, ] >= min_cells]
        y <- suppressMessages(DGEList(y, remove.zeros = TRUE))
        y <- calcNormFactors(y)
        d <- design[colnames(y), ]
        # skip cluster if there are less than 2 samples in any group
        if (any(colSums(d) < 2)) return(NULL)
        y <- estimateDisp(y, d)
        fit <- glmQLFit(y, d)
        if (!is.null(contrast)) {
            df <- do.call(rbind, lapply(colnames(contrast), function(c) {
                qlf <- glmQLFTest(fit, contrast = contrast[, c])
                tt <- topTags(qlf, n = Inf, p.value = Inf)$table
                gs <- rownames(tt)
                data.frame(
                    gene = gs, 
                    cluster_id = k, tt,
                    p_cells[[k]][gs, ], 
                    contrast = c,
                    row.names = NULL, 
                    stringsAsFactors = FALSE)
            }))
        } else {
            qlf <- glmQLFTest(fit, coef = coef)
            tt <- topTags(qlf, n = Inf, p.value = Inf)$table
            gs <- rownames(tt)
            df <- data.frame(
                gene = gs, 
                cluster_id = k, tt,
                p_cells[[k]][gs, ],
                row.names = NULL,
                stringsAsFactors = FALSE)
        }
        return(list(tt = df, dgel = y))
    })
    # remove skipped clusters
    res <- res[!sapply(res, is.null)]
    
    dgel <- lapply(res, "[[", "dgel")
    res <- lapply(res, "[[", "tt")
    res <- do.call(rbind, res)
    
    # return results
    list(res, 
        data = dgel,
        design = design,
        contrast = contrast,
        coef = coef)
}
