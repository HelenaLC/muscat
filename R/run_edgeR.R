#' @rdname run_edgeR
#' @title Cluster-specific DE analysis using \code{edgeR} 
#'
#' @description \code{run_edgeR} tests for cluster-specific 
#'   differential expression by aggregating single-cell 
#'   measurements and using \code{edgeR} for testing.
#'
#' @param x 
#'   a \code{[SingleCellExperiment]{SingleCellExperiment}}.
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
#' @author Helena Lucia Crowell \email{helena@crowells.eu}
#'
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @importFrom edgeR DGEList estimateDisp glmQLFit glmQLFTest topTags
#' @importFrom dplyr group_by_ select summarise_at ungroup %>%
#' @importFrom reshape2 dcast
#' @importFrom scater calculateCPM normalize
#' @importFrom tidyr complete
#'
#' @export

run_edgeR <- function(x, pb, design, contrast = NULL, coef = NULL, min_cells = 10, verbose = TRUE) {
    
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
    # w/ non-zero counts in each cluster & group
    p_cells <- df %>% 
        group_by_(~cluster_id, ~group) %>% 
        summarise_at(rownames(x), function(g) mean(g > 0))
    p_cells <- split(p_cells, p_cells$cluster_id)
    p_cells <- lapply(p_cells, function(x)
        x %>% ungroup() %>% select(-"cluster_id") %>% 
            data.frame(row.names = 1) %>% t())
    
    # for ea. cluster, run DEA w/ edgeR
    cluster_ids <- levels(colData(x)$cluster_id)
    res <- setNames(lapply(cluster_ids, function(k) {
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
            df <- apply(contrast, 2, function(c) {
                qlf <- glmQLFTest(fit, contrast = c)
                tt <- topTags(qlf, n = Inf, p = Inf)$table
                df <- data.frame(
                    gene = rownames(tt), 
                    cluster_id = k, tt,
                    row.names = NULL,
                    stringsAsFactors = FALSE)
                cbind(df, p_cells[[k]][df$gene, ])
            })
        } else {
            qlf <- glmQLFTest(fit, coef = coef)
            tt <- topTags(qlf, n = Inf, p = Inf)$table
            df <- data.frame(
                gene = rownames(tt), 
                cluster_id = k, tt,
                row.names = NULL,
                stringsAsFactors = FALSE)
            df <- cbind(df, p_cells[[k]][df$gene, ])
        }
        return(list(tt = df, dgel = y))
    }), cluster_ids)
    # remove skipped clusters
    res <- res[!sapply(res, is.null)]
    
    dgel <- lapply(res, "[[", "dgel")
    res <- lapply(res, "[[", "tt")
    if (!is.null(contrast)) {
        res <- lapply(colnames(contrast), function(c)
            data.frame(
                do.call(rbind, lapply(res, "[[", c)),
                contrast = c, row.names = NULL))
        res <- do.call(rbind, res)
    } else {
        res <- do.call(rbind, res)
    }
    
    # return results
    list(res, 
        data = dgel,
        design = design,
        contrast = contrast,
        coef = coef)
}
