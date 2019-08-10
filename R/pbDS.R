#' @rdname pbDS
#' @title pseudobulk DS analysis
#'
#' @description \code{pbDS} tests for DS after aggregating single-cell 
#'   measurements to pseudobulk data, by applying bulk RNA-seq DE methods, 
#'   such as \code{edgeR}, \code{DESeq2} and \code{limma}.
#'
#' @param pb a \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   containing pseudobulks as returned by \code{\link{aggregateData}}.
#' @param design For methods \code{"edegR"} and \code{"limma"}, a design matrix 
#'   with row & column names(!) created with \code{\link[stats]{model.matrix}}; 
#'   For \code{"DESeq2"}, a formula with variables in \code{colData(pb)}.
#'   Defaults to \code{~ group_id} or the corresponding \code{model.matrix}.
#' @param contrast a matrix of contrasts to test for
#'   created with \code{\link[limma]{makeContrasts}}.
#' @param coef passed to \code{\link[edgeR]{glmQLFTest}}.
#'   Ignored if \code{contrast} is not NULL.
#' @param method a character string.
#' @param min_cells a numeric. Specifies the minimum number of cells in a given 
#'   cluster-sample required to consider the sample for differential testing.
#' @param verbose logical. Should information on progress be reported?
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
#' # compute pseudobulk sum-counts & run DS analysis
#' pb <- aggregateData(sce)
#' res <- pbDS(pb, method = "limma-trend")
#'
#' names(res)
#' names(res$table)
#' head(res$table$`stim-ctrl`$`B cells`)
#' 
#' # count nb. of DE genes by cluster
#' vapply(res$table$`stim-ctrl`, function(u) 
#'   sum(u$p_adj.loc < 0.05), numeric(1))
#' 
#' # get top 5 hits for ea. cluster w/ abs(logFC) > 1
#' library(dplyr)
#' lapply(res$table$`stim-ctrl`, function(u)
#'   filter(u, abs(logFC) > 1) %>% 
#'     arrange(p_adj) %>% 
#'     slice(seq_len(5)))
#'
#' @author Helena L Crowell & Mark D Robinson
#'
#' @references 
#' Crowell, HL, Soneson, C, Germain, P-L, Calini, D, 
#' Collin, L, Raposo, C, Malhotra, D & Robinson, MD: 
#' On the discovery of population-specific state transitions from 
#' multi-sample multi-condition single-cell RNA sequencing data. 
#' \emph{bioRxiv} \strong{713412} (2018). 
#' doi: \url{https://doi.org/10.1101/713412}
#'
#' @importFrom DESeq2 DESeq results
#' @importFrom edgeR calcNormFactors DGEList 
#'   estimateDisp glmQLFit glmQLFTest topTags
#' @importFrom dplyr last rename
#' @importFrom limma makeContrasts contrasts.fit eBayes lmFit topTable voom
#' @importFrom stats model.matrix
#' @importFrom SummarizedExperiment colData
#' @importFrom tibble add_column
#' @export

pbDS <- function(pb, 
    method = c("edgeR", "DESeq2", "limma-trend", "limma-voom"),
    design, coef, contrast, min_cells = 10, verbose = TRUE) {

    # check validity of input arguments
    method <- match.arg(method)
    .check_pbs(pb, check_by = TRUE)
    .check_args_pbDS(as.list(environment()))
    
    if (missing("design")) {
        formula <- ~ 0 + group_id
        cd <- as.data.frame(colData(pb))
        design <- model.matrix(formula, cd)
        colnames(design) <- levels(pb$group_id)
    }
    
    if (missing("coef") && missing("contrast")) {
        c <- colnames(design)[c(ncol(design), 1)]
        c <- paste(c, collapse = "-")
        contrast <- makeContrasts(contrasts = c, levels = design)
        coef <- NULL
    }
    
    if (!is.null(contrast)) {
        ct <- "contrast"
        names(cs) <- cs <- colnames(contrast)
    } else {
        ct <- "coef"
        cs <- vapply(coef, function(i) 
            paste(colnames(design)[i], collapse = "-"),
            character(1))
        names(cs) <- names(coef) <- cs
    }
    
    # compute cluster-sample counts
    n_cells <- metadata(pb)$n_cells
    kids <- assayNames(pb)
    names(kids) <- kids
    
    # for ea. cluster, run DEA
    res <- lapply(kids, function (k) {
        if (verbose) cat(k, "..", sep = "")
        y <- assays(pb)[[k]]
        rmv <- n_cells[k, ] < min_cells
        y <- y[, !rmv]
        d <- design[colnames(y), ]
        if (method == "DESeq2") {
            mode(y) <- "integer"
            cd <- colData(pb)[!rmv, , drop = FALSE]
            y <- DESeqDataSetFromMatrix(y, cd, d)
            y <- suppressMessages(DESeq(y))
            tt <- lapply(cs, function(c) {
                res <- results(y, contrast = contrast[, c])
                .res_df(k, res, ct, c) %>% 
                    rename(logFC = "log2FoldChange",
                        p_val = "pvalue", p_adj.loc = "padj")
            })
        } else {
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
                    .res_df(k, tt, ct, c) %>% 
                        rename(p_val = "PValue", p_adj.loc = "FDR")
                })
            } else {
                if (method == "limma-trend") {
                    trend <- robust <- TRUE
                } else if (method == "limma-voom") {
                    trend <- robust <- FALSE
                    y <- suppressMessages(DGEList(y, remove.zeros = TRUE))
                    y <- calcNormFactors(y)
                    y <- voom(y, d)
                }
                w <- n_cells[k, !rmv]
                fit <- lmFit(y, d, weights = w)
                tt <- lapply(cs, function(c) {
                    cfit <- contrasts.fit(fit, contrast[, c], coef[[c]])
                    efit <- eBayes(cfit, trend = trend, robust = robust)
                    tt <- topTable(efit, number = Inf, sort.by = "none")  
                    .res_df(k, tt, ct, c) %>% 
                        rename(p_val = "P.Value", p_adj.loc = "adj.P.Val")
                })
            }
        }
        return(list(tt = tt, data = y))
    })
    # remove empty clusters
    skipped <- vapply(res, is.null, logical(1))
    if (any(skipped) & verbose)
        message(paste("Cluster(s)", dQuote(kids[skipped]), "skipped due to an",
            "insufficient number of cells in at least 2 samples per group."))
    res <- res[!skipped]
    kids <- kids[names(res)]
    
    # re-organize by contrast & 
    # do global p-value adjustment
    tt <- map(res, "tt")
    if (!is.null(cs)) {
        tt <- lapply(cs, map, .x = tt)
        tt <- .p_adj_global(tt)
    } else {
        tt <- lapply(tt, function(u) 
            add_column(u, .after = "p_adj.loc",
                p_adj.glb = p.adjust(u$p_adj.loc)))
    }

    # return results
    list(table = tt, data = map(res, "data"), method = method, 
        design = design, contrast = contrast, coef = coef)
}
