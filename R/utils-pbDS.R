#' @importFrom Matrix rowMeans rowSums
#' @importFrom matrixStats rowMedians
#' @importFrom purrr map_depth
#' @importFrom SummarizedExperiment assays
.pb <- function(x, cs, assay, fun) {
    fun <- switch(fun,
        rowSums = Matrix::rowSums,
        rowMeans = Matrix::rowMeans,
        rowMedians = function(u) 
            matrixStats::rowMedians(as.matrix(u)))
    pb <- map_depth(cs, -1, function(i) {
        if (length(i) == 0) return(numeric(nrow(x)))
        fun(assays(x)[[assay]][, i, drop = FALSE])
    })
    map_depth(pb, -2, function(u)
        as.matrix(data.frame(u,
            row.names = rownames(x),
            check.names = FALSE)))
}

# wrapper to create output tables
#   k:  cluster ID
#   tt: topTable data.frame
#' @importFrom tibble add_column
.res_df <- function(tbl, k, ct, c) {
    df <- data.frame(
        gene = rownames(tbl), cluster_id = k, tbl,
        row.names = NULL, stringsAsFactors = FALSE)
    df[[ct]] <- c; df
}

#' @importFrom DESeq2 DESeq results
#' @importFrom dplyr rename
#' @importFrom edgeR calcNormFactors DGEList estimateDisp
#'   filterByExpr glmQLFit glmQLFTest topTags
#' @importFrom scater isOutlier
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors metadata
.edgeR <- function(x, k, design, coef, contrast, ct, cs) {
    y <- assay(x, k)
    y <- suppressMessages(DGEList(y, 
        group = x$group_id[colnames(y)], 
        remove.zeros = TRUE))
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    # tests <- list(
    #     qlt = glmQLFTest, # genewise NB GLM with quasi-likelihood test
    #     treat = glmTreat) # test for DE relative to logFC threshold
    # tbl <- lapply(cs, function(c) {
    #     fits <- lapply(tests, function(fun) fun(fit, coef[[c]], contrast[, c]))
    #     tbls <- lapply(fits, topTags, n = Inf, sort.by = "none")
    #     # combine tables & reformat
    #     tbl_qlt <- tbls$qlt$table; tbl_treat <- tbls$treat$table
    #     tbl <- rename(tbl_qlt, p_val.qlt = "PValue", p_adj.qlt = "FDR")
    #     tbl$p_val.treat <- tbl_treat$PValue
    #     tbl$p_adj.treat <- tbl_treat$FDR
    #     tbl <- add_column(tbl, .after = "logFC", 
    #         unshrunk.logFC = tbl_treat$unshrunk.logFC)
    #     tbl <- .res_df(tbl, k, ct, c)
    # })
    tbl <- lapply(cs, function(c) {
        fit <- glmQLFTest(fit, coef[[c]], contrast[, c], )
        tbl <- topTags(fit, number = Inf, sort.by = "none")
        tbl <- .res_df(tbl, k, ct, c)
        rename(tbl, p_val = "PValue", p_adj.loc = "FDR")
    })
    list(table = tbl, data = y, fit = fit)
}

#' @importFrom dplyr rename
#' @importFrom edgeR calcNormFactors DGEList
#' @importFrom limma contrasts.fit eBayes lmFit topTable voom
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors metadata
.limma <- function(x, k, design, coef, contrast, ct, cs, method) {
    y <- assay(x, k)
    trend <- robust <- TRUE
    if (method == "voom") {
        trend <- robust <- FALSE
        y <- suppressMessages(DGEList(y, remove.zeros = TRUE))
        y <- calcNormFactors(y)
        y <- voom(y, design)
    } 
    w <- metadata(x)$n_cells[k, colnames(x)]   
    fit <- lmFit(y, design, weights = w)
    # tests <- list(
    #     bayes = eBayes, # eBayes moderated t-stat testing each contrast equal to 0
    #     treat = treat)  # eBayes moderated-t p-val relative to min logFC threshold
    # tbl <- lapply(cs, function(c) {
    #     fit <- contrasts.fit(fit, contrast[, c], coef[[c]])
    #     fits <- lapply(tests, function(fun) 
    #         fun(fit, trend = trend, robust = robust))
    #     tbls <- list(
    #         bayes = topTable(fits$bayes, number = Inf, sort.by = "none"),
    #         treat = topTreat(fits$treat, number = Inf, sort.by = "none"))
    #     tbl <- rename(tbls$bayes, t.bayes = "t",
    #         p_val.bayes = "P.Value", p_adj.bayes = "adj.P.Val")
    #     tbl$p_val.treat <- tbls$treat$P.Value
    #     tbl$p_adj.treat <- tbls$treat$adj.P.Val
    #     tbl$t.treat <- tbls$treat$t
    #     tbl <- .res_df(tbl, k, ct, c)
    # })
    tbl <- lapply(cs, function(c) {
        fit <- contrasts.fit(fit, contrast[, c], coef[[c]])
        fit <- topTable(fit, number = Inf, sort.by = "none")
        tbl <- .res_df(tbl, k, ct, c)
        rename(tbl, p_val = "P.Value", p_adj.loc = "adj.P.Val")
    })
    list(table = tbl, data = y, fit = fit)
}

.limma_trend <- function(x, k, design, coef, contrast, ct, cs)
    .limma(x, k, design, coef, contrast, ct, cs, method = "trend")

.limma_voom <- function(x, k, design, coef, contrast, ct, cs)
    .limma(x, k, design, coef, contrast, ct, cs, method = "voom")
  
#' @importFrom dplyr rename
#' @importFrom DESeq2 DESeq DESeqDataSetFromMatrix results 
#' @importFrom SummarizedExperiment assay colData
.DESeq2 <- function(x, k, design, contrast, ct, cs) {
    cd <- colData(x)
    y <- assay(x, k)
    mode(y) <- "integer"
    y <- DESeqDataSetFromMatrix(y, cd, design)
    y <- suppressMessages(DESeq(y))
    tbl <- lapply(cs, function(c) {
        tbl <- results(y, contrast[, c])
        tbl <- .res_df(tbl, k, ct, c)
        rename(tbl, logFC = "log2FoldChange", 
            p_val = "pvalue", p_adj.loc = "padj")
    })
    list(table = tbl, data = y)
}

