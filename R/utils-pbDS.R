#' @importFrom BiocParallel SerialParam
#' @importFrom scuttle summarizeAssayByGroup
#' @importFrom SummarizedExperiment assay colData
.pb <- function(x, by, assay, fun, BPPARAM = SerialParam()) {
  # compute pseudobulks
  suppressWarnings( 
    # temporarily suppressing warnings b/c 'median' 
    # warns about unspecified 'useNames' argument
  y <- summarizeAssayByGroup(x,
        assay.type = assay, 
        ids = (ids <- colData(x)[by]),
        statistics = fun,
        BPPARAM = BPPARAM))
  colnames(y) <- y[[by[length(by)]]]
  
  if (length(by) == 1) 
      return(assay(y))
  
  # reformat into one assay per 'by[1]'
  if (is.factor(ids <- y[[by[1]]]))
      ids <- droplevels(ids)
  is <- split(seq_len(ncol(y)), ids)
  ys <- lapply(is, \(.) assay(y)[, ., drop=FALSE])
  
  # fill in missing combinations
  for (i in seq_along(ys)) {
      fill <- setdiff(
          unique(y[[by[2]]]), 
          colnames(ys[[i]]))
      if (length(fill != 0)) {
          foo <- matrix(0, nrow(x), length(fill))
          colnames(foo) <- fill
          foo <- cbind(ys[[i]], foo)
          o <- paste(sort(unique(y[[by[2]]])))
          ys[[i]] <- foo[, o]
      }
  }
  return(ys)
}

# extract table of cell counts from 'int_colData'
# of pseudobulks as returned by 'aggregateData'
#' @importFrom S4Vectors metadata
#' @importFrom SingleCellExperiment int_colData
.n_cells <- function(x) {
    y <- int_colData(x)$n_cells
    if (is.null(y)) return(NULL)
    if (length(metadata(x)$agg_pars$by) == 2)
        y <- as.matrix(data.frame(y, check.names = FALSE))
    return(as.table(y))
}

# wrapper to create output tables
#   k:  cluster ID
#   tt: topTable data.frame
.res_df <- function(tbl, k, ct, c) {
    df <- data.frame(
        gene = rownames(tbl), cluster_id = k, tbl,
        row.names = NULL, stringsAsFactors = FALSE)
    df[[ct]] <- c; df
}

#' @importFrom dplyr rename
#' @importFrom edgeR normLibSizes DGEList estimateDisp
#'   filterByExpr glmQLFit glmQLFTest glmTreat topTags
#' @importFrom scater isOutlier
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors metadata
.edgeR <- function(x, k, design, coef, contrast, ct, cs, treat) {
    y <- assay(x, k)
    y <- suppressMessages(DGEList(y, 
        group = x$group_id[colnames(y)], 
        remove.zeros = TRUE))
    y <- normLibSizes(y)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    # treat: test for DE relative to logFC threshold
    # else:  genewise NB GLM with quasi-likelihood test
    .fun <- ifelse(treat, glmTreat, glmQLFTest)
    tbl <- lapply(cs, function(c) {
        fit <- .fun(fit, coef[[c]], contrast[, c])
        tbl <- topTags(fit, n = Inf, sort.by = "none")
        # combine tables & reformat
        tbl <- rename(tbl$table, p_val = "PValue", p_adj.loc = "FDR")
        tbl <- .res_df(tbl, k, ct, c)
    })
    list(table = tbl, data = y, fit = fit)
}

#' @importFrom matrixStats rowMedians
.edgeR_NB <- \(x, k, design, coef, contrast, ct, cs, nc) {
    y <- assay(x, k)
    # Gene_level filtering to remove genes detected in
    # almost all cells of almost all pseudobulk samples
    med_detection <- rowMedians(sweep(y, 2, nc, "/"))
    gene_filter <- med_detection < 0.9
    # Normalization offset to remove systematic differences between pseudobulk 
    # samples that are due to technical or nuisance biological variability. 
    # Idea obtained from cellular detection rate (CDR) normalization from MAST. 
    # Note that this normalization is used instead of 'edgeR::normLibSizes()'.
    of <- colMeans(sweep(y[gene_filter, ], 2, nc, "/"))
    # construct 'DGEList'
    y <- suppressMessages(DGEList(
        counts = y[gene_filter, ],
        group = x$group_id[colnames(y)],
        remove.zeros = TRUE))
    # add offsets to 'DGEList'
    y$offset <- log(nc * of) 
    # run an 'edgeR' analysis
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust = TRUE)
    tbl <- lapply(cs, function(c) {
        fit <- glmQLFTest(fit,
            coef[[c]],
            contrast[, c],
            poisson.bound = FALSE) 
        tbl <- topTags(fit, n = Inf, sort.by = "none")
        tbl <- rename(tbl$table, p_val = "PValue", p_adj.loc = "FDR") 
        tbl <- .res_df(tbl, k, ct, c)
    })
    list(table = tbl, data = y, fit = fit)
}

#' @importFrom dplyr rename
#' @importFrom edgeR normLibSizes DGEList
#' @importFrom limma contrasts.fit eBayes lmFit topTable topTreat voom treat
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors metadata
.limma <- function(x, k, design, coef, contrast, ct, cs, method, treat) {
    y <- assay(x, k)
    trend <- robust <- TRUE
    if (method == "voom") {
        trend <- robust <- FALSE
        y <- suppressMessages(DGEList(y, remove.zeros = TRUE))
        y <- normLibSizes(y)
        y <- voom(y, design)
    } 
    w <- .n_cells(x)[k, colnames(x)]   
    fit <- lmFit(y, design, weights = w)
    # treat: eBayes moderated-t p-val relative to min logFC threshold
    # else:  eBayes moderated t-stat testing each contrast equal to 0 
    .fun <- ifelse(treat, limma::treat, eBayes)
    .tbl <- ifelse(treat, topTreat, topTable)
    tbl <- lapply(cs, function(c) {
        fit <- contrasts.fit(fit, contrast[, c], coef[[c]])
        fit <- .fun(fit, trend = trend, robust = robust)
        tbl <- .tbl(fit, number = Inf, sort.by = "none")
        tbl <- rename(tbl, p_val = "P.Value", p_adj.loc = "adj.P.Val")
        tbl <- .res_df(tbl, k, ct, c)
    })
    list(table = tbl, data = y, fit = fit)
}

.limma_trend <- function(x, k, design, coef, contrast, ct, cs, treat)
    .limma(x, k, design, coef, contrast, ct, cs, method = "trend", treat)

.limma_voom <- function(x, k, design, coef, contrast, ct, cs, treat)
    .limma(x, k, design, coef, contrast, ct, cs, method = "voom", treat)
  
#' @importFrom SummarizedExperiment assay colData
.DESeq2 <- function(x, k, design, contrast, ct, cs) {
    if (!require("DESeq2", quietly=TRUE))
        stop("Install 'DESeq2' to use this method.")
    cd <- colData(x)
    y <- as.matrix(assay(x, k)); mode(y) <- "integer"
    y <- DESeq2::DESeqDataSetFromMatrix(y, cd, design)
    y <- suppressMessages(DESeq2::DESeq(y))
    tbl <- lapply(cs, function(c) {
        tbl <- DESeq2::results(y, contrast[, c])
        tbl <- .res_df(tbl, k, ct, c)
        old <- c("log2FoldChange", "pvalue", "padj")
        new <- c("logFC", "p_val", "p_adj.loc")
        idx <- match(old, names(tbl))
        names(tbl)[idx] <- new; tbl
    })
    list(table = tbl, data = y)
}

