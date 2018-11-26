#' @rdname run_diffcyt
#' @title ...
#'
#' @description ...
#'
#' @param x a \code{daFrame}.
#' @param group a character string specifying the grouping variable.
#' @param data character string specifying the data to use.
#'   Should be one of \code{assayNames(x)}.
#' @param log logical. Should log2 be taken before computing means/medians.
#' @param method a character string.
#' @param fun character string specifying the summary statistic to use.
#'
#' @author Helena Lucia Crowell
#'
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @importFrom diffcyt calcCounts calcMedians testDS_limma testDS_LMM
#'   createContrast createDesignMatrix createFormula
#' @importFrom scater normalize
#' @importFrom S4Vectors metadata
#'
#' @export

run_diffcyt <- function(x, group = "group", k = NULL,
    data = c("counts", "normcounts"), log = FALSE,
    method = c("limma", "LMM"), fun = c("mean", "median")) {
    
    stopifnot(class(x) == "SummarizedExperiment")
    stopifnot(is.character(group), length(group) == 1, group %in% colnames(colData(x)))
    stopifnot(is.logical(log), length(log) == 1)
    data <- match.arg(data)
    method <- match.arg(method)
    fun <- match.arg(fun)
    
    # flip object dimensions such that 
    # rows = cells & columns = genes
    x <- SummarizedExperiment(
        assays = lapply(assays(x), t),
        rowData = colData(x),
        colData = rowData(x),
        metadata = metadata(x))
    
    cs <- calcCounts(x)
    
    if (data == "normcounts" & !"normcounts" %in% assayNames(x)) {
        sce <- SingleCellExperiment(assays = list(counts = assays(x)$counts))
        suppressWarnings(sce <- normalize(sce, return_log = FALSE))
        assays(x)$normcounts <- assays(sce)$normcounts
    }
    stopifnot(data %in% assayNames(x))
    assay(x) <- assays(x)[[data]]
    if (log) assay(x) <- log2(assay(x) + 1)
    
    ms <- switch(fun,
        mean = calcMeans(x),
        median = calcMedians(x))
    
    contrast <- createContrast(c(0, 1))
    sample_ids <- levels(factor(rowData(x)$sample_id))
    m <- match(sample_ids, rowData(x)$sample_id)
    groups <- rowData(x)[[group]][m]
    df <- data.frame(sample_ids, groups)
    
    res <- switch(method,
        limma = {
            design <- createDesignMatrix(df, cols_design = 2)
            testDS_limma(cs, ms, design, contrast)
        },
        LMM = {
            formula <- createFormula(df, cols_fixed = 2)
            testDS_LMM(cs, ms, formula, contrast)
         
        })
    
    # reformat output
    res <- rowData(res)
    colnames(res)[2] <- "gene"
    return(res)
}
