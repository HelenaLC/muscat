#' @rdname runMAST
#' @title Run DS analysis using MAST
#' 
#' @description ...
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param formula model formula.
#' @param contrast a contrast matrix 
#'   created with \code{\link[limma]{makeContrasts}}.
#' @param assay character string specifying 
#'   which assay to use as input data.
#' 
#' @return a list containing \itemize{
#' \item \code{table}: a list of differential testing results 
#' for each contrast and cluster, and 
#' \item the \code{formula} and \code{contrast} used.}
#' 
#' @examples 
#' 
#' @author Helena L. Crowell \email{helena.crowell@uzh.ch}
#' 
#' @importFrom data.table data.table
#' @importFrom MAST FromMatrix lrTest SceToSingleCellAssay zlm
#' @importFrom methods is
#' @importFrom purrr map
#' @importFrom stats p.adjust
#' @importFrom SummarizedExperiment assays colData colData<- rowData rowData<- 
#' @importFrom tibble add_column
#' @export

runMAST <- function(x, formula, contrast, assay = "logcpm") {
    
    # validity checks
    .check_sce(x, req_group = TRUE)
    .check_arg_assay(x, assay)
    
    cs <- colnames(contrast)
    names(cs) <- cs
    stopifnot(!is.null(cs), is(contrast, "matrix"),  
        length(cs) == length(unique(cs)))
    
    # split cells by cluster
    cells_by_cluster <- .split_cells(x, by = "cluster_id")
    
    kids <- levels(x$cluster_id)
    names(kids) <- kids
    
    # need this for MAST to be happy
    colData(x)$wellKey <- colnames(x)
    rowData(x)$primerid <- rownames(x)
    
    tt <- lapply(kids, function(k) {
        if (verbose) cat(k, "..", sep = "")
        y <- x[, cells_by_cluster[[k]]]
        sca <- SceToSingleCellAssay(y, check_sanity = FALSE)
        suppressMessages(fit <- zlm(formula, sca))
        lapply(cs, function(c) {
            cm <- as.matrix(contrast[, c])
            suppressMessages(lrt <- lrTest(fit, cm))
            p_val <- lrt[, "hurdle", "Pr(>Chisq)"]
            p_adj.loc <- p.adjust(p_val)
            data.frame(gene = rownames(x), cluster_id = k, p_val, p_adj.loc,
                contrast = c, row.names = NULL, stringsAsFactors = FALSE)
        })
    })
    # re-organize by contrast & 
    # do global p-value adjustment
    tt <- lapply(cs, function(c) map(tt, c))
    tt <- .p_adj_global(tt)
    
    # return results
    list(table = tt,
        formula = formula,
        contrast = contrast)
}
