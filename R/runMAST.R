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
#' data(kang)
#' sce <- simData(kang, n_genes = 200, n_cells = 100, 
#'     p_dd = c(0.8, 0, 0.2, 0, 0, 0), fc = 4)
#'     
#' library(edgeR)
#' library(SingleCellExperiment)
#' 
#' dge <- DGEList(counts = assay(sce))
#' dge <- calcNormFactors(dge)
#' cpm <- cpm(dge)
#' assays(sce)$logcpm <- log2(cpm+1)
#' 
#' formula <- ~ 0 + group_id
#' contrast <- limma::makeContrasts("B-A", levels = c("A", "B"))
#' 
#' res <- runMAST(sce, formula, contrast)
#' names(res)
#' names(res$table)       # one list per contrast
#' names(res$table$`B-A`) # one table per cluster
#' 
#' # access results for a specific contrast & cluster
#' head(res$table$`B-A`[[1]])
#' 
#' # filter results
#' tbl_fil <- purrr::modify_depth(res$table, 
#'     depth = 2, dplyr::filter, p.adj < 0.05)
#' max(res$table$`B-A`[[1]]$p.adj) # before
#' max(tbl$`B-A`[[1]]$p.adj)       # after
#' 
#' @author Helena L. Crowell \email{helena.crowell@uzh.ch}
#' 
#' @importFrom data.table data.table
#' @importFrom MAST FromMatrix lrTest SceToSingleCellAssay zlm
#' @importFrom methods is
#' @importFrom stats p.adjust
#' @importFrom SummarizedExperiment assays colData colData<- rowData rowData<- 
#' @importFrom tibble add_column
#' 
#' @export

runMAST <- function(x, formula, contrast, assay = "logcpm") {
    
    cs <- colnames(contrast)
    names(cs) <- cs
    
    # validity checks
    .check_sce(x, req_group = TRUE)
    .check_arg_assay(x, assay)
    stopifnot(is(contrast, "matrix"), !is.null(cs), length(cs) == length(unique(cs)))
    
    # split cells by cluster
    cells_by_cluster <- .split_cells(x, by = "cluster_id")
    
    kids <- levels(colData(x)$cluster_id)
    names(kids) <- kids
    nk <- length(kids)
    
    # need this for MAST to be happy
    colData(x)$wellKey <- colnames(x)
    rowData(x)$primerid <- rownames(x)
    
    res <- lapply(kids, function(k) {
        y <- x[, cells_by_cluster[[k]]]
        sca <- SceToSingleCellAssay(y, check_sanity = FALSE)
        suppressMessages(fit <- zlm(formula, sca))
        lapply(cs, function(c) {
            cm <- as.matrix(contrast[, c])
            suppressMessages(lrt <- lrTest(fit, cm))
            p.val <- lrt[, "hurdle", "Pr(>Chisq)"]
            data.frame(gene = rownames(x), cluster_id = k, p.val,
                contrast = c, row.names = NULL, stringsAsFactors = FALSE)
        })
    })
    # re-organize by contrast
    res <- lapply(cs, function(c) lapply(res, "[[", c))
    p_val <- modify_depth(res, 2, "p.val")
    
    # p-value adjustment (across all test in ea. cluster)
    p_adj <- vapply(p_val, function(u) 
        p.adjust(unlist(u), "BH"),
        numeric(nrow(x) * nk))
    
    # re-split by cluster
    p_adj <- apply(p_adj, 2, split, 
        rep(kids, each = nrow(x)))
    
    # insert adjusted p-values into results
    res <- lapply(cs, function(c) 
        lapply(kids, function(k) 
            res[[c]][[k]] %>% add_column(
                p_adj = p_adj[[c]][[k]], 
                .before = "contrast")))
    
    # return results
    list(table = res,
        formula = formula,
        contrast = contrast)
}
