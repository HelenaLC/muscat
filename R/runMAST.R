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
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot(is(formula, "formula"))
    stopifnot(is.character(assay), length(assay) == 1, assay %in% assayNames(x))
    stopifnot(is(contrast, "matrix"), !is.null(cs), length(cs) == length(unique(cs)))
    
    # split cells by cluster
    cells_by_cluster <- .split_cells(x, by = "cluster_id")
    
    cluster_ids <- levels(colData(x)$cluster_id)
    names(cluster_ids) <- cluster_ids
    n_clusters <- length(cluster_ids)
    
    # need this for MAST to be happy
    colData(x)$wellKey <- colnames(x)
    rowData(x)$primerid <- rownames(x)
    
    res <- lapply(cluster_ids, function(k) {
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
    p.val <- modify_depth(res, 2, "p.val")
    
    # p-value adjustment (across all test in ea. cluster)
    p.adj <- vapply(p.val, function(u) 
        p.adjust(unlist(u), "BH"),
        numeric(nrow(x) * n_clusters))
    
    # re-split by cluster
    p.adj <- apply(p.adj, 2, split, 
        rep(cluster_ids, each = nrow(x)))
    
    # insert adjusted p-values into results
    res <- lapply(cs, function(c) 
        lapply(cluster_ids, function(k) 
            res[[c]][[k]] %>% add_column(
                p.adj = p.adj[[c]][[k]], 
                .before = "contrast")))
    
    # return results
    list(table = res,
        formula = formula,
        contrast = contrast)
}
