#' runMAST
#' Run MAST
#' 
#' ...
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' 
#' @examples 
#' data(kang)
#' sce <- simData(kang, n_genes = 200, n_cells = 200, 
#'     p_dd = c(0.8, 0, 0.2, 0, 0, 0), fc = 4)
#'     
#' library(edgeR)
#' library(SingleCellExperiment)
#' dge <- DGEList(counts = assay(sce))
#' dge <- calcNormFactors(dge)
#' cpm <- cpm(dge)
#' assays(sce)$logcpm <- log2(cpm+1)
#' 
#' @importFrom data.table data.table split
#' @importFrom MAST FromMatrix lrTest zlm.SingleCellAssay
#' @importFrom methods is
#' @importFrom SummarizedExperiment assays colData
#' 
#' @export

runMAST <- function(x, formula, contrast,  assay = "logcpm") {
    
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot(assay %in% assayNames(x))
    
    # split by cluster
    dt <- data.table(data.frame(colData(x)), cell = colnames(x))
    dt_split <- split(dt, by = "cluster_id", flatten = FALSE)
    cells_by_cluster <- lapply(dt_split, "[[", "cell")
    
    cluster_ids <- levels(colData(x)$cluster_id)
    names(cluster_ids) <- cluster_ids
    
    colData(x)$wellKey <- colnames(x)
    rowData(x)$primerid <- rownames(x)
    
    cs <- colnames(contrast)
    names(cs) <- cs
    
    res <- lapply(cluster_ids, function(k) {
        sce <- x[, cells_by_cluster[[k]]]
        sca <- SceToSingleCellAssay(sce, check_sanity = FALSE)
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
    lapply(cs, function(c) lapply(res, "[[", c))
}
# x <- sce0
# groups <- colData(x)$group_id
# groups <- factor(groups, levels = c("A", "B", "C"))
# groups[1:10] <- "C"
# colData(x)$group_id <- groups
# design <- model.matrix(~0 + groups)
# colnames(design) <- levels(groups)
# contrast <- makeContrasts(contrasts = c("B-A", "C-A"), levels = design)
# res <- runMAST(x, ~group_id, contrast)
# names(res)
# names(res[[1]])
# head(res[[1]])
