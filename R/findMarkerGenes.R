#' @rdname findMarkerGenes
#' @title Identification of marker-genes via pairwise DE analysis
#' 
#' @description Test for DE expression between all pairs of clusters.
#'   For each cluster, cluster-specific genes can then be identified 
#'   by filtering for genes that are DE against most other clusters. 
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param cluster_id a character string specifying the clustering to use.
#'   Must match a \code{x@meta.data} column.
#' @param mc.cores passed to \code{\link[parallel]{mclapply}}.
#' @param assay,block,logfc.threshold,test.use,min.pct 
#'   passed to \code{\link[Seurat]{FindMarkers}}.
#' 
#' @importFrom parallel mclapply
#' @importFrom Seurat FindMarkers
#' @export

findMarkerGenes <- function(x, cluster_id, mc.cores = 1, assay = "RNA", block = NULL,
                            logfc.threshold = 1, test.use = "wilcox", min.pct = 0.25) {
  
  stopifnot(class(x) == "Seurat")
  stopifnot(is.numeric(mc.cores), length(mc.cores) == 1, mc.cores > 0)
  stopifnot(is.character(assay), length(assay) == 1, assay %in% names(x@assays))
  stopifnot(is.null(block) | (is.character(block) & length(block) == 1 & block %in% names(x@meta.data)))
  stopifnot(is.numeric(logfc.threshold), length(logfc.threshold) == 1, logfc.threshold >= 0)
  stopifnot(is.character(test.use), length(test.use) == 1,
            test.use %in% c("wilcox", "bimod", "roc", "t", "negbinom",
                            "poisson", "LR", "MAST", "DESeq2"))
  stopifnot(is.numeric(min.pct), length(min.pct) == 1, min.pct > 0, min.pct < 1)
  
  if (!is.null(block))
    block <- x@meta.data[[block]]
  x@active.assay <- assay
  x@active.ident <- setNames(x@meta.data[[cluster_id]], rownames(x@meta.data))
  cluster_ids <- levels(x@active.ident)
  
  do.call(rbind, lapply(cluster_ids, function(c1) {
    do.call(rbind, parallel::mclapply(cluster_ids, function(c2) {
      if (c1 != c2) {
        tryCatch({
          f <- FindMarkers(x, ident.1 = c1, ident.2 = c2, assay = assay, block = block,
                           logfc.threshold = logfc.threshold, test.use = test.use, min.pct = min.pct, only.pos = TRUE)
          f$gene <- rownames(f)
          f$cluster1 <- c1
          f$cluster2 <- c2
          f
        }, error = function(e) NULL)
      }
    }, mc.preschedule = FALSE, mc.cores = mc.cores))
  }))
}
