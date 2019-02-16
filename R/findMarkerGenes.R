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
#' @param block,lfc,test.use,min.pct passed to \code{\link[scran]{findMarkers}}.
#' 
#' @importFrom methods is
#' @importFrom parallel mclapply
#' @importFrom scran findMarkers

# pairwiseDE <- function(x, assay = "logcounts", 
#     lfc = 1, block = NULL, n_cores = 1) {
#     
#     # check validity of input arguments
#     # .check_sce(x)
#     # .check_arg_assay(x, assay)
#     # stopifnot(
#     #     is.numeric(n_cores), length(n_cores) == 1, 
#     #     n_cores > 0, as.integer(n_cores) == n_cores)
#     
#     if (!is.null(block))
#         block <- colData(x)[[block]]
#     kids <- colData(x)$cluster_id
#     res <- pairwiseTTests(
#         assays(x)$logcounts, 
#         clusters = kids, block = block,
#         direction = "up", lfc = lfc)
#     
#     kids <- levels(kids)
#     names(kids) <- kids
#     
#     tbl <- res$statistics %>% map_depth(1, data.frame)
#     names(tbl) <- res$pairs$first
#     
#     c1 <- res$pairs$first
#     idx <- split(seq_along(c1), c1)
#     tbl <- lapply(idx, function(i) 
#         bind_rows(tbl[i]) %>% add_column(.before = "logFC", 
#             gene = rep(rownames(x), length(i)),
#             cluster_id = rep(res$pairs$second[i], each = nrow(x))) %>% 
#             rename(p_val = "p.value", p_adj = "FDR"))
#     
#     tbl <- lapply(tbl, function(i)
#         data.frame(
#             cluster_id = res$pairs$second[i], 
#             res$statistics[[i]]))
# })
