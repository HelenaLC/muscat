#' @rdname prepSCE
#' @title Prepare SCE for DS analysis
#' 
#' @description ...
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param cluster_id,sample_id,group_id character strings specifying
#'   the \code{colData(x)} columns containing cluster assignments,
#'   unique sample identifiers, and group IDs (e.g., treatment).
#' @param drop logical. Specifies whether \code{colData(x)} columns
#'   besides those specified as \code{cluster_id,sample_id,group_id}
#'   should be retained (default \code{drop = FALSE}) 
#'   or removed (\code{drop = TRUE}). 
#' 
#' @examples
#' # generate random counts
#' ng <- 50
#' nc <- 200
#' counts <- matrix(sample(ng * nc), nrow = ng, ncol = nc)
#'     
#' # generate some cell metadata
#' gids <- sample(c("groupA", "groupB"), nc, TRUE)   
#' sids <- sample(paste0("sample", seq_len(3)), nc, TRUE) 
#' kids <- sample(paste0("cluster", seq_len(5)), nc, TRUE) 
#' batch <- sample(seq_len(3), nc, TRUE)
#' 
#' # construct SCE
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(
#'   assays = list(counts = counts),
#'   colData = data.frame(group = gids, id = sids, cluster = kids, batch))
#'     
#' # prep. for workflow
#' sce <- prepSCE(sce,
#'   group_id = "group",
#'   sample_id = "id", 
#'   cluster_id = "cluster")
#'     
#' head(colData(sce))
#' metadata(sce)$experiment_info
#' 
#' @return a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' 
#' @importFrom dplyr mutate_all
#' @importFrom magrittr set_colnames
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom SingleCellExperiment reducedDims SingleCellExperiment
#' @importFrom SummarizedExperiment assays colData rowData
#' @export

prepSCE <- function(x, cluster_id, sample_id, group_id, drop = FALSE) {
    
    stopifnot(is(x, "SingleCellExperiment"))

    ids <- as.list(match.call())
    ids <- ids[grep("id", names(ids))]
    ids <- unlist(ids)
    
    stopifnot(is.character(ids))
    stopifnot(all(ids %in% colnames(colData(x))))
    
    cd0 <- colData(x)
    cd <- cd0[ids] %>% 
        data.frame(check.names = FALSE) %>% 
        mutate_all(as.factor) %>% 
        set_colnames(names(ids))
    
    if (!drop)
        cd <- data.frame(cd,
            cd0[setdiff(colnames(cd0), ids)], 
            check.names = FALSE)
        
    # construct SCE
    sce <- SingleCellExperiment(
        assays = assays(x),
        colData = DataFrame(cd),
        rowData = rowData(x),
        reducedDims = reducedDims(x))
    
    # construct metadata
    ei <- .make_ei(sce)
    metadata(sce)$experiment_info <- ei
    
    return(sce)
}
