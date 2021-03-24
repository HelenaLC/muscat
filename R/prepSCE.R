#' @rdname prepSCE
#' @title Prepare SCE for DS analysis
#' 
#' @description ...
#' 
#' @param x a \linkS4class{SingleCellExperiment}.
#' @param kid,sid,gid character strings specifying
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
#'     
#' # generate some cell metadata
#' gids <- sample(c("groupA", "groupB"), nc, TRUE)   
#' sids <- sample(paste0("sample", seq_len(3)), nc, TRUE) 
#' kids <- sample(paste0("cluster", seq_len(5)), nc, TRUE) 
#' batch <- sample(seq_len(3), nc, TRUE)
#' cd <- data.frame(group = gids, id = sids, cluster = kids, batch)
#' 
#' # construct SCE
#' library(scuttle)
#' sce <- mockSCE(ncells = nc, ngenes = ng)
#' colData(sce) <- cbind(colData(sce), cd)
#' 
#' # prep. for workflow
#' sce <- prepSCE(sce, kid = "cluster", sid = "id", gid = "group")
#' head(colData(sce))
#' metadata(sce)$experiment_info
#' sce
#' 
#' @author Helena L Crowell
#' 
#' @return a \linkS4class{SingleCellExperiment}.
#' 
#' @importFrom dplyr mutate_all
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom SingleCellExperiment reducedDims SingleCellExperiment
#' @importFrom SummarizedExperiment assays colData rowData
#' @export

prepSCE <- function(x, 
    kid = "cluster_id", 
    sid = "sample_id", 
    gid = "group_id", 
    drop = FALSE) {

    stopifnot(is(x, "SingleCellExperiment"))

    args <- as.list(environment())
    ids <- args[grep("[a-z]id", names(args))]
    ids <- unlist(ids)
    
    stopifnot(is.character(ids))
    stopifnot(all(ids %in% colnames(colData(x))))
    
    cd0 <- colData(x)
    cd <- data.frame(cd0[ids], check.names = FALSE)
    cd <- mutate_all(cd, as.factor)
    colnames(cd) <- unlist(formals()[names(ids)])
    
    if (!drop)
        cd <- data.frame(cd,
            cd0[setdiff(colnames(cd0), ids)], 
            check.names = FALSE)

    # replace colData in SCE
    colData(x) <- DataFrame(cd)

    # construct metadata
    ei <- .make_ei(x)
    metadata(x)$experiment_info <- ei
    
    return(x)
}
