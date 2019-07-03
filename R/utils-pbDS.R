#' @importFrom Matrix rowMeans rowSums
#' @importFrom matrixStats rowMedians
#' @importFrom purrr map_depth
#' @importFrom SummarizedExperiment assays
#' @importFrom utils getFromNamespace
.pb <- function(x, cs, assay, fun) {
    fun <- switch(fun,
        rowMedians = function(u) 
            matrixStats::rowMedians(as.matrix(u)),
        getFromNamespace(fun, "Matrix"))
    pb <- map_depth(cs, -1, function(i) {
        if (length(i) == 0) return(numeric(nrow(x)))
        fun(assays(x)[[assay]][, i, drop = FALSE])
    })
    map_depth(pb, -2, function(u)
        as.matrix(data.frame(u,
            row.names = rownames(x),
            check.names = FALSE)))
}

# wrapper to create output tables
#   k:  cluster ID
#   tt: topTable data.frame
#   ct: comparison type; "contrast" or "coef"
#   c:  character string specifying the comparison
res_df <- function(k, tt, ct, c) {
    df <- data.frame(
        gene = rownames(tt), cluster_id = k, tt,
        row.names = NULL, stringsAsFactors = FALSE)
    df[[ct]] <- c
    return(df)
}