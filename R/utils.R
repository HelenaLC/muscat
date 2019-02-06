

# ------------------------------------------------------------------------------
# generate experimental design metadata table 
# for an input SCE or colData data.frame
# ------------------------------------------------------------------------------
#' @importFrom SummarizedExperiment colData
.make_ei <- function(x) {
    if (is(x, "SingleCellExperiment"))
        x <- colData(x)
    m <- match(levels(x$sample_id), x$sample_id)
    data.frame(
        sample_id = levels(x$sample_id),
        group_id = x$group_id[m],
        n_cells = as.numeric(table(x$sample_id)))
}

# ------------------------------------------------------------------------------
# compute pseudo-bulks
# ------------------------------------------------------------------------------
#' @importFrom purrr map_depth
#' @importFrom SummarizedExperiment assays
.pb <- function(x, cs, assay, fun) {
    fun <- getFromNamespace(fun, "Matrix")
    pb <- map_depth(cs, -1, function(i) {
        if (length(i) == 0) return(numeric(nrow(x)))
        fun(assays(x)[[assay]][, i, drop = FALSE])
    })
    map_depth(pb, -2, function(u) 
        data.frame(u, 
            row.names = rownames(x),
            check.names = FALSE))
}

# ------------------------------------------------------------------------------
# split cells by cluster-sample
# ------------------------------------------------------------------------------
#   x:  a SingleCellExperiment or colData
#   by: character vector specifying colData column(s) to split by
# > If length(by) == 1, a list of length nlevels(colData$by), else,
#   a nested list with 2nd level of length nlevels(colData$by[2])
# ------------------------------------------------------------------------------
#' @importFrom data.table data.table
#' @importFrom purrr map_depth
.split_cells <- function(x, 
    by = c("cluster_id", "sample_id")) {
    if (is(x, "SingleCellExperiment"))
        x <- colData(x)
    cd <- data.frame(x[by], check.names = FALSE)
    cd <- data.table(cd, cell = rownames(cd)) %>% 
        split(by = by, sorted = TRUE, flatten = FALSE)
    map_depth(cd, length(by), "cell")
}

# ------------------------------------------------------------------------------
# removes lowest contribution row/column until all entries >= n
# used by prepSim() to filter clusters & samples to use for simulation
# ------------------------------------------------------------------------------
.filter_matrix <- function(m, n = 100) {
    while (any(m < n)) {
        s <- sum(m)
        rows <- rowSums(m) / s
        cols <- colSums(m) / s
        x <- c(rows, cols)
        rmv <- names(which.min(x))
        y <- TRUE
        while (y) {
            if (rmv %in% rownames(m)) {
                if (all(m[rmv, ] >= n)) {
                    x <- x[names(x) != rmv]
                    rmv <- names(which.min(x)) 
                } else {
                    m <-  m[rownames(m) != rmv, ]
                    y <- FALSE
                }
            } else {
                if (all(m[, rmv] >= n)) {
                    x <- x[names(x) != rmv]
                    rmv <- names(which.min(x))
                } else {
                    m <- m[, colnames(m) != rmv]
                    y <- FALSE
                }
            }
        }
    }
    return(m)
}








