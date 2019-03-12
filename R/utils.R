cluster_colors <- c(
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

# ==============================================================================
# scale values b/w 0 and 1 using 
# low (1%) and high (99%) quantiles as boundaries
# ------------------------------------------------------------------------------
#' @importFrom matrixStats rowQuantiles
.scale <- function(x) {
    qs <- rowQuantiles(as.matrix(x), probs = c(.01, .99))
    x <- (x - qs[, 1]) / (qs[, 2] - qs[, 1])
    x[x < 0] <- 0
    x[x > 1] <- 1
    return(x)
}

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
#' @importFrom utils getFromNamespace
.pb <- function(x, cs, assay, fun) {
    fun <- switch(fun,
        rowMedians = getFromNamespace(fun, "matrixStats"),
        getFromNamespace(fun, "Matrix"))
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
# global p-value adjustment
# ------------------------------------------------------------------------------
#   x: results table; a nested list w/ 
#      1st level = comparisons and 2nd level = clusters
# > adds 'p_adj.glb' column to the result table of ea. comparison & cluster
# ------------------------------------------------------------------------------
.p_adj_global <- function(x) {
    cs <- names(x)
    ks <- names(x[[1]])
    names(cs) <- cs
    names(ks) <- ks
    lapply(cs, function(c) {
        # get p-values
        p_val <- map(x[[c]], "p_val")
        # adjust for each comparison
        p_adj <- p.adjust(unlist(p_val))
        # re-split by cluster
        ns <- vapply(p_val, length, numeric(1))
        p_adj <- split(p_adj, rep.int(ks, ns))
        # insert into results tables
        lapply(ks, function(k) x[[c]][[k]] %>% add_column(
            p_adj.glb = p_adj[[k]], .after = "p_adj.loc"))
    })
}
