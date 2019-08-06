# ==============================================================================
# filters rows/columns from input matrix `m` until all entries >= `n`, 
# such that ea. iteration removes the row/column w/ the smallest summed value.
# ------------------------------------------------------------------------------
.filter_matrix <- function(m, n = 100) {
    while (any(m < n)) {
        # get candidate rows/cols for removal
        i <- m < n
        r <- apply(i, 1, any)
        c <- apply(i, 2, any)
        # get smallest row/col
        rs <- rowSums(m)
        cs <- colSums(m)
        r <- which(r)[which.min(rs[r])]
        c <- which(c)[which.min(cs[c])]
        # priorities removal of rows over cols
        if (rs[r] <= cs[c]) {
            m <- m[-r, , drop = FALSE]
        } else {
            m <- m[, -c, drop = FALSE]
        }
        if (any(dim(m) == 1)) 
            break
    }
    return(m)
}

.update_sce <- function(sce) {
    # update colData
    cd <- as.data.frame(colData(sce))
    cd <- mutate_if(cd, is.factor, droplevels) 
    colData(sce) <- DataFrame(cd, row.names = colnames(sce))
    # update metadata
    ei <- metadata(sce)$experiment_info
    ei <- ei[ei$sample_id %in% levels(sce$sample_id), ]
    ei <- mutate_if(ei, is.factor, droplevels)
    metadata(sce)$experiment_info <- ei
    return(sce)
}

.filter_sce <- function(sce, kids, sids) {
    cs1 <- sce$cluster_id %in% kids
    cs2 <- sce$sample_id %in% sids
    sce <- sce[, cs1 & cs2]
    sce <- .update_sce(sce)
    return(sce)
}

.cluster_colors <- c(
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

# ==============================================================================
# wrapper for z-normalization
# ------------------------------------------------------------------------------
.z_norm <- function(x, th = 2.5) {
    x <- as.matrix(x)
    sds <- rowSds(x, na.rm = TRUE)
    sds[sds == 0] <- 1
    x <- t(t(x - rowMeans(x, na.rm = TRUE)) / sds)
    #x <- (x - rowMeans(x, na.rm = TRUE)) / sds
    x[x >  th] <-  th
    x[x < -th] <- -th
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
    sids <- unique(x$sample_id)
    m <- match(sids, x$sample_id)
    df <- data.frame(
        stringsAsFactors = FALSE,
        sample_id = sids,
        group_id = x$group_id[m],
        n_cells = as.numeric(table(x$sample_id)[sids]))
    for (i in c("sample_id", "group_id"))
        if (is.factor(x[[i]]))
            df <- mutate_at(df, i, factor, levels = levels(x[[i]]))
    return(df)
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

# ------------------------------------------------------------------------------
# toy SCE for unit-testing
# ------------------------------------------------------------------------------
#' @importFrom SingleCellExperiment SingleCellExperiment
.toySCE <- function() {
    ngs <- 300
    ncs <- 2e3
    
    gs <- paste0("gene", seq_len(ngs))
    cs <- paste0("cell", seq_len(ncs))
    
    y <- sample(seq(0, 100, 1), ngs * ncs, TRUE)
    y <- matrix(y,
        nrow = ngs, ncol = ncs,
        dimnames = list(gs, cs))
    
    kids <- sample(paste0("k", seq_len(5)), ncs, TRUE)
    sids <- sample(paste0("s", seq_len(4)), ncs, TRUE)
    gids <- sample(paste0("g", seq_len(3)), ncs, TRUE)
    sids <- paste(sids, gids, sep = ".")
    
    cd <- data.frame(
        sample_id = sids,
        group_id = gids,
        cluster_id = kids)
    
    SingleCellExperiment(
        assay = list(counts = y), colData = cd,
        metadata = list(experiment_info = .make_ei(cd)))
}
