#' simDD
#' 
#' Simulation of complex scRNA-seq data 
#' 
#' \code{simDD} simulates multiple clusters and samples 
#' across 2 experimental conditions from a real scRNA-seq data set.
#' 
#' @param x a list of length 2 containing a
#'   \code{\linkS4class{SummarizedExperiment}} & a NB fit.
#' @param n_genes # of genes to simulate. 
#' @param n_cells # of cells to simulate. 
#'   Either a single numeric or a range to sample from.
#' @param ns nb. of genes common to 1, 2, ..., all clusters.
#' @param p_dd numeric vector of length 6.
#'   Specifies the probability of a gene being
#'   EE, EP, DE, DP, DM, or DB, respectively.
#' @param seed random seed. 
#' 
#' @examples
#' data(kang_se, kang_fit)
#' simDD(kang_se, kang_fit,
#'     n_genes = 10, n_cells = 10,
#'     p_dd = c(1,0,0,0,0,0), seed = 1)
#' 
#' @import SummarizedExperiment 
#' @importFrom data.table data.table
#' @importFrom edgeR DGEList estimateDisp glmFit
#' @importFrom stats model.matrix rnbinom
#' @importFrom S4Vectors split
#' @importFrom zeallot %<-%
#' 
#' @export

simDD <- function(x, fit = NULL,
    n_genes = 200, n_cells = c(50, 100), ns = NULL,
    p_dd = c(0.4, 0.35, 0.12, 0.08, 0.04, 0.01), seed = 1) {
    
    fit <- x[[2]]
    x <- x[[1]]
    
    # validty checks
    stopifnot(class(x) == "SummarizedExperiment")
    stopifnot(is.numeric(n_genes), length(n_genes) == 1)
    stopifnot(is.numeric(n_cells), length(n_cells) == 1 | length(n_cells) == 2)
    stopifnot(is.numeric(p_dd), length(p_dd) == 6, sum(p_dd) == 1)
    stopifnot(is.numeric(seed), length(seed) == 1)
    
    cluster_ids <- levels(factor(colData(x)$cluster_id))
    sample_ids <- levels(factor(colData(x)$sample_id))
    
    n_clusters <- length(cluster_ids)
    n_samples <- length(sample_ids)
    
    set.seed(seed)
    
    # split cells by cluster-sample
    dt <- data.table(
        t(assays(x)$counts), 
        cell = colnames(x), 
        cluster_id = colData(x)$cluster_id, 
        sample_id = colData(x)$sample_id)
    dt <- split(dt, by = c("cluster_id", "sample_id"), 
        sorted = TRUE, keep.by = FALSE, flatten = FALSE)
    cells <- sapply(dt, sapply, "[[", "cell")
    
    # get nb. of cells to simulate per cluster-sample
    if (length(n_cells) == 1) {
        n_cells <- list(rep(n_cells, 2))
    } else {
        n_cells <- replicate(n_clusters * n_samples, 
            list(sample(n_cells[1]:n_cells[2], 2)))
    }
    n_cells <- matrix(n_cells, 
        nrow = n_samples, ncol = n_clusters, 
        dimnames = list(sample_ids, cluster_ids))
    
    # sample nb. of genes to simulate per category
    cats <- c("ee", "ep", "de", "dp", "dm", "db")
    ndd <- replicate(n_clusters, {
        ns <- sample(cats, n_genes, replace = TRUE, prob = p_dd)
        factor(ns, levels = cats) 
        }, simplify = FALSE)
    ndd <- sapply(ndd, table)
    colnames(ndd) <- cluster_ids
    
    # initialize count matrix
    y <- matrix(0, 
        nrow = n_genes, 
        ncol = sum(unlist(n_cells)),
        dimnames = list(
            paste0("gene", seq_len(n_genes)),
            paste0("cell", seq_len(sum(unlist(n_cells))))))
    
    # sample gene indices
    is <- sapply(cluster_ids, function(c, gs = rownames(y))
        sapply(cats, function(cat) { 
            n <- ndd[cat, c]
            x <- sample(gs, n)
            gs <<- setdiff(gs, x)
            return(x) }))
    
    # sample cell indices
    cs <- colnames(y)
    js <- sapply(cluster_ids, function(c) 
        setNames(lapply(sample_ids, function(s)
            lapply(n_cells[[s, c]], function(n) {
                x <- sample(cs, n)
                cs <<- setdiff(cs, x)
                return(x) })), sample_ids))
    colnames(is) <- colnames(js) <- cluster_ids

    # sample genes to simulate from
    if (is.null(ns)) {
        ns <- sample(n_genes, n_clusters)
        while (sum(ns) != n_genes)
            ns <- floor(ns / sum(ns) * n_genes) + sample(c(0, 1), n_clusters, replace = TRUE)
    } 

    gs <- rownames(x)
    gs_used_for_sim <- matrix( 
        sample(gs, n_clusters * n_genes, replace = TRUE), 
        nrow = n_genes, ncol = n_clusters, 
        dimnames = list(rownames(y), cluster_ids))
    idx_gs <- c(0, cumsum(ns))
    idx_gs <- sapply(seq_len(n_clusters), function(i) (idx_gs[i]+1):idx_gs[i+1])
    idx_cs <- sapply(seq_len(n_clusters), function(n) sample(n_clusters, n))
    for (n in seq_len(n_clusters)) {
        if (ns[n] != 0) {
            i <- idx_gs[[n]]
            j <- idx_cs[[n]]
            # sample common genes
            gs_n <- sample(gs, ns[n], replace = TRUE)
            # store genes used for simulation
            gs_used_for_sim[i, j] <- replicate(n, gs_n)
        }
    }
    # count genes overlapping between 1, 2, ..., all clusters
    #table(apply(gs_used_for_sim, 1, function(gs) max(table(gs))))

    for (c in cluster_ids) {
        # genes to simulate from
        gs <- gs_used_for_sim[, c]
        
        # get NB parameters
        m <- fit$coefficients[gs, 1]
        d <- setNames(fit$dispersion, rownames(x))[gs]
        d <- 1 / d
        
        # get gene indices & nb. of genes by category
        c(iee, iep, ide, idp, idm, idb) %<-% sapply(cats, function(i) is [i, c])
        c(nee, nep, nde, ndp, ndm, ndb) %<-% vapply(cats, function(i) ndd[i, c], numeric(1))
        
        for (s in sample_ids) {
            # cells to simulate from
            cs <- cells[[s, c]]

            # compute mus
            o <- setNames(fit$offset[1, ], colnames(x))[cs]
            mu <- sapply(exp(o), "*", exp(m))
            
            # get cell indices & nb. of cells by group
            ng1 <- length(g1 <- js[[s, c]][[1]])
            ng2 <- length(g2 <- js[[s, c]][[2]])

            # simulate data
            if (nee > 0) y[iee, c(g1, g2)] <- simdd("ee", gs[iee, c], cs, ng1, ng2, mu, d)
            if (nep > 0) y[iep, c(g1, g2)] <- simdd("ep", gs[iep, c], cs, ng1, ng2, mu, d)
            if (nde > 0) y[ide, c(g1, g2)] <- simdd("de", gs[ide, c], cs, ng1, ng2, mu, d)
            if (ndp > 0) y[idp, c(g1, g2)] <- simdd("dp", gs[idp, c], cs, ng1, ng2, mu, d)
            if (ndm > 0) y[idm, c(g1, g2)] <- simdd("dm", gs[idm, c], cs, ng1, ng2, mu, d)
            if (ndb > 0) y[idb, c(g1, g2)] <- simdd("db", gs[idb, c], cs, ng1, ng2, mu, d)
        }
    }
    
    # construct SummarizedExperiment -------------------------------------------

    gene_info <- lapply(cluster_ids, function(c) 
        do.call(rbind, lapply(cats, function(cat) {
            genes <- unlist(is[[cat, c]])
            if (length(genes != 0))
                data.frame(gene = genes, cluster_id = c, category = cat)
        })))
    gene_info <- do.call(rbind, gene_info)
    gene_info <- gene_info[order(as.numeric(gsub("[a-z]", "", gene_info$gene))), ]
    rownames(gene_info) <- NULL
    
    row_data <- data.frame(
        row.names = rownames(y), 
        marker_name = rownames(y),
        marker_class = factor("state", levels = c("none", "type", "state")))
    
    col_data <- do.call(rbind, lapply(cluster_ids, function(c) 
        do.call(rbind, lapply(sample_ids, function(s) 
            data.frame(row.names = 1,
                unlist(js[[s, c]]), 
                cluster_id = c, sample_id = s, 
                group = rep.int(c("A", "B"), n_cells[[s, c]]))))))
    col_data <- col_data[colnames(y), ]
    col_data$sample_id <- paste0(col_data$group, col_data$sample_id)
    
    ei <- data.frame(
        sample_id = paste0(rep(c("A", "B"), each = length(sample_ids)), sample_ids),
        group = rep(c("A", "B"), each = length(sample_ids)))
    md <- list(
        experiment_info = ei,
        n_cells = table(col_data$sample_id),
        gene_info = gene_info, 
        sim_genes = gs_used_for_sim)
    
    SummarizedExperiment(assays = list(counts = y),
        rowData = row_data, colData = col_data, metadata = md)
}