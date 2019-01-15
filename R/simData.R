#' simData
#' 
#' Simulation of complex scRNA-seq data 
#' 
#' \code{simData} simulates multiple clusters and samples 
#' across 2 experimental conditions from a real scRNA-seq data set.
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
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
#' data(kang)
#' simData(kang,
#'     n_genes = 10, n_cells = 10,
#'     p_dd = c(1,0,0,0,0,0), seed = 1)
#' 
#' @import SingleCellExperiment 
#' @importFrom data.table data.table
#' @importFrom edgeR DGEList estimateDisp glmFit
#' @importFrom stats model.matrix rnbinom setNames
#' @importFrom S4Vectors split
#' @importFrom zeallot %<-%
#' 
#' @export

simData <- function(x, n_genes, n_cells, p_dd, fc = 2, seed = 1) {
    
    # check validity of input arguments
    stopifnot(class(x) == "SingleCellExperiment")
    stopifnot(is.numeric(n_genes), length(n_genes) == 1)
    stopifnot(is.numeric(n_cells), length(n_cells) == 1 | length(n_cells) == 2)
    stopifnot(is.numeric(p_dd), length(p_dd) == 6, sum(p_dd) == 1)
    stopifnot(is.numeric(fc), is.numeric(fc), fc > 1)
    stopifnot(is.numeric(seed), length(seed) == 1)
    
    cluster_ids <- levels(colData(x)$cluster_id)
    sample_ids <- levels(colData(x)$sample_id)
    n_clusters <- length(cluster_ids)
    n_samples <- length(sample_ids)
    
    # split cells by cluster-sample
    dt <- data.table(
        cell = colnames(x), 
        cluster_id = colData(x)$cluster_id,
        sample_id = colData(x)$sample_id)
    dt_split <- split(dt,
        by = c("cluster_id", "sample_id"), 
        keep.by = FALSE, flatten = FALSE)
    cells_by_cluster_sample <- sapply(dt_split, sapply, "[[", "cell")
    
    # sample nb. of cells to simulate per cluster-sample
    if (length(n_cells) == 1) {
        n_cells <- list(rep(n_cells, 2))
    } else {
        n_cells <- replicate(n_clusters * n_samples, 
            list(sample(n_cells[1]:n_cells[2], 2)))
    }
    n_cells <- matrix(n_cells, 
        nrow = n_samples, ncol = n_clusters, 
        dimnames = list(sample_ids, cluster_ids))
    
    # initialize count matrix
    y <- matrix(0, 
        nrow = n_genes, 
        ncol = sum(unlist(n_cells)),
        dimnames = list(
            paste0("gene", seq_len(n_genes)),
            paste0("cell", seq_len(sum(unlist(n_cells))))))
    
    # sample nb. of genes to simulate per category
    ndd <- replicate(n_clusters, {
        ns <- sample(cats, n_genes, replace = TRUE, prob = p_dd)
        factor(ns, levels = cats) 
    }, simplify = FALSE)
    ndd <- sapply(ndd, table)
    colnames(ndd) <- cluster_ids
    
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
    
    # sample genes to simulate from
    gs <- replicate(n_clusters, sample(rownames(x), n_genes, replace = TRUE))
    rownames(gs) <- rownames(y)
    colnames(gs) <- cluster_ids
    
    # sample fold-changes
    fcs <- sapply(cluster_ids, function(k) 
        sapply(cats, function(c) { 
            n <- ndd[c, k]
            signs <- sample(c(-1, 1), size = n, replace = TRUE)
            fcs <- 2 ^ ( rgamma(n, 4, 4 / fc) * signs )
            names(fcs) <- gs[is[[c, k]], k]
            return(fcs)
    }))
    
    for (c in cluster_ids) {
        # get NB parameters
        m <- rowData(x)[gs[, c], ]$beta
        d <- rowData(x)[gs[, c], ]$dispersion
        names(m) <- names(d) <- gs[, c]
        
        for (s in sample_ids) {
            # cells to simulate from
            cs <- cells_by_cluster_sample[[s, c]]
            
            # compute mus
            o <- setNames(colData(x)[cs, ]$offset, cs)
            mu <- sapply(exp(o), "*", exp(m))
            
            # get cell indices & nb. of cells by group
            ng1 <- length(g1 <- js[[s, c]][[1]])
            ng2 <- length(g2 <- js[[s, c]][[2]])
            
            # simulate data
            for (cat in cats)
                if (ndd[cat, c] > 0) y[is[[cat, c]], c(g1, g2)] <- 
                simdd(cat, gs[is[[cat, c]], c], cs, ng1, ng2, mu, d, fcs[[cat, c]])
        }
    }
    
    # construct SCE
    gi <- do.call(rbind, lapply(cluster_ids, function(c)
        do.call(rbind, lapply(cats, function(cat) if (ndd[cat, c] != 0)
            data.frame(gene = is[[cat, c]], cluster_id = c, category = cat)))))
    gi <- gi[order(as.numeric(gsub("[a-z]", "", gi$gene))), ]
    gi$category <- factor(gi$category, levels = ddSingleCell:::cats)
    rownames(gi) <- NULL
    
    col_data <- do.call(rbind, lapply(cluster_ids, function(c) 
        do.call(rbind, lapply(sample_ids, function(s) 
            data.frame(
                row.names = 1,
                unlist(js[[s, c]]), 
                cluster_id = c, sample_id = s, 
                group = rep.int(c("A", "B"), n_cells[[s, c]]))))))
    col_data <- col_data[colnames(y), ]
    col_data$sample_id <- factor(paste(col_data$group, col_data$sample_id, sep = "."))
    
    sample_id <- levels(col_data$sample_id)
    group <- gsub("(A|B)[.].*", "\\1", sample_id)
    ei <- data.frame(sample_id, group)
    md <- list(
        experiment_info = ei,
        n_cells = table(col_data$sample_id),
        gene_info = gi, sim_genes = gs)
    
    SingleCellExperiment(
        assays = list(counts = y),
        colData = col_data, 
        metadata = md)
}
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    