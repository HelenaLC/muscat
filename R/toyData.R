toyData <- function(seed) {
    set.seed(seed)
    
    # specify ranges to sample from
    n_groups   <- 2:4
    n_samples  <- 2:5
    n_clusters <- 5:10
    n_cells    <- 100:500
    n_genes    <- 50:500
    counts     <- 0:1000
    
    # sample nb. of genes
    n_genes <- sample(n_genes, 1)
    genes <- paste0("gene", seq_len(n_genes))
    
    # sample nb. of groups
    n_groups <- sample(n_groups, 1)
    groups <- paste0("group", seq_len(n_groups))

    # sample nb. of samples per group
    n_samples <- sample(n_samples, n_groups)
    samples <- paste0("sample", seq_len(sum(n_samples)))
    
    # sample nb. of clusters
    n_clusters <- sample(n_clusters, 1)
    clusters <- paste0("cluster", seq_len(n_clusters))
    
    # sample nb. of cells per cluster-sample
    n_cells <- sample(n_cells, n_clusters * sum(n_samples), replace = TRUE)
    n_cells <- matrix(n_cells, n_clusters, sum(n_samples), dimnames = list(clusters, samples))
    cells <- paste0("cell", seq_len(sum(n_cells)))
        
    # sample sample and cluster IDs
    col_data <- matrix(NA, nrow = sum(n_cells), ncol = 2,
        dimnames = list(cells, c("sample_id", "cluster_id")))
    cells_tmp <- cells
    for (i in clusters) {
        for (j in samples) {
        n <- n_cells[i, j]
        idx <- sample(cells_tmp, n)
        col_data[idx, "cluster_id"] <- i
        cells_tmp <- setdiff(cells_tmp, idx)
        }
    }
    cells_tmp <- cells
    for (i in samples) {
        n <- sum(n_cells[, i])
        idx <- sample(cells_tmp, n)
        col_data[idx, "sample_id"] <- i
        cells_tmp <- setdiff(cells_tmp, idx)
    }
        
    # generate matrix of random integers
    counts <- sample(counts, sum(n_cells), replace = TRUE)
    counts <- matrix(counts, n_genes, sum(n_cells),
        dimnames = list(genes, cells))
    
    # construct SingleCellExperiment
    SingleCellExperiment(
        assays = list(counts = counts),
        colData = data.frame(col_data))
}
