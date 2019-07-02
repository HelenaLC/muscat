.toy_sce <- function() {
    n_cells <- 2e3
    n_genes <- 500

    cs <- paste0("cell", seq_len(n_cells))
    gs <- paste0("gene", seq_len(n_genes))
    y <- replicate(n_cells, rnbinom(n_genes, 1, 0.1))
    dimnames(y) <- list(gs, cs)
    
    kids <- paste0("k", seq_len(3))
    sids <- paste0("s", seq_len(2))
    gids <- paste0("g", seq_len(3))
    cd <- data.frame(
        cluster_id = sample(kids, n_cells, replace = TRUE),
        sample_id = sample(sids, n_cells, replace = TRUE),
        group_id = sample(gids, n_cells, replace = TRUE))
    cd$sample_id <- factor(with(cd, 
        paste(sample_id, group_id, sep = ".")))
    
    sce <- SingleCellExperiment(
        assays = list(counts = y),
        colData = cd)
    
    ei <- .make_ei(sce)
    metadata(sce)$experiment_info <- ei
    
    return(sce)
}