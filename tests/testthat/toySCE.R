toySCE <- function() {
    ngs <- 300
    ncs <- 2e3
    
    gs <- paste0("gene", seq_len(ngs))
    cs <- paste0("cell", seq_len(ncs))
    
    counts <- sample(0:100, ngs * ncs, TRUE)
    counts <- matrix(counts,
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
        assay = list(counts = counts), colData = cd,
        metadata = list(experiment_info = .make_ei(cd)))
}
