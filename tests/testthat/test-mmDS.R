context("DS analysis using mixed-models")

# load packages
suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(SingleCellExperiment)
})

# generate toy dataset
seed <- as.numeric(format(Sys.time(), "%s"))
set.seed(seed)
sce <- toySCE()

kids <- sce$cluster_id
sids <- sce$sample_id
gids <- sce$group_id

# randomly select n_de DE genes & multiply counts by 100 for group 2
n_de <- 10
de_gs <- sample(rownames(sce), n_de)
g3 <- gids == "g3"
assay(sce[de_gs, g3]) <- assay(sce[de_gs, g3]) * 100

for (vst in c("DESeq2", "sctransform")) {
    test_that(paste0("mmDS-vst-", vst), {
        res <- mmDS(sce, method = "vst", vst = vst, verbose = TRUE)
        
        expect_is(res, "list")
        expect_identical(names(res), levels(kids))
        
        p_adj <- map(res, pull, "p_adj.loc")
        n_de_res <- unlist(map(p_adj, function(u) sum(u < 1e-12)))
        expect_true(all(n_de_res == n_de))
        
        os <- map(p_adj, order)
        de_gs_res <- map(os, function(o) rownames(sce)[o][seq_len(n_de)])
        expect_true(all(unlist(map(de_gs_res, setequal, de_gs))))
    })
}



