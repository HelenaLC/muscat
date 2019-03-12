context("DS analysis using pseudo-bulks")

# load packages
suppressPackageStartupMessages({
    library(purrr)
    library(SummarizedExperiment)
})

# generate toy dataset
seed <- as.numeric(format(Sys.time(), "%s"))
set.seed(seed)
sce <- toyData()
sce <- prepData(sce, "cluster_id", "sample_id", "group_id")

kids <- sce$cluster_id
sids <- sce$sample_id
gids <- sce$group_id

# randomly select n_de DE genes & multiply counts by 100 for group 2
n_de <- 10
de_gs <- sample(rownames(sce), n_de)
g23 <- gids %in% c("g2", "g3")
assay(sce[de_gs, g23]) <- assay(sce[de_gs, g23]) * 100

# compute pseudo-bulks
pb <- aggregateData(sce, assay = "counts", fun = "sum")

# specify design & contrast matrices
ei <- metadata(sce)$experiment_info
design <- model.matrix(~ 0 + ei$group_id)
dimnames(design) <- list(ei$sample_id, levels(ei$group_id))
contrast <- limma::makeContrasts("g1-g2", "g1-g3", levels = design)

for (method in c("edgeR", "limma-trend", "limma-voom")) {
    test_that(paste("pbDS", method, sep = "_"), {
        # test for cluster-wise differential expression
        res <- pbDS(sce, pb, design, contrast, method = method, verbose = FALSE)
        # check that nb. of DE genes is n_de in ea. comparison & cluster
        p_adj <- map_depth(res$table, 2, "p_adj.loc")
        n_de_res <- unlist(map_depth(p_adj, 2, function(u) sum(u < 1e-6)))
        expect_true(all(n_de_res == n_de))
        # check that DE genes are correct
        de_gs_res <- map_depth(res$table, 2, function(u)
            u$gene[order(u$p_adj.loc)][seq_len(n_de)])
        expect_true(all(unlist(map_depth(de_gs_res, 2, setequal, de_gs))))
    })
}
