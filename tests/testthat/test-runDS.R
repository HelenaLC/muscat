context("DS analysis using edgeR & limma")

# load packages
suppressPackageStartupMessages({
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

# randomly select 10 DE genes & multiply counts by 100 for group 2
gs <- sample(rownames(sce), 10)
g23 <- gids %in% c("g2", "g3")
assay(sce[gs, g23]) <- assay(sce[gs, g23]) * 10

# compute pseudo-bulks
pb <- aggregateData(sce, assay = "counts", fun = "sum")

# specify design & contrast matrices
ei <- metadata(sce)$experiment_info
design <- model.matrix(~ 0 + ei$group_id)
dimnames(design) <- list(ei$sample_id, levels(ei$group_id))
contrast <- limma::makeContrasts("g1-g2", "g1-g3", levels = design)

for (method in c("edgeR", "limma-trend", "limma-voom")) {
    test_that(paste("runDS", method, sep = "_"), {
        # test for cluster-wise differential expression
        res <- runDS(sce, pb, design, contrast, method = method, verbose = FALSE)
        # check that nb. of DE genes is 10 in ea. comparison & cluster
        n_de <- map_depth(res$table, 2, function(u) sum(u$p_adj < 1e-6))
        expect_true(all(unlist(n_de) == 10))
        # check that DE genes are correct
        de_gs <- map_depth(res$table, 2, function(u)
            u$gene[order(u$p_adj)][seq_len(10)])
        expect_true(all(unlist(map_depth(de_gs, 2, setequal, gs))))
    })
}
