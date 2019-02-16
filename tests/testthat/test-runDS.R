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

kids <- colData(sce)$cluster_id
sids <- colData(sce)$sample_id
gids <- colData(sce)$group_id

# randomly select 10 DE genes & multiply counts by 1000 for group 2
gs <- sample(rownames(sce), 10)
g2 <- gids == "g2"
assay(sce[gs, g2]) <- assay(sce[gs, g2]) * 10

# compute pseudo-bulks
pb <- aggregateData(sce, assay = "counts", fun = "sum")

# specify design & contrast matrices
ei <- metadata(sce)$experiment_info
design <- model.matrix(~ 0 + ei$group)
dimnames(design) <- list(ei$sample_id, levels(ei$group))
contrast <- limma::makeContrasts("g1-g2", levels = design)

for (method in c("edgeR", "limma")) {
    test_that(paste("runDS", method, sep = "_"), {
        # test for cluster-wise differential expression
        res <- runDS(sce, pb, design, contrast, method = method, verbose = FALSE)
        res <- res$table[[1]]   
        # check that nb. of DE genes is 10 in ea. cluster
        n_de <- vapply(res, function(u) sum(u$p_adj < 1e-3), numeric(1))
        expect_true(all(n_de == 10))
        # check that DE genes are correct
        de_gs <- lapply(res, function(u)
            u$gene[order(u$p_adj)][seq_len(10)])
        expect_true(all(vapply(levels(kids), function(k) 
            setequal(gs, de_gs[[k]]), logical(1))))
    })
}
