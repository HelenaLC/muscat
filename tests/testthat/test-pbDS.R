context("DS analysis using pseudo-bulks")

# load packages
suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(SummarizedExperiment)
})

# generate toy dataset
seed <- as.numeric(format(Sys.time(), "%s"))
set.seed(seed)
sce <- toyData()

kids <- sce$cluster_id
sids <- sce$sample_id
gids <- sce$group_id

# randomly select n_de DE genes & multiply counts by 100 for group 2
n_de <- 10
de_gs <- sample(rownames(sce), n_de)
g23 <- gids %in% c("g2", "g3")
assay(sce[de_gs, g23]) <- assay(sce[de_gs, g23]) * 100

# pbDS() -----------------------------------------------------------------------
pb <- aggregateData(sce, assay = "counts", fun = "sum")
ei <- metadata(sce)$experiment_info
design <- model.matrix(~ 0 + ei$group_id)
dimnames(design) <- list(ei$sample_id, levels(ei$group_id))
contrast <- limma::makeContrasts("g2-g1", "g3-g1", levels = design)

for (method in c("edgeR", "limma-trend", "limma-voom")) {
    test_that(paste("pbDS", method, sep = "."), {
        # test for cluster-wise differential expression
        res <- pbDS(sce, pb, design, contrast, method = method, verbose = FALSE)
        
        expect_is(res, "list")
        expect_identical(length(res[[1]]), ncol(contrast))
        expect_identical(names(res[[1]]), colnames(contrast))
        expect_true(all(vapply(map(res[[1]], names), "==", 
            levels(kids), FUN.VALUE = logical(nlevels(kids)))))
        
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

# mmDS() -----------------------------------------------------------------------
# test_that("mmDS;method='dream'", {
#     res <- mmDS(sce, method = "dream", verbose = TRUE)
# 
#     expect_is(res, "list")
#     expect_identical(names(res), levels(kids))
# 
#     p_adj <- map(res, pull, "p_adj.loc")
#     n_de_res <- unlist(map(p_adj, function(u) sum(u < 1e-6)))
#     expect_true(all(n_de_res == n_de))
# 
#     os <- map(p_adj, order)
#     de_gs_res <- map(os, function(o) rownames(sce)[o][seq_len(n_de)])
#     expect_true(all(unlist(map(de_gs_res, setequal, de_gs))))
# })

# global p-value adjustment ----------------------------------------------------
test_that(".p_adj_global", {
    cs <- paste0("c", seq_len(5))
    ks <- paste0("k", seq_len(8))
    ns <- sample(1e3, length(ks))
    names(cs) <- cs
    names(ks) <- ks
    names(ns) <- ks
    df <- lapply(cs, function(c)
        lapply(ks, function(k)
            data.frame(
                p_val = rgamma(ns[k], 0.1),
                p_adj.loc = rgamma(ns[k], 0.1))))
    df_adj <- .p_adj_global(df)
    for (c in cs)
        expect_identical(
            p.adjust(bind_rows(df[[c]])$p_val),
            bind_rows(df_adj[[c]])$p_adj.glb)
})
