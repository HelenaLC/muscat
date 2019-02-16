context("DS analysis using MAST")

# load packages
suppressPackageStartupMessages({
    library(limma)
    library(SummarizedExperiment)
})

# generate toy dataset
seed <- as.numeric(format(Sys.time(), "%s"))
set.seed(seed)
sce <- toyData()

# randomly select 10 DE genes & multiply counts by 1000 for groups 2 & 3
cs_by_g <- .split_cells(sce, "group_id")
g23 <- unlist(cs_by_g[c("g2", "g3")])
gs <- sample(rownames(sce), 10)
assay(sce[gs, g23]) <- assay(sce[gs, g23]) * 1e3

# compute CPM
cpm <- edgeR::cpm(assay(sce))
assays(sce)$logcpm <- log2(cpm + 1)

# run MAST
ei <- metadata(sce)$experiment_info
contrast <- makeContrasts("g1-g2", "g1-g3", levels = levels(ei$group_id))
cs <- colnames(contrast)

res <- runMAST(sce, ~ 0 + group_id, contrast, assay = "logcpm")
tbl <- res$table

# ------------------------------------------------------------------------------

test_that("# & identity of DE genes is correct.", {
    # check that nb. of DE genes is 10 in ea. cluster
    n_de <- map_depth(tbl, 2, function(u)
        sum(u$p_adj < 1e-3))
    expect_true(all(unlist(n_de) == 10))
    # check that DE genes are correct
    de_gs <- modify_depth(tbl, 2, function(u)
        u$gene[order(u$p_adj)][seq_len(10)])
    expect_true(all(unlist(map_depth(de_gs, 2, setequal, gs))))
})









