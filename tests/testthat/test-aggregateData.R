context("Aggregation of single-cell to pseudo-bulk data")

# load packages
suppressPackageStartupMessages({
    library(SummarizedExperiment)
})

# generate toy dataset
seed <- as.numeric(format(Sys.Date(), "%s"))
set.seed(seed)
sce <- toyData()

kids <- colData(sce)$cluster_id
sids <- colData(sce)$sample_id

pb <- aggregateData(sce, assay = "counts", by = c("cluster_id", "sample_id"), fun = "sum")

test_that("aggregateData", {
    expect_error(aggregateData(sce, by = "xxx"))
    expect_error(aggregateData(sce, assay = "xxx"))
    expect_error(aggregateData(sce, fun = "xxx"))
    
    expect_is(pb, "SingleCellExperiment")
    expect_identical(assayNames(pb), levels(kids))
    
    expect_identical(nrow(pb), nrow(sce))
    expect_identical(ncol(pb), nlevels(sids))
    
    expect_identical(rownames(pb), rownames(sce))
    expect_identical(colnames(pb), levels(sids))
    
    # random spot check
    k <- sample(levels(kids), 1)
    s <- sample(levels(sids), 1)
    g <- sample(rownames(sce), 1)
    i <- sids == s & kids == k
    expect_equal(assays(pb)[[k]][g, s], sum(assay(sce)[g, i]))
})


