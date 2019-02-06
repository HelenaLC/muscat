context("Aggregation of single-cell to pseudo-bulk data")

# load packages
suppressPackageStartupMessages({
    library(SummarizedExperiment)
})

# generate toy dataset
seed <- as.numeric(format(Sys.Date(), "%s"))
set.seed(seed)
sce <- toyData()

cids <- colData(sce)$cluster_id
sids <- colData(sce)$sample_id

# compute pseudobulks


test_that("aggregateData", {
  expect_error(aggregateData(sce, by = "xxx"))
  expect_error(aggregateData(sce, assay = "xxx"))
  expect_error(aggregateData(sce, fun = "xxx"))
  
  pb <- aggregateData(sce, assay = "counts", by = c("cluster_id", "sample_id"), fun = "sum")

  expect_is(pb, "SingleCellExperiment")
  expect_true(identical(levels(cids), assayNames(pb)))
  expect_true(identical(levels(sids), colnames(pb)))
  expect_true(identical(rownames(sce), rownames(pb)))
    
   
    expect_equal(length(pb), nlevels(cluster_ids))
    expect_true(all(sapply(pb, ncol) == nlevels(sample_ids)))
    expect_true(all(sapply(pb, nrow) == nrow(sce)))
    expect_equal(sum(sapply(pb, sum)), sum(assay(sce)))

    # random spot check
    c <- sample(levels(cluster_ids), 1)
    s <- sample(levels(sample_ids), 1)
    g <- sample(rownames(sce), 1)
    cells_keep <- sample_ids == s & cluster_ids == c
    expect_equal(pb[[c]][g, s], sum(assay(sce)[g, cells_keep]))
})


