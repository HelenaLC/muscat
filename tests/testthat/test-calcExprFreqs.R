context("Expression frequencies by sample & group")

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
gids <- colData(sce)$group_id

# put in 50% of random 0s
n <- length(assay(sce))
i <- sample(assay(sce), round(n * 0.5))
assay(sce)[i] <- 0

# calculate expr. freqs.
sce <- prepData(sce, "cluster_id", "sample_id", "group_id")
x <- calcExprFreqs(sce, assay = "counts", th = 0)

test_that("Output is correctly structured SE", {
    expect_is(x, "SummarizedExperiment")
    expect_identical(assayNames(x), levels(kids))
    
    expect_identical(nrow(x), nrow(sce))
    expect_identical(ncol(x), nlevels(sids) + nlevels(gids))
    
    expect_identical(rownames(x), rownames(sce))
    expect_identical(colnames(x), c(levels(sids), levels(gids)))
})

test_that("Frequencies lie in [0, 1] w/o NAs", {
    expect_true(all(!vapply(assays(x), function(u) any(is.na(u)), logical(1))))
    r <- vapply(assays(x), range, numeric(2))
    expect_true(all(r[1, ] >= 0))
    expect_true(all(r[2, ] <= 1))
})

test_that("Spot check", {
    # sample cluster
    k <- sample(levels(kids), 1)
    ki <- kids == k
    # sample sample & group
    s <- sample(levels(sids), 1)
    g <- sample(levels(gids), 1)
    si <- sids == s & ki
    gi <- gids == g & ki
    # sample gene & check frequencies vs. truth
    gene <- sample(rownames(sce), 1)
    expect_identical(sum(counts(sce)[gene, si] > 0) / sum(si), assays(x)[[k]][gene, s])
    expect_identical(sum(counts(sce)[gene, gi] > 0) / sum(gi), assays(x)[[k]][gene, g])
})

