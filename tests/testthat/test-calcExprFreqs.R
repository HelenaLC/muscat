# load packages
suppressMessages({
    library(SingleCellExperiment)
})

# generate toy dataset
seed <- as.numeric(format(Sys.Date(), "%s"))
set.seed(seed)
sce <- .toySCE()

# put in 50% random 0s
n <- length(assay(sce))
i <- sample(n, round(n * 0.5))
assay(sce)[i] <- 0

# store nb. / IDs of clusters, samples, groups
nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
ng <- length(gids <- levels(sce$group_id))

# calculate expression frequencies
x <- calcExprFreqs(sce, assay = "counts", th = 0)

test_that("Output is correctly structured SCE", {
    expect_is(x, "SingleCellExperiment")
    expect_identical(assayNames(x), kids)
    
    expect_identical(nrow(x), nrow(sce))
    expect_identical(ncol(x), ns + ng)

    expect_identical(rownames(x), rownames(sce))
    expect_identical(colnames(x), c(sids, gids))
})

test_that("Frequencies lie in [0, 1] w/o any NAs", {
    v <- unlist(assays(x))
    expect_true(all(v >= 0))
    expect_true(all(v <= 1))
    expect_true(!any(is.na(v)))
})

test_that("10x random spot checks", {
    replicate(10, {
        # sample cluster, sample & group
        k <- sample(kids, 1)
        s <- sample(sids, 1)
        g <- sample(gids, 1)
        ki <- sce$cluster_id == k
        si <- sce$sample_id == s & ki
        gi <- sce$group_id == g & ki
        # sample gene & check frequencies vs. truth
        i <- sample(rownames(sce), 1)
        expect_identical(mean(counts(sce)[i, si] > 0), assays(x)[[k]][i, s])
        expect_identical(mean(counts(sce)[i, gi] > 0), assays(x)[[k]][i, g])
    })
})
