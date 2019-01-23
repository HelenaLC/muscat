context("Expression frequencies by sample & group")

# generate toy dataset
seed <- as.numeric(format(Sys.time(), "%s"))
sce <- toyData(seed = seed)

# randomly split into 2 groups
ei <- data.frame(sample_id = levels(colData(sce)$sample_id))
g2 <- sample(ei$sample_id, round(nlevels(ei$sample_id) / 2))
ei$group_id <- "A"
m <- match(g2, ei$sample_id)
ei$group_id[m] <- "B"

metadata(sce)$experiment_info <- ei
m <- match(colData(sce)$sample_id, ei$sample_id)
colData(sce)$group_id <- factor(ei$group_id[m])

# calculate expr. freqs.
x <- calcExprFreqs(sce)

test_that("Frequencies lie in [0, 1] w/o NAs", {
    expect_true(all(!sapply(x, function(u) any(is.na(u)))))
    r <- sapply(x, range)
    expect_true(all(r[1, ] >= 0))
    expect_true(all(r[2, ] <= 0))
})

test_that("Output is list w/ names corresponding to cluster IDs", {
    expect_is(x, "list")
    cluster_ids <- levels(colData(sce)$cluster_id)
    expect_true(all(names(x) %in% cluster_ids))
})

test_that("Row names = genes & column names = group/sample IDs", {
    expect_true(all(sapply(x, function(u)
        all.equal(rownames(u), rownames(sce)))))
    group_ids <- levels(colData(sce)$group_id)
    sample_ids <- levels(colData(sce)$sample_id)
    expect_true(all(sapply(x, function(u)
        all(colnames(u) %in% c(group_ids, sample_ids)))))
})

test_that("Spot check", {
    # sample cluster
    cluster_ids <- colData(sce)$cluster_id
    k <- sample(levels(cluster_ids), 1)
    ki <- cluster_ids == k
    # sample sample & group
    s <- sample(ei$sample_id, 1)
    g <- sample(ei$group_id, 1)
    si <- colData(sce)$sample_id == s & ki
    gi <- colData(sce)$group_id == g & ki
    # sample gene & check frequencies vs. truth
    gene <- sample(rownames(sce), 1)
    expect_identical(sum(counts(sce)[gene, si] > 0) / sum(si), x[[k]][gene, s])
    expect_identical(sum(counts(sce)[gene, gi] > 0) / sum(gi), x[[k]][gene, g])
})