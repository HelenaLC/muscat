# generate toy dataset
seed <- as.numeric(format(Sys.Date(), "%s"))
data <- toyData(seed = seed)

cluster_ids <- colData(data)$cluster_id
sample_ids <- colData(data)$sample_id

# compute pseudobulks
pb <- aggregateData(data, data = "counts", fun = "sum")

test_that("aggregateData", {
    expect_error(aggregateData(data, data = "missing_assay", fun = "sum"))
    expect_error(aggregateData(data, data = "counts", fun = "missing_function"))
    
    expect_is(pb, "list")
    expect_equal(length(pb), nlevels(cluster_ids))
    expect_true(all(sapply(pb, ncol) == nlevels(sample_ids)))
    expect_true(all(sapply(pb, nrow) == nrow(data)))
    expect_equal(sum(sapply(pb, sum)), sum(assay(data)))

    # random spot check
    c <- sample(levels(cluster_ids), 1)
    s <- sample(levels(sample_ids), 1)
    g <- sample(rownames(data), 1)
    cells_keep <- sample_ids == s & cluster_ids == c
    expect_equal(pb[[c]][g, s], sum(assay(data)[g, cells_keep]))
})


