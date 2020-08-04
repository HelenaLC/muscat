# generate toy dataset
seed <- as.numeric(format(Sys.time(), "%s"))
set.seed(seed); x <- .toySCE()

nk <- length(kids <- levels(x$cluster_id))
ns <- length(sids <- levels(x$sample_id))
ng <- length(gids <- levels(x$group_id))

test_that("pbMDS()", {
    y <- x; class(y) <- "x"
    expect_error(pbMDS(y))
    y <- aggregateData(x)
    expect_is(p <- pbMDS(y), "ggplot")
    expect_true(nrow(p$data) == nk*ns)
    expect_true(all(table(p$data$cluster_id) == ns))
    expect_true(all(table(p$data$group_id) == ns*nk/ng))
    cs1 <- x$group_id != gids[1]   # remove group
    cs2 <- x$sample_id  != sids[1] # remove sample
    cs3 <- x$cluster_id != kids[1] # remove cluster
    cs4 <- cs1 & cs2 # remove sample in single group
    cs5 <- cs1 & cs3 # remove cluster in single group
    cs6 <- cs2 & cs3 # remove cluster-sample instance
    for (cs in paste0("cs", seq_len(6))) 
        expect_silent(pbMDS(aggregateData(x[, get(cs)])))
})
