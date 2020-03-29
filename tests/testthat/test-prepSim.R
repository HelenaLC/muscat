# generate toy dataset
seed <- as.numeric(format(Sys.time(), "%s"))
set.seed(seed); x <- .toySCE()

nk <- length(kids <- levels(x$cluster_id))
ns <- length(sids <- levels(x$sample_id))
ng <- length(gids <- levels(x$group_id))

test_that("prepSim()", {
    z <- x; class(z) <- "x"
    expect_error(prepSim(z))
    expect_error(prepSim(x, group_keep = "x"))
    y <- prepSim(x, 0, 0, 0, 0)
    expect_is(y, "SingleCellExperiment")
    expect_true(nrow(y) == nrow(x))
    expect_true(ncol(y) == sum(x$group_id == gids[1]))
    g <- sample(setdiff(gids, gids[1]), 1)
    y <- prepSim(x, 0, 0, 0, 0, group_keep = g)
    expect_true(ncol(y) == sum(x$group_id == g))
})
