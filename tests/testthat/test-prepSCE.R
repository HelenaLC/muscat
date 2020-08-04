# generate toy dataset
seed <- as.numeric(format(Sys.time(), "%s"))
set.seed(seed); x <- .toySCE()

nk <- length(kids <- levels(x$cluster_id))
ns <- length(sids <- levels(x$sample_id))
ng <- length(gids <- levels(x$group_id))

test_that("prepSCE()", {
    y <- x; class(y) <- "x"
    expect_error(prepSCE(y))
    expect_is(y <- prepSCE(x), "SingleCellExperiment")
    expect_identical(dim(y), dim(x))
    expect_identical(dimnames(y), dimnames(x))
    ids <- formals("prepSCE")
    ids <- ids[grep("[a-z]id", names(ids))]
    ids <- unname(unlist(ids))
    expect_identical(colnames(colData(y)), ids)
    x$foo <- sample(ncol(x))
    y <- prepSCE(x, drop = TRUE)
    expect_true(is.null(y$foo))
    y <- prepSCE(x, drop = FALSE)
    expect_identical(y$foo, x$foo)
})
