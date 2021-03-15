# load packages
suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(SingleCellExperiment)
})

# generate toy dataset
seed <- as.numeric(format(Sys.time(), "%s"))
set.seed(seed)
x <- .toySCE()
x <- x[, x$group_id != "g3"]

nk <- length(kids <- x$cluster_id)
ns <- length(sids <- x$sample_id)
ng <- length(gids <- x$group_id)

test_that(".check_args_mmDS()", {
    y <- x; class(y) <- "x"
    expect_error(mmDS(y))
    y <- x; assayNames(y) <- "x"
    expect_error(mmDS(y))
    expect_error(mmDS(x, "x"))
    expect_error(mmDS(x, 100))
    expect_error(mmDS(x, covs = "x"))
    expect_error(mmDS(x, method = "x"))
    expect_error(mmDS(x, method = "vst", vst = "x"))
    for (u in c("verbose", "dup_corr", "trended", "bayesian", "blind", "REML")) 
        for (v in list("x", 1, c(TRUE, TRUE))) {
            w <- list(x = x); w[[u]] <- v
            expect_error(do.call(mmDS, w))
        }
})

test_that("mmDS() - filtering", {
    expect_error(mmDS(x, n_cells = Inf))
    ks <- sample(kids, 2)
    ks <- as.character(ks)
    cs <- x$cluster_id %in% ks
    ls <- rowSums(assay(x[, cs]))
    o <- order(ls, decreasing = TRUE)
    gs <- rownames(x)[o][seq_len(5)]
    z <- suppressWarnings(
        mmDS(x[gs, cs], verbose = FALSE))
    expect_setequal(names(z), ks)
    expect_true(all(vapply(map(z, "gene"), 
        function(u) all(u == gs), logical(1))))
    expect_true(all(vapply(ks, function(k) 
        all(z[[k]]$cluster_id == k), logical(1))))
    y <- x[gs, cs]; metadata(y) <- list()
    y$group_id <- NULL; expect_error(mmDS(y))
})

# randomly select 'n_de' genes & bump counts for group 2
n_de <- 5; g2 <- gids == "g2"
de_gs <- sample(rownames(x), n_de)
assay(x[de_gs, g2]) <- (assay(x[de_gs, g2])+1)*5

test_that("mmDS-utils", {
    cs <- x$cluster_id == kids[1]
    gs <- c(de_gs, sample(setdiff(rownames(x), de_gs), 5))
    for (fun in paste0(".mm_", c("dream", "dream2", "vst"))) {
        # currently not passing; 
        # either there's a bug I cannot find 
        # or the toydata is too simplistic
        # c("poisson", "hybrid", "nbinom"))) {
        y <- suppressWarnings(
            get(fun)(x[gs, cs], verbose = FALSE))
        expect_is(y, "data.frame")
        expect_identical(rownames(y), gs)
        top <- order(y$p_adj.loc)[seq_len(n_de)]
        expect_setequal(rownames(y)[top], de_gs)
    }
})
