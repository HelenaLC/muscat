suppressPackageStartupMessages({
    library(SingleCellExperiment)
})

x <- .toySCE()
kids <- levels(x$cluster_id)
sids <- levels(x$sample_id)

test_that(".filter_matrix()", {
    replicate(5, {
        x <- matrix(sample(seq_len(200)), 10, 20)
        y <- .filter_matrix(x, (n <- sample(200, 1)))
        expect_true(all(y >= n) | any(dim(y) == 1))
    })
})
test_that(".update_sce()", {
    replicate(5, {
        ks <- sample(kids, (nk <- 2))
        ss <- sample(sids, (ns <- 3))
        cs <- x$cluster_id %in% ks & x$sample_id %in% ss
        y <- .update_sce(x[, cs])
        expect_equal(ncol(y), sum(cs))
        expect_true(setequal(levels(y$cluster_id), ks))
        expect_true(setequal(levels(y$sample_id), ss))
    })
})
test_that(".filter_sce()", {
    replicate(5, {
        ks <- sample(kids, (nk <- 2))
        ss <- sample(sids, (ns <- 3))
        y <- .filter_sce(x, ks, ss)
        expect_true(setequal(levels(y$cluster_id), ks))
        expect_true(setequal(levels(y$sample_id), ss))
        cs <- x$cluster_id %in% ks & x$sample_id %in% ss
        expect_equal(ncol(y), sum(cs))
        expect_equal(colnames(y), colnames(x)[cs])
        ei <- metadata(y)$experiment_info
        expect_true(setequal(ei$sample_id, ss))
    })
})
test_that(".scale()", {
    replicate(5, {
        y <- .scale(x <- matrix(runif(200), (n <- 10), (m <- 20)))
        expect_true(all(apply(y, 1, min) == 0))
        expect_true(all(apply(y, 1, max) == 1))
        qs <- quantile(x[1, ], c(0.01, 0.99))
        z <- (x[1, ] - qs[1]) / diff(qs)
        z[z < 0] <- 0; z[z > 1] <- 1
        expect_true(all(y[1, ] == z))
        i <- sample(n, 1); j <- sample(m, 1)
        x[i, j] <- NA; z <- .scale(x)
        expect_true(sum(is.na(z)) == 1)
        os <- lapply(list(y, z), function(u) {
            rngs <- colRanges(u[-i, -j], na.rm = TRUE)
            apply(rngs, 2, order)
        })
        expect_identical(os[[1]], os[[2]])
    })
})
