# load packages
suppressMessages({
    library(dplyr)
    library(purrr)
    library(SingleCellExperiment)
})

# generate toy dataset
seed <- as.numeric(format(Sys.time(), "%s"))
set.seed(seed)
x <- .toySCE()

nk <- length(kids <- levels(x$cluster_id))
ns <- length(sids <- levels(x$sample_id))
ng <- length(gids <- levels(x$group_id))

# sample 'n_de' genes & multiply counts by 10 for 'g2/3'-cells
g23 <- x$group_id != "g1"
de_gs <- sample(rownames(x), (n_de <- 5))
assay(x[de_gs, g23]) <- assay(x[de_gs, g23]) * 10

# aggregate & run pseudobulk DS analysis
nc <- length(cs <- list(2, 3))
y <- aggregateData(x, assay = "counts", fun = "sum")
y <- pbDS(y, coef = cs, verbose = FALSE)
    
test_that("resDS()", {
    v <- list(col = list(nr = nrow(x)*nk, ng = nrow(x), nk = nk))
    v$row <- lapply(v$col, "*", nc)
    v$col$char_cols <- c("gene", "cluster_id")
    v$row$char_cols <- c(v$col$char_cols, "coef")
    for (bind in c("row", "col")) {
        z <- resDS(x, y, bind, frq = FALSE, cpm = FALSE)
        expect_is(z, "data.frame")
        expect_identical(nrow(z), v[[bind]]$nr)
        expect_true(all(table(z$gene) == v[[bind]]$nk))
        expect_true(all(table(z$cluster_id) == v[[bind]]$ng))
        is_char <- colnames(z) %in% v[[bind]]$char_cols
        expect_true(all(apply(z[, !is_char], 2, class) == "numeric"))
        expect_true(all(apply(z[,  is_char], 2, class) == "character"))
    }
})
test_that("resDS() - 'frq = TRUE'", {
    z <- resDS(x, y, frq = TRUE)
    u <- z[, grep("frq", colnames(z))]
    expect_true(ncol(u) == ns + ng)
    expect_true(all(u <= 1 & u >= 0 | is.na(u)))
    # remove single cluster-sample instance
    s <- sample(sids, 1); k <- sample(kids, 1)
    x_ <- x[, !(x$sample_id == s & x$cluster_id == k)]
    y_ <- aggregateData(x_, assay = "counts", fun = "sum")
    y_ <- pbDS(y_, coef = cs, verbose = FALSE)
    z <- resDS(x_, y_, frq = TRUE)
    u <- z[, grep("frq", colnames(z))]
    expect_true(ncol(u) == ns + ng)
    expect_true(all(u <= 1 & u >= 0 | is.na(u)))
    expect_true(all(z[z$cluster_id == k, paste0(s, ".frq")] == 0))
})
test_that("resDS() - 'cpm = TRUE'", {
    z <- resDS(x, y, cpm = TRUE)
    u <- z[, grep("cpm", colnames(z))]
    expect_true(ncol(u) == ns)
    expect_true(all(u %% 2 == 0 | is.na(u)))
})
