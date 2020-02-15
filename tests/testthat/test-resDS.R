context("DS analysis results reformatting")

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
nc <- length(cs <- c(2, 3))
y <- aggregateData(x, assay = "counts", fun = "sum")
y <- pbDS(y, coef = cs, verbose = FALSE)
    
test_that("resDS()", {
    v <- list(col = list(nr = nrow(x)*nk, ng = nk, nk = nrow(x)))
    v$row <- lapply(v$col, "*", nc)
    v$col$char_cols <- c("gene", "cluster_id")
    v$row$char_cols <- c(v$col$char_cols, "coef")
    for (bind in c("row", "col")) {
        z <- resDS(x, y, bind, frq = FALSE, cpm = FALSE)
        expect_is(z, "data.frame")
        expect_identical(nrow(z), v[[bind]]$nr)
        expect_true(all(table(z$gene) == v[[bind]]$ng))
        expect_true(all(table(z$cluster_id) == v[[bind]]$nk))
        is_char <- colnames(z) %in% v[[bind]]$char_cols
        expect_true(all(apply(z[, !is_char], 2, class) == "numeric"))
        expect_true(all(apply(z[,  is_char], 2, class) == "character"))
    }
    # TODO: with expression frequencies
    
    # TODO: with cluster-sample wise CPMs
})
