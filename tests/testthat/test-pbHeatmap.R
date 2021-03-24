# load packages
suppressMessages({
    library(purrr)
    library(scater)
    library(SingleCellExperiment)
})

# generate toy dataset
set.seed(as.numeric(format(Sys.time(), "%s")))
sce <- .toySCE()

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
ng <- length(gids <- levels(sce$group_id))

g3 <- sce$group_id == "g3"
g23 <- sce$group_id %in% c("g2", "g3")

# sample 'nde' genes & multiply counts by 10 for 'g2'- & 'g3'-cells
degs <- sample(rownames(sce), (nde <- 5))
assay(sce[degs, g23]) <- assay(sce[degs, g23]) * 10
sce <- logNormCounts(computeLibraryFactors(sce))

# run DS analysis using 'edgeR' on pseudobulks
pb <- aggregateData(sce, assay = "counts", fun = "sum")
res <- pbDS(pb, method = "edgeR", min_cells = 0, filter = "none", verbose = FALSE)

# compute pseudobulk mean of logcounts
pb <- aggregateData(sce, assay = "logcounts", fun = "mean")

test_that("pbHeatmap() - input arguments", {
    # get list of default parameters
    defs <- as.list(eval(formals(pbHeatmap)))
    defs$x <- sce; defs$y <- res
    expect_is(do.call(pbHeatmap, defs), "Heatmap")
    # check that invalid arguments throw error
    fail <- list(g = "x", k = "x", c = "x",
        sort_by = "x", sort_by = "gene",
        sort_by = c("logFC", "p_val"),
        lfc = "x", lfc = c(1, 2),
        assay = 1, assay = "x",
        decreasing = "x", decreasing = c(TRUE, FALSE),
        row_anno = "x", row_anno = c(TRUE, FALSE), 
        col_anno = "x", col_anno = c(TRUE, FALSE),
        normalize = "x", normalize = c(TRUE, FALSE))
    for (i in seq_along(fail)) {
        args <- defs; args[[names(fail)[i]]] <- fail[[i]]
        expect_error(do.call(pbHeatmap, args))
    }
})
    
test_that("pbHeatmap() - subset of clusters", {
    ks <- sample(kids, (nks <- 3))
    p <- pbHeatmap(sce, res, k = ks, 
        lfc = 0, fdr = Inf, sort_by = "none",
        top_n = (nds <- 10), normalize = FALSE)
    expect_is(p, "Heatmap")
    gs <- rownames(y <- p@matrix)[seq_len(nds)]
    expect_equal(dim(y), c(nks*nds, ns))
    expect_identical(rownames(y), rep(rownames(sce)[seq_len(nds)], nks))
    expect_identical(colnames(y), sids)
    for (i in seq_len(nks)) {
        k <- kids[match(ks, kids)][i]
        expect_equal(y[seq_len(nds)+nds*(i-1), ], assay(pb, k)[gs, ])
    }
})

test_that("pbHeatmap() - subset of genes", {
    gs <- sample(rownames(sce), (ngs <- 20))
    p <- pbHeatmap(sce, res, g = gs, 
        lfc = 0, fdr = Inf, sort_by = "none",
        top_n = Inf, normalize = FALSE)
    expect_is(p, "Heatmap")
    gs <- rownames(y <- p@matrix)[seq_len(ngs)]
    expect_equal(dim(y), c(nk*ngs, ns))
    expect_identical(rownames(y), rep(gs, nk))
    expect_identical(colnames(y), sids)
})
