context("Aggregation of single-cell to pseudobulk data")

# load packages
suppressMessages({
    library(SingleCellExperiment)
})

# generate toy dataset
seed <- as.numeric(format(Sys.Date(), "%s"))
set.seed(seed)
sce <- .toySCE()
nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
ng <- length(gids <- levels(sce$group_id))

test_that("aggregation across 2 factors", {
    for (fun in c("sum", "mean", "median")) {
        pb <- aggregateData(sce, by = c("cluster_id", "sample_id"), fun = fun)
        expect_error(aggregateData(sce, assay = "x"))
        expect_error(aggregateData(sce, fun = "x"))
        expect_error(aggregateData(sce, by = "x"))
        
        pars <- metadata(pb)$agg_pars
        expect_identical(pars$by, c("cluster_id", "sample_id"))
        expect_identical(pars$assay, "counts")
        expect_identical(pars$fun, fun)
        expect_false(pars$scale)
        
        expect_is(pb, "SingleCellExperiment")
        expect_identical(assayNames(pb), kids)
        
        expect_identical(nrow(pb), nrow(sce))
        expect_identical(ncol(pb), ns)
        
        expect_identical(rownames(pb), rownames(sce))
        expect_identical(colnames(pb), sids)
        
        expect_equivalent(metadata(pb)$n_cells,
            table(sce$cluster_id, sce$sample_id))
        
        # 10x random spot check
        replicate(10, {
            k <- sample(kids, 1)
            s <- sample(sids, 1)
            g <- sample(rownames(sce), 1)
            i <- sce$sample_id == s & sce$cluster_id == k
            expect_equal(assays(pb)[[k]][g, s], get(fun)(assay(sce)[g, i]))
        })
    }
})

test_that("aggregation across 1 factor", {
    for (fun in c("sum", "mean", "median")) {
        pb <- aggregateData(sce, by = "cluster_id", fun = fun)
        expect_is(pb, "SingleCellExperiment")
        expect_identical(nrow(pb), nrow(sce))
        expect_identical(ncol(pb), nk)
        expect_identical(rownames(pb), rownames(sce))
        expect_identical(colnames(pb), kids)
        expect_equivalent(table(sce$cluster_id), metadata(pb)$n_cells)
        # random spot check
        k <- sample(kids, 1)
        g <- sample(rownames(sce), 1)
        i <- sce$cluster_id == k
        expect_equal(assay(pb)[g, k], get(fun)(assay(sce)[g, i]))
    }
})

test_that("pbFlatten()", {
    x <- aggregateData(sce, by = "cluster_id")
    expect_error(pbFlatten(x))
    x <- aggregateData(sce)
    cd <- colData(y <- pbFlatten(x))
    expect_is(y, "SingleCellExperiment")
    expect_true(length(assays(y)) == 2)
    expect_identical(rownames(y), rownames(x))
    expect_true(ncol(y) == nk*ncol(x))
    a <- do.call(cbind, as.list(assays(x)))
    expect_equivalent(assay(y), a)
    expect_true(all(table(cd$sample_id) == nk))
    expect_true(all(table(cd$cluster_id) == ns))
    expect_true(all(table(cd$group_id) == ns/ng*nk))
    expect_equivalent(y$n_cells, c(table(sce$sample_id, sce$cluster_id)))
    # without normalization
    y <- pbFlatten(x, normalize = FALSE)
    expect_true(assayNames(y) == "counts")
})
