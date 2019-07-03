context("Aggregation of single-cell to pseudo-bulk data")

# load packages
suppressMessages(library(SummarizedExperiment))

# generate toy dataset
seed <- as.numeric(format(Sys.Date(), "%s"))
set.seed(seed)
sce <- toySCE()
kids <- sce$cluster_id
sids <- sce$sample_id

test_that("aggregation across 2 factors", {
    for (fun in c("sum", "mean", "median")) {
        pb <- aggregateData(sce, by = c("cluster_id", "sample_id"), fun = fun)
        expect_error(aggregateData(sce, assay = "x"))
        expect_error(aggregateData(sce, fun = "x"))
        expect_error(aggregateData(sce, by = "x"))
        
        expect_is(pb, "SingleCellExperiment")
        expect_identical(assayNames(pb), levels(kids))
        
        expect_identical(nrow(pb), nrow(sce))
        expect_identical(ncol(pb), nlevels(sids))
        
        expect_identical(rownames(pb), rownames(sce))
        expect_identical(colnames(pb), levels(sids))
        
        expect_identical(as.numeric(table(sids)), pb$n_cells)
        
        # random spot check
        k <- sample(levels(kids), 1)
        s <- sample(levels(sids), 1)
        g <- sample(rownames(sce), 1)
        i <- sids == s & kids == k
        expect_equal(assays(pb)[[k]][g, s], get(fun)(assay(sce)[g, i]))
    }
})

test_that("aggregation across 1 factor", {
    for (fun in c("sum", "mean", "median")) {
        pb <- aggregateData(sce, by = "cluster_id", fun = fun)
        expect_is(pb, "SingleCellExperiment")
        
        expect_identical(nrow(pb), nrow(sce))
        expect_identical(ncol(pb), nlevels(kids))
        
        expect_identical(rownames(pb), rownames(sce))
        expect_identical(colnames(pb), levels(kids))
        
        expect_identical(as.numeric(table(kids)), pb$n_cells)
        
        # random spot check
        k <- sample(levels(kids), 1)
        g <- sample(rownames(sce), 1)
        i <- kids == k
        expect_equal(assay(pb)[g, k], get(fun)(assay(sce)[g, i]))
    }
})
