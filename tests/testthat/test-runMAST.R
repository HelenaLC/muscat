context("DS analysis using MAST")

# load packages
suppressPackageStartupMessages({
    library(limma)
    library(SummarizedExperiment)
})

# generate toy dataset
seed <- as.numeric(format(Sys.time(), "%s"))
set.seed(seed)
sce <- toyData()
gs <- rownames(sce)

# split into 2 groups
sample_ids <- colData(sce)$sample_id
ei <- data.frame(sample_id = levels(sample_ids))
g2 <- sample(nrow(ei), round(nrow(ei) / 2))
ei$group_id <- "A"
ei$group_id[g2] <- "B"
ei$group_id <- factor(ei$group_id)
n_cells <- table(sample_ids)
colData(sce)$group_id <- rep(ei$group_id, n_cells)

# impute 10% DE genes
de_gs <- sample(gs, round(nrow(sce) / 10))
g2 <- colData(sce)$group_id == "B"
assay(sce)[de_gs, g2] <- assay(sce)[de_gs, g2] * 100

# compute CPM
cpm <- cpm(assay(sce))
assays(sce)$logcpm <- log2(cpm + 1)

# run MAST
contrast <- makeContrasts("B-A", levels = levels(ei$group_id))
res <- runMAST(sce, ~ 0 + group_id, contrast, assay = "logcpm")
tbl <- res$table$`B-A`

# ------------------------------------------------------------------------------

test_that("# & identity of genes is correct.", {
    # get significant hits
    p.adj <- lapply(tbl, "[[", "p.adj")
    de <- lapply(p.adj, "<", 1e-3)
    
    # check that nb. of DE genes is right
    n_de <- vapply(de, sum, numeric(1))
    expect_true(all(n_de == length(de_gs)))

    # check that DE genes are right
    expect_true(all(vapply(de, function(u) 
        identical(sort(gs[u]), sort(de_gs)), 
        logical(1))))
})









