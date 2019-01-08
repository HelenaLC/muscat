# generate toy dataset
seed <- as.numeric(format(Sys.time(), "%s"))
data <- toyData(seed = seed)

cluster_ids <- colData(data)$cluster_id
sample_ids <- colData(data)$sample_id

# compute pseudobulks
pb <- aggregateData(data, data = "counts", fun = "sum")

# randomly select 10 DE genes & multiply counts by 100 in half of samples
g2 <- sample(levels(sample_ids), round(nlevels(sample_ids) / 2))
idx <- replicate(nlevels(cluster_ids), sample(nrow(data), 10))
colnames(idx) <- levels(cluster_ids)
for (k in levels(cluster_ids))
    pb[[k]][idx[, k], g2] <- 100 * pb[[k]][idx[, k], g2]

# specify design & contrast matrices
ei <- data.frame(sample_id = levels(sample_ids))
ei$group <- factor(c("A", "B")[as.numeric(ei$sample_id %in% g2) + 1])
design <- model.matrix(~ 0 + ei$group)
dimnames(design) <- list(ei$sample_id, levels(ei$group))
contrast <- makeContrasts("B-A", levels = design)

# test for cluster-wise differential expression
res <- run_edgeR(data, pb, design, contrast, verbose = FALSE)[[1]]
res_by_cluster <- split(res, res$cluster_id)

test_that("run_edgeR", {
    # check that nb. of DE genes is 10 in ea. cluster
    n_de <- sapply(res_by_cluster, function(x) sum(x$FDR < 1e-3))
    expect_true(all(n_de == 10))
    
    # check that DE genes are correct
    de_gs <- sapply(res_by_cluster, function(x) filter(x, FDR < 1e-3)$gene)
    expect_true(all(sapply(levels(cluster_ids), function(k)
        all(rownames(data)[idx[, k]] %in% de_gs[, k]))))
})