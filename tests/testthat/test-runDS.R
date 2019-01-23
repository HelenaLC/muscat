context("DS analysis using edgeR & limma")

# generate toy dataset
seed <- as.numeric(format(Sys.time(), "%s"))
sce <- toyData(seed = seed)

cluster_ids <- colData(sce)$cluster_id
sample_ids <- colData(sce)$sample_id

# compute pseudobulks
pb <- aggregateData(sce, data = "counts", fun = "sum")

# randomly select 10 DE genes & multiply counts by 100 for half the samples
g2 <- sample(levels(sample_ids), round(nlevels(sample_ids) / 2))
idx <- replicate(nlevels(cluster_ids), sample(nrow(sce), 10))
colnames(idx) <- levels(cluster_ids)
for (k in levels(cluster_ids))
    pb[[k]][idx[, k], g2] <- 10 * pb[[k]][idx[, k], g2]

# specify design & contrast matrices
ei <- data.frame(sample_id = levels(sample_ids))
ei$group <- factor(c("A", "B")[as.numeric(ei$sample_id %in% g2) + 1])
design <- model.matrix(~ 0 + ei$group)
dimnames(design) <- list(ei$sample_id, levels(ei$group))
contrast <- limma::makeContrasts("B-A", levels = design)

for (method in c("edgeR", "limma")) {
    test_that(paste("runDS", method, sep = "_"), {
        # test for cluster-wise differential expression
        res <- runDS(sce, pb, design, contrast, method = method, verbose = FALSE)
        res <- res$table$`B-A`   
        # check that nb. of DE genes is 10 in ea. cluster
        i <- ifelse(method == "edgeR", "FDR", "adj.P.Val")
        n_de <- sapply(res, function(x) sum(x[[i]] < 1e-6))
        expect_true(all(n_de == 10))
        # check that DE genes are correct
        de_gs <- lapply(res, function(x) 
            x$gene[order(x[[i]])][seq_len(10)])
        expect_true(all(sapply(levels(cluster_ids), function(k)
            all(rownames(sce)[idx[, k]] %in% de_gs[[k]]))))
    })
}
