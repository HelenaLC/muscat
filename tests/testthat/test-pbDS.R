context("DS analysis via pseudobulks")

# load packages
suppressMessages({
    library(dplyr)
    library(purrr)
    library(SingleCellExperiment)
})

# generate toy dataset
seed <- as.numeric(format(Sys.time(), "%s"))
set.seed(seed)
sce <- .toySCE()

kids <- sce$cluster_id
sids <- sce$sample_id
gids <- sce$group_id

# sample 'n_de' DS genes & multiply counts by 100 for group IDs "g2" and "g3"
n_ds <- 5
ds_gs <- sample(rownames(sce), n_ds)
g23 <- gids %in% c("g2", "g3")
assay(sce[ds_gs, g23]) <- assay(sce[ds_gs, g23]) * 100

pb <- aggregateData(sce, assay = "counts", fun = "sum")
ei <- metadata(sce)$experiment_info
design <- model.matrix(~ 0 + ei$group_id)
dimnames(design) <- list(ei$sample_id, levels(ei$group_id))
contrast <- limma::makeContrasts("g2-g1", "g3-g1", levels = design)

for (method in eval(as.list(args(pbDS))$method)) {
    test_that(paste("pbDS", method, sep = "."), {
        res <- pbDS(pb, 
            design = design, contrast = contrast,
            method = method, verbose = FALSE)
        
        expect_is(res, "list")
        expect_identical(length(res$table), ncol(contrast))
        expect_identical(names(res$table), colnames(contrast))
        expect_true(all(vapply(map(res$table, names), "==", 
            levels(kids), FUN.VALUE = logical(nlevels(kids)))))
        
        # check that top genes equal 'ds_gs' in ea. comparison & cluster
        top <- map_depth(res$table, 2, function(u) {
            dplyr::arrange(u, p_adj.loc) %>% 
                dplyr::slice(seq_len(n_ds)) %>% 
                pull("gene")
        }) %>% Reduce(f = "c")
        expect_true(all(unlist(map(top, setequal, ds_gs))))
    })
}

# global p-value adjustment ----------------------------------------------------
test_that(".p_adj_global", {
    cs <- paste0("c", seq_len(5))
    ks <- paste0("k", seq_len(8))
    ns <- sample(1e3, length(ks))
    names(cs) <- cs
    names(ks) <- ks
    names(ns) <- ks
    df <- lapply(cs, function(c)
        lapply(ks, function(k)
            data.frame(
                p_val = rgamma(ns[k], 0.1),
                p_adj.loc = rgamma(ns[k], 0.1))))
    df_adj <- .p_adj_global(df)
    for (c in cs)
        expect_identical(
            p.adjust(bind_rows(df[[c]])$p_val),
            bind_rows(df_adj[[c]])$p_adj.glb)
})
