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

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
ng <- length(gids <- levels(sce$group_id))

g3 <- sce$group_id == "g3"
g23 <- sce$group_id %in% c("g2", "g3")

# sample 'nde' genes & multiply counts by 10 for 'g3'-cells
foo <- sce
degs <- sample(rownames(foo), (nde <- 5))
assay(foo[degs, g3]) <- assay(foo[degs, g3]) * 10
pb <- aggregateData(foo, assay = "counts", fun = "sum")

for (method in eval(as.list(args(pbDS))$method)) {
    test_that(paste("defaults - pbDS", method, sep = "."), {
        res <- pbDS(pb, method = method, verbose = FALSE)
        tbl <- res$table[[1]]
        expect_identical(names(tbl), kids)
        top <- map(tbl, function(u)
            dplyr::arrange(u, p_adj.loc) %>% 
                dplyr::slice(seq_len(5)) %>% 
                pull("gene"))
        expect_true(all(vapply(top, setequal, y = degs, logical(1))))
    })
}

# sample 'nde' genes & multiply counts by 10 for 'g2'- & 'g3'-cells
foo <- sce
degs <- sample(rownames(foo), (nde <- 5))
assay(foo[degs, g23]) <- assay(foo[degs, g23]) * 10
pb <- aggregateData(foo, assay = "counts", fun = "sum")

# specify design & contrast matrix
ei <- metadata(sce)$experiment_info
design <- model.matrix(~ 0 + ei$group_id)
dimnames(design) <- list(ei$sample_id, levels(ei$group_id))
contrast <- limma::makeContrasts("g2-g1", "g3-g1", levels = design)

for (method in eval(as.list(args(pbDS))$method)) {
    test_that(paste("pbDS", method, sep = "."), {
        res <- pbDS(pb, 
            method = method, verbose = FALSE,
            design = design, contrast = contrast)
            
        expect_is(res, "list")
        expect_identical(length(res$table), ncol(contrast))
        expect_identical(names(res$table), colnames(contrast))
        expect_true(all(vapply(map(res$table, names), "==",
            levels(kids), FUN.VALUE = logical(nlevels(kids)))))
        
        # check that top genes equal 'degs' in ea. comparison & cluster
        top <- map_depth(res$table, 2, function(u) {
            dplyr::arrange(u, p_adj.loc) %>% 
                dplyr::slice(seq_len(nde)) %>% 
                pull("gene")
        }) %>% Reduce(f = "c")
        expect_true(all(unlist(map(top, setequal, degs))))
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
