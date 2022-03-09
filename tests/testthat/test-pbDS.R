# load packages
suppressMessages({
    library(dplyr)
    library(purrr)
    library(SingleCellExperiment)
})

# generate toy dataset
seed <- as.numeric(format(Sys.time(), "%s"))
set.seed(seed)
sce <- .toySCE(c(100, 2e3))

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
ng <- length(gids <- levels(sce$group_id))

g3 <- sce$group_id == "g3"
g23 <- sce$group_id %in% c("g2", "g3")

# sample 'nde' genes & multiply counts by 10 for 'g2'- & 'g3'-cells
degs <- sample(rownames(sce), (nde <- 5))
assay(sce[degs, g23]) <- assay(sce[degs, g23]) * 10
pb <- aggregateData(sce, assay = "counts", fun = "sum")

# specify design & contrast matrix
ei <- metadata(sce)$experiment_info
design <- model.matrix(~ 0 + ei$group_id)
dimnames(design) <- list(ei$sample_id, levels(ei$group_id))
contrast <- limma::makeContrasts("g2-g1", "g3-g1", levels = design)

# default method settings ------------------------------------------------------
for (method in eval(as.list(args(pbDS))$method)) {
    test_that(paste("defaults - pbDS", method, sep = "."), {
        res <- pbDS(pb, 
            method = method, 
            min_cells = 0, 
            filter = "none", 
            verbose = FALSE)
        tbl <- res$table[[1]]
        expect_identical(names(tbl), kids)
        top <- map(tbl, function(u)
            dplyr::arrange(u, p_adj.loc) %>% 
                dplyr::slice(seq_len(5)) %>% 
                pull("gene"))
        expect_true(all(vapply(top, setequal, y = degs, logical(1))))
    })
}

# multiple contrast w/ & w/o 'treat' -------------------------------------------
for (method in eval(as.list(args(pbDS))$method)) {
    if (grepl("edgeR|limma", method)) {
        treat <- c(FALSE, TRUE)
    } else treat <- FALSE
    test_that(paste("pbDS", method, sep = "."), {
        for (t in treat) {
            res <- pbDS(pb, 
                treat = t, filter = "none",
                method = method, verbose = FALSE,
                design = design, contrast = contrast)
            
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
        }
    })
}

# global p-value adjustment ----------------------------------------------------
test_that(".p_adj_global", {
    names(cs) <- cs <- paste0("c", seq_len(5))
    names(ks) <- ks <- paste0("k", seq_len(8))
    ns <- sample(1e3, length(ks))
    names(ns) <- ks
    df <- lapply(ks, function(k) lapply(cs, function(c) {
        df <- data.frame(p_val = rgamma(ns[k], 0.1))
        df$p_adj.loc <- p.adjust(df$p_val, method = "BH")
        return(df) 
    })) %>% .p_adj_global
    for (c in cs)
        expect_identical(
            p.adjust(bind_rows(df[[c]])$p_val, method = "BH"),
            bind_rows(df[[c]])$p_adj.glb)
})
