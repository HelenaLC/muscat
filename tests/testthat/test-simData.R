suppressMessages({
    library(dplyr)
    library(SingleCellExperiment)
})

data(example_sce)
ref <- prepSim(example_sce, verbose = FALSE)
ng <- 200; nc <- 2e3; ns <- 3; nk <- 2

test_that("no groups, no samples, no clusters", {
    s <- sample(levels(ref$sample_id), 1)
    k <- sample(levels(ref$cluster_id), 1)
    sub <- ref[, ref$sample_id == s & ref$cluster_id == k]
    sim <- simData(sub, ng, nc, dd = FALSE, force = TRUE)
    expect_null(sim$group_id)
    expect_null(sim$sample_id)
    expect_null(sim$cluster_id)
})

test_that("no groups, yes samples and clusters", {
    sids <- levels(ref$sample_id)
    kids <- levels(ref$cluster_id)
    sim <- simData(ref, ng, nc, dd = FALSE, force = TRUE)
    expect_null(sim$group_id)
    expect_true(setequal(metadata(sim)$ref_sids, sids))
    expect_true(setequal(metadata(sim)$ref_kids, kids))
})

test_that("only samples", {
    k <- sample(levels(ref$cluster_id), 1)
    sub <- ref[, ref$cluster_id == k]
    sim <- simData(sub, ng, nc, dd = FALSE, force = TRUE)
    expect_null(sim$group_id)
    expect_null(sim$culster_id)
    expect_setequal(metadata(sim)$ref_sids, levels(sub$sample_id))
})

test_that("only clusters", {
    s <- sample(levels(ref$sample_id), 1)
    sub <- ref[, ref$sample_id == s]
    sim <- simData(sub, ng, nc, dd = FALSE, force = TRUE)
    expect_null(sim$group_id)
    expect_null(sim$sample_id)
    expect_setequal(metadata(sim)$ref_kids, levels(sub$cluster_id))
})


test_that("pbDS() gets at least 50% right for 10% DE genes", {
    replicate(3, {
        sim <- simData(ref, ng, nc, ns, nk, 
            p_dd = c(0.9,0,0.1,0,0,0), force = TRUE)
        gi <- metadata(sim)$gene_info
        pb <- aggregateData(sim)
        re <- pbDS(pb, verbose = FALSE)
        re <- bind_rows(re$table[[1]])
        re <- arrange(re, p_adj.loc)
        gi$id <- with(gi, paste0(gene, cluster_id))
        re$id <- with(re, paste0(gene, cluster_id))
        n_de <- sum(de <- gi$category == "de")
        expect_gte(sum(re$id[seq_len(n_de)] %in% gi$id[de]), round(n_de/2))
    })
})

test_that("Single group simulation", {
    gs <- c("A", "B")
    ps <- list(c(1, 0), c(0, 1))
    names(ps) <- gs
    for (g in gs) {
        x <- simData(ref, 
            nc, ns, nk, ng = 10, force = TRUE,
            probs = list(NULL, NULL, ps[[g]]))
        expect_identical(levels(x$group_id), g)
        sids <- sprintf("sample%s.%s", seq_len(ns), g)
        expect_identical(levels(x$sample_id), sids)
        gi <- metadata(x)$gene_info
        ms <- paste0("sim_mean.", setdiff(gs, g))
        expect_true(all(is.na(gi[[ms]])))
    }
})

test_that("Pure simulations give single DD category", {
    for (c in cats) {
        sim <- simData(ref, nc, ns, nk, 
            p_dd = as.numeric(cats == c),
            ng = ng, force = TRUE)
        gi <- metadata(sim)$gene_info
        expect_equal(table(gi$category)[[c]], ng*nk)
    }
})

test_that("simData() - input arguments", {
    # number of genes mismatch b/w simulation & reference
    expect_error(simData(ref, ng = 100, force = FALSE))
    expect_silent(simData(ref, ng = 100, force = TRUE))
    nk <- length(kids <- levels(ref$cluster_id))
    # named 'rel_lfc's (mis)match cluster names
    lfc <- rep(1, nk); names(lfc) <- kids; lfc2 <- lfc; names(lfc2)[1] <- "x"
    expect_silent(simData(ref, nk = nk, rel_lfc = lfc, ng = 10, force = TRUE))
    expect_error(simData(ref, nk = nk, rel_lfc = lfc2, ng = 10, force = TRUE))
})

test_that("Type genes & cluster phylogeny", {
    ng <- 1e3; pt <- 0.05
    t <- "(('cluster1':0.1,'cluster2':0.1):0.4,'cluster3':0.5);"
    args <- list(
        list(pt = pt, t = NULL), # type genes but no phylogeny
        list(t = t, pt = 0))     # both, type genes and phylogeny
    cs <- c("type", "state")     # possible gene classes
    for (i in seq_along(args)) {
        rd <- rowData(x <- simData(ref, ng = ng, force = TRUE,
            p_type = args[[i]]$pt, phylo_tree = args[[i]]$t))
        expect_is(rd, "DataFrame")
        expect_equal(dim(rd), c(ng, 2))
        expect_equal(colnames(rd), c("class", "specs"))
        ns <- as.list(table(rd$class))
        is_state <- rd$class == "state"
        expect_true(all(is.na(unlist(rd$specs[is_state]))))
        expect_true(all(unlist(rd$specs[!is_state]) %in% levels(x$cluster_id)))
    }
    # missing clusters in input phylogeny
    t <- "(('cluster1':0.1,'cluster2':0.1):0.4,'cluster4':0.5);"
    expect_error(simData(ref, phylo_tree = t))
    # specification of both, 'p_type' and 'phylo_tree'
    expect_error(simData(ref, p_type = pt, phylo_tree = t))
})
