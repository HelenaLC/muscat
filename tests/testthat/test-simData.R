context("Simulation of 'complex' scRNA-seq data")

suppressMessages({
    library(dplyr)
    library(SingleCellExperiment)
})

k <- paste0("cluster", seq_len(5))
s <- paste0("sample", seq_len(4))
g <- paste0("group", seq_len(3))

test_that(".sample_n_cells", {
    n <- 50
    x <- .sample_n_cells(n, k, s)
    expect_true(all(unlist(x) == n))
    expect_identical(rownames(x), s)
    expect_identical(colnames(x), k)
    expect_true(all(vapply(x, length, numeric(1)) == 2))
    
    n <- c(10, 100)
    x <- .sample_n_cells(n, k, s)
    rng <- vapply(x, range, numeric(2))
    expect_true(all(rng[1, ] >= n[1]))
    expect_true(all(rng[2, ] <= n[2]))
})

test_that(".split_cells", {
    n_cells <- 1e3
    cells <- paste0("cell", seq_len(n_cells))
    x <- matrix(0, 
        nrow = 1, ncol = n_cells, 
        dimnames = list(NULL, cells))
    cd <- data.frame(
        row.names = cells,
        cluster_id = sample(k, n_cells, TRUE),
        sample_id = sample(s, n_cells, TRUE))
    
    cs <- .split_cells(cd, "cluster_id")
    expect_identical(names(cs), k)
    expect_identical(
        as.numeric(vapply(cs, length, numeric(1))),
        as.numeric(table(cd$cluster_id)))
    
    cs <- .split_cells(cd, c("cluster_id", "sample_id"))
    expect_identical(names(cs), k)
    nms_lvl2 <- vapply(cs, names, character(length(s)))
    expect_true(all(apply(nms_lvl2, 2, identical, s)))
    
    cs <- .split_cells(cd, c("sample_id", "cluster_id"))
    expect_identical(names(cs), s)
    nms_lvl2 <- vapply(cs, names, character(length(k)))
    expect_true(all(apply(nms_lvl2, 2, identical, k)))
})

test_that(".sample_cell_md", {
    n <- 1e3
    ids <- list(k, s, g)
    md <- .sample_cell_md(n, ids)
    ns <- apply(md, 2, table)
    ms <- vapply(ns, mean, numeric(1))
    expect_true(all(vapply(seq_along(ids), function(i) 
        ms[[i]] == n/length(ids[[i]]), logical(1))))
})

data(sce)
ref <- prepSim(sce, verbose = FALSE)
ng <- 200; nc <- 2e3; ns <- 3; nk <- 2

test_that("simData() - 'paired = TRUE/FALSE'", {
    replicate(5, {
        md <- metadata(simData(ref, paired = TRUE, ng = ng, force = TRUE))
        expect_true(all(apply(md$ref_sids, 1, 
            function(u) length(unique(u)) == 1)))
    })
    foo <- replicate(5, {
        md <- metadata(simData(ref, paired = FALSE, ng = ng, force = TRUE))
        all(apply(md$ref_sids, 1, function(u) length(unique(u)) == 1))
    })
    expect_true(!all(foo))
})

test_that("pbDS() gets at least 50% right for 10% DE genes", {
    replicate(5, {
        sim <- simData(ref, nc, ns, nk, 
            p_dd = c(0.9,0,0.1,0,0,0),
            ng = ng, force = TRUE)
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

test_that("Pure simulations give single DD category", {
    for (c in cats) {
        sim <- simData(ref, nc, ns, nk, 
            p_dd = as.numeric(cats == c),
            ng = ng, force = TRUE)
        gi <- metadata(sim)$gene_info
        expect_equal(table(gi$category)[[c]], ng*nk)
    }
})

test_that("Number of genes mismatch b/w simulation & reference", {
    expect_error(simData(ref, ng = 100, force = FALSE))
    expect_silent(simData(ref, ng = 100, force = TRUE))
})

test_that("Type genes & cluster phylogeny", {
    ng <- 1e3; pt <- 0.05
    t <- "(('cluster1':0.1,'cluster2':0.1):0.4,'cluster3':0.5);"
    args <- list(
        list(pt = pt, t = NULL), # type genes but no phylogeny
        list(t = t, pt = 0),     # phylogeny but no type genes
        list(t = t, pt = pt))    # both, type genes and phylogeny
    cs <- c("type", "state", "shared") # possible gene classes
    cs_ex <- c("shared", "type", NULL) # these shouldn't appear
    for (i in seq_along(args)) {
        rd <- rowData(x <- simData(ref, ng = ng, force = TRUE,
            p_type = args[[i]]$pt, phylo_tree = args[[i]]))
        expect_is(rd, "DataFrame")
        expect_equal(dim(rd), c(ng, 2))
        expect_equal(colnames(rd), c("class", "specs"))
        ns <- as.list(table(rd$class))
        expect_true(all(unlist(ns[[cs_ex[i]]]) == 0))
        is_state <- rd$class == "state"
        expect_true(all(is.na(unlist(rd$specs[is_state]))))
        expect_true(all(unlist(rd$specs[!is_state]) %in% levels(x$cluster_id)))
    }
    # missing clusters in input phylogeny
    t <- "(('cluster1':0.1,'cluster2':0.1):0.4,'cluster4':0.5);"
    expect_error(simData(ref, phylo_tree = t))
})
