context("Simulation of 'complex' scRNA-seq data")

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
    n_cells <- 1000
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
    expect_true(all(vapply(1:3, function(i) 
        ms[[i]] == n / length(ids[[i]]), logical(1))))
})
    