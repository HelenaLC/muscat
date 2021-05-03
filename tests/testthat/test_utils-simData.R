nk <- length(k <- paste0("cluster", seq_len(5)))
ns <- length(s <- paste0("sample", seq_len(4)))
ng <- length(g <- paste0("group", seq_len(3)))

test_that(".get_ns", {
    replicate(10, {
        # default to using all samples when 
        # dd = FALSE or dd = TRUE, paired = TRUE
        ns_ref <- sample(100, 1)
        args <- list(
            list(dd = TRUE, paired = TRUE),
            list(dd = FALSE, paired = TRUE),
            list(dd = FALSE, paired = FALSE))
        for (. in args) {
            .$ns_ref <- ns_ref
            .$force <- FALSE
            . <- c(., list(ns_sim = NULL))
            ns <- do.call(.get_ns, .)
            expect_identical(ns_ref, ns)
        }
        # error when force = FALSE and 
        # desired number of samples to simulate 
        # exceeds available reference samples
        ns_sim <- ns_ref + 1
        args <- list(
            list(dd = TRUE, paired = TRUE),
            list(dd = TRUE, paired = FALSE),
            list(dd = FALSE, paired = TRUE),
            list(dd = FALSE, paired = FALSE))
        for (. in args) {
            .$ns_ref <- ns_ref
            .$ns_sim <- ns_sim
            .$force <- FALSE
            expect_error(do.call(.get_ns, .))
            # pass when force = TRUE
            .$force <- TRUE
            ns <- do.call(.get_ns, .)
            expect_identical(ns_sim, ns)
        }
    })
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
    ids <- list(k, s, g)
    md <- .sample_cell_md((n <- 1e3), ids)
    ms <- vapply(apply(md, 2, table), mean, numeric(1))
    expect_true(all(vapply(seq_along(ids), function(i) 
        ms[[i]] == n/length(ids[[i]]), logical(1))))
    set.seed(1); a <- .sample_cell_md(n, ids)
    set.seed(1); b <- .sample_cell_md(n, ids, 
        list(rep(1/nk,nk),rep(1/ns,ns),rep(1/ng,ng)))
    expect_identical(a, b)
})
