context("Check validity of input arguments")

library(SingleCellExperiment)
sce <- .toySCE()

test_that(".check_sce", {
    # object is not 'SingleCellExperiment'
    x <- sce; class(x) <- "x"
    expect_error(.check_sce(x))
    # "group_id" 'colData' column not required
    x <- sce; x$group_id <- NULL
    expect_silent(.check_sce(x, FALSE))
    # missing "x_id" 'colData' column;
    # x = "sample", "cluster", "group"
    for (i in colnames(colData(sce))) {
        x <- sce; x[[i]] <- NULL
        expect_error(.check_sce(x))
    }
})

test_that(".check_arg_assay", {
    expect_error(.check_arg_assay(sce, 1))
    expect_error(.check_arg_assay(sce, "x"))
    expect_error(.check_arg_assay(sce, c(assayNames(sce)[1], "x")))
    expect_silent(.check_arg_assay(sce, assayNames(sce)[1]))
})

test_that(".check_args_simData", {
    u <- list(x = sce, ng = 10, nc = 100, ns = 3, nk = 2, probs = NULL, 
        p_dd = diag(6)[1, ], p_type = 0.1, lfc = 1, rel_lfc = NULL, 
        p_ep = 0.1, p_dp = 0.1, p_dm = 0.1, paired = FALSE, force = TRUE,
        phylo_tree = NULL, phylo_pars = c(0, 3))
    expect_silent(.check_args_simData(u))
    # 'ng', 'nc', 'ns', 'nk' should be single numerics > 0
    for (arg in c("ng", "nc", "ns", "nk")) {
        v <- u; v[[arg]] <- 1; expect_silent(.check_args_simData(v))
        for (val in list(numeric(2)+1, numeric(1), NA, "a"))
            v[[arg]] <- val; expect_error(.check_args_simData(v))
    }
    # 'p_dd' should be length of DD categories, sum to 1, and be in [0,1]
    v <- u; for (val in list(rep(1,5)/5, rep(1,6)/5, c(-1,2,rep(0,4))))
        v$p_dd <- val; expect_error(.check_args_simData(v))
    # 'paired' should be length-one logical
    v <- u; for (val in list("x", 123, c(TRUE, FALSE)))
        v$paired <- val; expect_error(.check_args_simData(v))
    # 'rel_lfc' should be >0 and be of length 'nk'
    u$rel_lfc <- rep(1, u$nk)
    expect_silent(.check_args_simData(u))
    v <- u; for (val in list(-1, rep(1,u$nk+1)))
        v$rel_lfc <- val; expect_error(.check_args_simData(v))
})
