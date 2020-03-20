#' @describeIn mmDS
#'
#' see details.
#'
#' @param dup_corr logical; whether to use
#'   \code{\link[limma:dupcor]{duplicateCorrelation}}.
#' @param trended logical; whether to use expression-dependent variance priors
#'  in \code{\link[limma]{eBayes}}.
#' @param ddf character string specifying the method for estimating
#'  the effective degrees of freedom. For \code{method = "dream"},
#'  either \code{"Satterthwaite"} (faster) or \code{"Kenward-Roger"}
#'  (more accurate); see \code{?variancePartition::dream} for details.
#'  For \code{method = "vst"}, method \code{"lme4"} is also valid;
#'  see \code{\link[lmerTest]{contest.lmerModLmerTest}}.
#'
#' @details
#' The \code{.mm_*} functions (e.g. \code{.mm_dream}) expect cells from a single
#'  cluster, and do not perform filtering or handle incorrect parameters well.
#' Meant to be called by \code{mmDS} with \code{method = c("dream", "vst")} and
#' \code{vst = c("sctransform", "DESeq2")} to be applied across all clusters.
#' \describe{
#' \item{\code{method = "dream2"}}{
#'   \code{variancePartition}'s (>=1.14.1) voom-lme4-implementation
#'   of mixed models for RNA-seq data; function \code{dream}.}
#' \item{\code{method = "dream"}}{
#'   \code{variancePartition}'s older voom-lme4-implementation
#'   of mixed models for RNA-seq data; function \code{dream}.}
#' \item{\code{method = "vst"}}{
#'   \describe{
#'   \item{\code{vst = "sctransform"}}{
#'     \code{lmer} or \code{blmer} mixed models on
#'     \code{\link[sctransform]{vst}} transformed counts.}
#'   \item{\code{vst = "DESeq2"}}{
#'     \code{\link[DESeq2]{varianceStabilizingTransformation}}
#'     followed by \code{lme4} mixed models.}}}}
#'
#' @importFrom BiocParallel MulticoreParam SerialParam
#' @importFrom edgeR DGEList
#' @importFrom dplyr %>% last mutate_at rename
#' @importFrom limma duplicateCorrelation eBayes topTable voom
#' @importFrom matrixStats rowSds
#' @importFrom scater computeLibraryFactors
#' @importFrom SingleCellExperiment counts sizeFactors
#' @importFrom stats as.formula model.matrix
#' @importFrom variancePartition dream getContrast
.mm_dream <- function(x,
    coef = NULL, covs = NULL,
    dup_corr = FALSE, trended = FALSE,
    ddf = c("Satterthwaite", "Kenward-Roger"),
    n_threads = 1, verbose = FALSE) {

    if (is.null(sizeFactors(x)))
        x <- computeLibraryFactors(x)

    ddf <- match.arg(ddf)
    x <- x[rowSds(as.matrix(counts(x))) > 0, ]
    y <- DGEList(counts(x), norm.factors = 1/sizeFactors(x))

    cd <- .prep_cd(x, covs)

    formula <- paste0("~", paste(c(covs, "group_id"), collapse = "+"))
    mm <- model.matrix(as.formula(formula), data = cd)
    v <- voom(y, mm)

    if (dup_corr) {
        dc <- duplicateCorrelation(v, mm, block = x$sample_id)
        v <- voom(y, mm, block = x$sample_id, correlation = dc$consensus)
    }

    if (n_threads > 1) {
        bp <- MulticoreParam(n_threads, progressbar = verbose)
    } else {
        bp <- SerialParam(progressbar = verbose)
    }

    formula <- paste0(formula, "+(1|sample_id)")
    if (verbose) print(formula)

    if (is.null(coef)) {
        coef <- last(colnames(mm))
        if (verbose)
            message("Argument 'coef' not specified; ",
                "testing for ", dQuote(coef), ".")
    }

    contrast <- getContrast(v, as.formula(formula), cd, coef)
    .dream <- expression(dream(v, formula, cd, contrast, ddf, 
        BPPARAM = bp, suppressWarnings = !verbose, quiet  = !verbose))
    if (verbose) fit <- eval(.dream) else suppressWarnings(fit <- eval(.dream))
    #fit <- eBayes(fit, trend = trended, robust = TRUE)

    topTable(fit, coef, number = Inf, sort.by = "none") %>%
        rename(p_val = "P.Value", p_adj.loc = "adj.P.Val")
}

#' @describeIn mmDS
#'
#' see details.
#'
#' @importFrom edgeR DGEList
#' @importFrom BiocParallel MulticoreParam SerialParam
#' @importFrom dplyr %>% last rename
#' @importFrom limma topTable
#' @importFrom matrixStats rowSds
#' @importFrom scater computeLibraryFactors
#' @importFrom SingleCellExperiment counts sizeFactors
#' @importFrom stats as.formula
#' @importFrom variancePartition dream voomWithDreamWeights
.mm_dream2 <- function(x, coef = NULL, covs = NULL,
    ddf = c("Satterthwaite", "Kenward-Roger"),
    n_threads = 1, verbose = FALSE) {

    if (is.null(sizeFactors(x)))
        x <- computeLibraryFactors(x)

    ddf <- match.arg(ddf)
    x <- x[rowSds(as.matrix(counts(x))) > 0, ]
    y <- DGEList(counts(x), norm.factors = 1 / sizeFactors(x))

    cd <- .prep_cd(x, covs)
    formula <- paste0(c("~(1|sample_id)", covs, "group_id"), collapse = "+")

    if (n_threads > 1) {
        bp <- MulticoreParam(n_threads, progressbar = verbose)
    } else {
        bp <- SerialParam(progressbar = verbose)
    }

    if (is.null(coef)) {
        coef <- paste0("group_id", last(levels(cd$group_id)))
        if (verbose)
            message("Argument 'coef' not specified; ",
                    "testing for ", dQuote(coef), ".")
    }

    if (verbose) print(formula)
    formula <- as.formula(formula)
    
    .dream <- expression(voomWithDreamWeights(y, 
        formula, cd, BPPARAM = bp, quiet = !verbose))
    if (verbose) v <- eval(.dream) else suppressMessages(v <- eval(.dream))
    res <- dream(v, formula, cd, BPPARAM = bp, ddf = ddf, 
        suppressWarnings = !verbose, quiet = !verbose)
    tbl <- topTable(res, coef = coef, Inf, sort.by = "none")
    rename(tbl, p_val = "P.Value", p_adj.loc = "adj.P.Val")
}

#' @describeIn mmDS
#'
#' see details.
#'
#' @param vst method to use as variance-stabilizing transformations.
#'   \code{"sctransform"} for \code{\link[sctransform]{vst}}; \code{"DESeq2"}
#'   for \code{\link[DESeq2]{varianceStabilizingTransformation}}.
#' @param bayesian logical; whether to use bayesian mixed models.
#' @param blind logical; whether to ignore experimental design for the vst.
#' @param REML logical; whether to maximize REML instead of log-likelihood.
#'
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom dplyr last
#' @importFrom purrr set_names
#' @importFrom SingleCellExperiment counts
.mm_vst <- function(x,
    vst = c("sctransform", "DESeq2"),
    coef = NULL, covs = NULL,
    bayesian = FALSE, blind = TRUE, REML = TRUE,
    ddf = c("Satterthwaite", "Kenward-Roger", "lme4"),
    n_threads = 1, verbose = FALSE) {

    vst <- match.arg(vst)
    ddf <- match.arg(ddf)
    cd <- .prep_cd(x, covs)
    y <- counts(x)

    # get formula
    formula <- paste(c("~(1|sample_id)", covs, "group_id"), collapse = "+")
    if (verbose) print(formula)
    formula <- as.formula(paste("u", formula))

    # get coefficient to test
    if (is.null(coef)) {
        coef <- paste0("group_id", last(levels(cd$group_id)))
        if (verbose)
            message("Argument 'coef' not specified; ",
                "testing for ", dQuote(coef), ".")
    }

    # fit mixed models for ea. gene
    fits <- bplapply(seq_len(nrow(y)), function(i)
        .fit_lmer(cbind(u = y[i, ], cd), formula, coef, bayesian, REML, ddf),
        BPPARAM = MulticoreParam(n_threads)) %>%
        set_names(rownames(y))

    if (verbose) message("Applying empirical Bayes moderation..")
    fits <- .mm_eBayes(fits, coef)
    i <- which(colnames(fits) == "p_val")
    fits$p_adj.loc <- p.adjust(fits$p_val)
    fits[, c(seq_len(i), ncol(fits), seq_len(ncol(fits)-1)[-seq_len(i)])]
}

# helper to prepare colData for .mm_dream/vst
#' @importFrom dplyr %>% mutate_at mutate_if
#' @importFrom SummarizedExperiment colData
.prep_cd <- function(x, covs) {
    cd <- colData(x)[c("sample_id", "group_id", covs)]
    cd <- data.frame(cd, check.names = FALSE)
    cd <- mutate_if(cd, is.factor, droplevels)
    if (!is.null(covs))
        cd <- mutate_at(cd, covs, function(u) 
            if (is.numeric(u)) scale(u))
    rownames(cd) <- colnames(x); cd
}

# fits mixed models and returns fit information required for eBayes
#' @importFrom blme blmer
#' @importFrom lmerTest lmer contest
#' @importFrom lme4 .makeCC lmerControl
#' @importFrom purrr map set_names
#' @importFrom stats df.residual residuals sd
.fit_lmer <- function(df, formula, coef, bayesian, REML, ddf) {
    # here we should do some handling of convergence/singularity
    fun <- ifelse(bayesian, blmer, lmerTest::lmer)
    mod <- tryCatch(fun(formula, df, REML, control = lmerControl(
        check.conv.singular = .makeCC(action = "ignore", tol = 1e-4))),
        error = function(e) e)
    if (inherits(mod, "error")) mod

    tryCatch({
        coefs <- colnames(coef(mod)[[1]])
        re <- cbind(coef(summary(mod)), p_val = NA_real_)
        re[, "p_val"] <- ifelse(
            "Pr(>|t|)" %in% colnames(re),
            re[, "Pr(>|t|)"], NA_real_)
        cs <- as.numeric(coefs == coef)
        ct <- lmerTest::contest(mod, cs, ddf = ddf)
        re[cs == 1, ncol(re)] <- ct[, ncol(ct)]
        re <- re[, c(seq_len(3), ncol(re))]
        split(re, col(re)) %>%
            map(set_names, coefs) %>%
            set_names(c("beta", "SE", "stat", "p_val")) %>%
            c(list(
                Amean = mean(df$u),
                sigma = sd(residuals(mod)),
                df.residual = df.residual(mod)))
    }, error = function(e) e)
}

#' @describeIn mmDS
#'
#' see details.
#' 
#' @param family character string specifying which GLMM to fit:
#'   \code{"poisson"} for \code{\link[blme]{bglmer}},
#'   \code{"nbinom"} for \code{\link[glmmTMB]{glmmTMB}}.
#' @param moderate logical; whether to perform empirical Bayes moderation.
#' 
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom dplyr %>% bind_rows last
#' @importFrom purrr set_names
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment assay
.mm_glmm <- function(x, coef = NULL, covs = NULL, n_threads = 1,
    family = c("poisson","nbinom"), verbose = TRUE, moderate = FALSE) {

    family <- match.arg(family)
    cd <- .prep_cd(x, covs)
    y <- counts(x)

    if (is.null(sizeFactors(x))) {
        cd$ls <- log(colSums(y))
    } else {
        cd$ls <- sizeFactors(x)
    }

    # get formula
    str <- c("~(1|sample_id)+offset(ls)", covs, "group_id")
    formula <- paste(str, collapse = "+")
    if (verbose) print(formula)
    formula <- as.formula(paste("u", formula))

    # get coefficient to test
    if (is.null(coef)) {
        coef <- paste0("group_id", last(levels(cd$group_id)))
        if (verbose)
            message("Argument 'coef' not specified; ",
                "testing for ", dQuote(coef), ".")
    }

    # fit mixed model for ea. gene
    names(gs) <- gs <- rownames(y)
    fits <- bplapply(gs, function(g) {
        df <- data.frame(u = y[g, ], cd)
        if (moderate) {
            .fit <- switch(family,
                nbinom = .fit_nbinom,
                poisson = .fit_bglmer)
            .fit(df, formula, coef)
        } else {
            tryCatch({
                switch(family,
                    nbinom = {
                        mod <- glmmTMB(formula, df,
                            family = nbinom1, REML = FALSE)
                        coef(summary(mod))[[1]][coef, ] },
                    poisson = {
                        mod <- bglmer(formula, df, family = "poisson")
                        coef(summary(mod))[coef, ] })
            }, error=function(e) rep(NA_real_, 4))
        }
    }, BPPARAM = MulticoreParam(n_threads, progressbar=verbose))

    if (moderate){
        if (verbose) message("Applying empirical Bayes moderation..")
        fits <- .mm_eBayes(fits, coef)
    } else {
        fits <- as.data.frame(t(bind_rows(fits)))
        colnames(fits) <- c("beta", "SE", "stat", "p_val")
    }
    i <- which(colnames(fits) == "p_val")
    fits$p_adj.loc <- p.adjust(fits$p_val)
    fits[, c(seq_len(i), ncol(fits), seq_len(ncol(fits)-1)[-seq_len(i)])]
}

.mm_poisson <- function(x, coef = NULL, covs = NULL,
    n_threads = 1, verbose = TRUE, moderate = FALSE)
    .mm_glmm(x, coef, covs, n_threads, family = "poisson",
        verbose = verbose, moderate = moderate)

.mm_nbinom <- function(x, coef = NULL, covs = NULL,
    n_threads = 1, verbose = TRUE, moderate = FALSE)
    .mm_glmm(x, coef, covs, n_threads, family = "nbinom",
        verbose = verbose, moderate = moderate)

#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom blme bglmer
#' @importFrom dplyr %>% bind_rows last rename
#' @importFrom glmmTMB glmmTMB nbinom1
#' @importFrom Matrix colSums t
#' @importFrom SingleCellExperiment counts sizeFactors
#' @importFrom stats model.matrix
#' @importFrom SummarizedExperiment colData
.mm_hybrid <- function(x, coef = NULL, covs = NULL, n_threads = 1, 
    verbose = TRUE, fam = c("nbinom","poisson"), th = 0.1) {

    fam <- match.arg(fam)
    x$cluster_id <- droplevels(x$cluster_id)
    cd <- .prep_cd(x, covs)
    y <- counts(x)

    if (is.null(sizeFactors(x))) {
        cd$ls <- log(colSums(y))
    } else {
        cd$ls <- sizeFactors(x)
    }

    # compute pseudobulks (sum of counts)
    pb <- aggregateData(x)
    pb_cd <- as.data.frame(colData(pb))

    # get sample-level formula
    str <- c("~group_id", intersect(covs, names(pb_cd)))
    formula <- as.formula(paste(str, collapse = "+"))
    mm <- model.matrix(formula, pb_cd)

    # get cell-level formula
    str <- c("~(1|sample_id)+offset(ls)", covs, "group_id")
    formula <- paste(str, collapse = "+")
    if (verbose) print(formula)
    formula <- as.formula(paste("u", formula))

    # get coefficient to test
    if (is.null(coef)) {
        coef <- last(colnames(mm))
        if (verbose)
            message("Argument 'coef' not specified; ",
                "testing for ", dQuote(coef), ".")
    }

    # pseudobulk analysis
    res <- pbDS(pb, design = mm, 
        coef = which(colnames(mm) == coef), 
        method = "edgeR", verbose = verbose)
    res <- res$table[[1]][[1]]
    cols <- c("F", "p_adj.loc", "coef", "p_adj.glb")
    cols_keep <- setdiff(colnames(res), cols)
    res <- rename(res[, cols_keep], p_val.pb = "p_val")
    row.names(res) <- res$gene
    gs <- as.character(res$gene)
    gs_keep <- res$p_val.pb < th
    names(gs) <- gs <- as.character(gs)[gs_keep]

    # fit mixed model for ea. gene below threshold in PB analysis
    fits <- bplapply(gs, function(i) {
        df <- data.frame(u = y[i, ], cd)
        tryCatch({
            switch(fam,
                nbinom = {
                    mod <- glmmTMB(formula, df, family = nbinom1, REML = FALSE)
                    coef(summary(mod))[[1]][coef, ]
                },
                poisson = {
                    mod <- bglmer(formula, df, family = "poisson")
                    coef(summary(mod))[coef, ]
                })
        }, error = function(e) rep(NA_real_, 4))
    }, BPPARAM = MulticoreParam(n_threads, progressbar = verbose))

    res$glmm.est <- res$p_val.glmm <- NA_real_
    res[gs, c("glmm.est", "p_val.glmm")] <- t(bind_rows(fits))[, c(1,4)]
    res$geomean.p_val <- res$mean.p_val <- res$p_val <- res$p_val.pb

    mat <- as.matrix(res[gs, c("p_val.pb", "p_val.glmm")])
    res[gs, "geomean.p_val"] <- 10^-rowMeans(-log10(mat))
    res[gs, "mean.p_val"] <- rowMeans(mat)
    res[gs, "p_val"] <- res[gs, "p_val.glmm"]

    res <- res[, -c(1, 2)]
    i <- which(colnames(res) == "p_val")
    res$p_adj.loc <- p.adjust(res$p_val)
    res[, c(seq_len(i), ncol(res), seq_len(ncol(res)-1)[-seq_len(i)])]
}

# fits negative binomial mixed models and
# returns fit information required for 'eBayes'
#' @importFrom glmmTMB glmmTMB nbinom1
#' @importFrom purrr map set_names
#' @importFrom stats coef df.residual residuals sd
.fit_nbinom <- function(df, formula, coef){
    mod <- tryCatch(
        glmmTMB(formula, family = nbinom1, data = df, REML = FALSE),
        error = function(e) { print(e); NULL })
    if (is.null(mod)) return(mod)
    tryCatch({
        coefs <- colnames(coef(mod)[[1]][[1]])
        re <- coef(summary(mod))[[1]]
        re <- split(re, col(re)) %>%
            map(set_names, coefs) %>%
            set_names(c("beta", "SE", "stat", "p_val"))
        c(re, list(
            Amean = mean(df$u),
            sigma = sd(residuals(mod)),
            df.residual = df.residual(mod)))
    }, error = function(e) NULL)
}


# fits poisson mixed models and
# returns fit information required for 'eBayes'
#' @importFrom blme bglmer
#' @importFrom purrr map set_names
#' @importFrom stats coef df.residual residuals sd
.fit_bglmer <- function(df, formula, coef){
    mod <- tryCatch(
        bglmer(formula, family="poisson", data=df),
        error = function(e) { print(e); NULL })
    if (is.null(mod)) return(mod)
    tryCatch({
        coefs <- colnames(coef(mod)[[1]])
        re <- coef(summary(mod))
        re <- split(re, col(re)) %>%
            map(set_names, coefs) %>%
            set_names(c("beta", "SE", "stat", "p_val"))
        c(re, list(
            Amean = mean(df$u),
            sigma = sd(residuals(mod)),
            df.residual = df.residual(mod)))
    }, error = function(e) NULL)
}

# formats a list of .fit_lmer results into
# an eBayes compatible list & performs moderation
#' @importFrom dplyr %>% bind_cols pull
#' @importFrom limma eBayes
#' @importFrom purrr map set_names
.mm_eBayes <- function(fits, coef, trended = TRUE) {
    rmv <- vapply(fits, inherits, what = "error", logical(1))
    f <- fits[!rmv]
    if (length(f) > 0) {
        vars <- set_names(names(f[[1]]))
        res <- lapply(vars, map, .x = f) %>%
            map(data.frame) %>% map(t) %>%
            map(function(u) {if (ncol(u) == 1) c(u) else u})
        nms <- c("coefficients", "stdev.unscaled", "z", "PValue")
        names(res)[seq_len(4)] <- nms
        res <- eBayes(res, trend = trended, robust = TRUE)
        res <- res[c("coefficients", "z", "PValue", "p.value")] %>%
            set_names(c("beta", "stat", "p_val0", "p_val")) %>%
            map(as.data.frame) %>% map(pull, coef) %>%
            data.frame(row.names = names(f))
    } else {
        res <- matrix(NA, nrow = 0, ncol = 4) 
        colnames(res) <- c("beta", "stat", "p_val0", "p_val")
        as.data.frame(res)
    }
    if (any(rmv)) {
        res[names(which(rmv)), "error"] <-
            vapply(fits[rmv], as.character, character(1))
    }
    res[names(fits), ]
}

# wrappers for variance-stabilizing transformations
# using the 'sctransform' and 'DESeq2' package, respectively
# ==============================================================================
#' @importFrom sctransform vst
#' @importFrom SingleCellExperiment counts
.vst_sctransform <- function(x, verbose) {
    sctransform::vst(counts(x), show_progress = verbose)$y
}
# ------------------------------------------------------------------------------
#' @importFrom DESeq2 varianceStabilizingTransformation
#' @importFrom SingleCellExperiment counts sizeFactors
.vst_DESeq2 <- function(x, covs, blind) {
    if (is.null(sizeFactors(x)))
        x <- computeLibraryFactors(x)
    covs <- paste(c(covs, "sample_id"), collapse = "+")
    formula <- as.formula(paste("~", covs))
    y <- counts(x)
    mode(y) <- "integer"
    y <- DESeqDataSetFromMatrix(y, colData(x), formula)
    sizeFactors(y) <- sizeFactors(x)
    if (!blind)
        y <- estimateDispersions(y)
    y <- varianceStabilizingTransformation(y, blind)
    assay(y)
}
