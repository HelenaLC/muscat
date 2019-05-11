#' @describeIn mmDS
#' 
#' see details.
#' 
#' @param dup_corr logical; whether to use 
#'   \code{\link[limma]{duplicateCorrelation}}.
#' @param trended logical; whether to use expression-dependent variance priors
#'  in \code{\link[limma]{eBayes}}.
#' @param ddf character string specifying the method for estimating 
#'  the effective degrees of freedom. For \code{method = "dream"}, 
#'  either \code{"Satterthwaite"} (faster) or \code{"Kenward-Roger"} 
#'  (more accurate); see \code{\link[variancePartition]{dream}}.
#'  For \code{method = "vst"}, method \code{"lme4"} is also valid;
#'  see \code{\link[lmerTest]{contest.lmerModLmerTest}}.
#'   
#' @details 
#' \code{.mm_dream} and \code{.mm_vst} expect cells from a single cluster, 
#' and do not perform filtering or handle incorrect parameters well. 
#' Meant to be called by \code{mmDS} with \code{method = c("dream", "vst")} and
#' \code{vst = c("sctransform", "DESeq2")} to be applied across all clusters.
#' \describe{
#' \item{\code{method = "dream"}}{
#'   voom-lme4-implementation \code{\link[variancePartition]{dream}} 
#'   of mixed models for RNAseq data.}
#' \item{\code{method = "vst"}}{
#'   \describe{
#'   \item{\code{vst = "sctransform"}}{
#'     \code{lmer} or \code{blmer} mixed models on 
#'     \code{\link[sctransform]{vst}} transformed counts.}
#'   \item{\code{vst = "DESeq2"}}{
#'     \code{\link[DESeq2]{varianceStabilizingTransformation}} 
#'     followed by \code{lme4} mixed models.}}}}
#'   
#' @importFrom edgeR DGEList
#' @importFrom doParallel registerDoParallel
#' @importFrom dplyr %>% last mutate_at rename
#' @importFrom limma duplicateCorrelation eBayes topTable voom
#' @importFrom magrittr set_rownames
#' @importFrom matrixStats rowSds
#' @importFrom parallel makeCluster stopCluster
#' @importFrom scran computeSumFactors
#' @importFrom SingleCellExperiment counts sizeFactors
#' @importFrom stats as.formula model.matrix
#' @importFrom variancePartition dream getContrast
.mm_dream <- function(x, 
    coef, covs, n_threads, verbose, 
    dup_corr = FALSE, trended = TRUE, 
    ddf = c("Satterthwaite", "Kenward-Roger")) {
    
    if (is.null(sizeFactors(x)))
        x <- computeSumFactors(x)
    
    ddf <- match.arg(ddf)
    x <- x[rowSds(as.matrix(counts(x))) > 0, ]
    y <- DGEList(counts(x), norm.factors = 1 / sizeFactors(x))
    
    cd <- .prep_cd(x, covs)
    
    formula <- paste0("~", paste(c(covs, "group_id"), collapse = "+"))
    mm <- model.matrix(as.formula(formula), data = cd)
    v <- voom(y, mm)
    
    if (dup_corr) {
        dup_corr <- duplicateCorrelation(v, mm, block = x$sample_id)
        v <- voom(y, mm, block = x$sample_id, correlation = dup_corr$consensus)
    }
    
    if (n_threads > 1) {
        cl <- makeCluster(n_threads)
        registerDoParallel(cl)
    }
    
    formula <- paste0(formula, "+(1|sample_id)")
    if (verbose) print(formula)
    
    if (is.null(coef)) {
        coef <- last(colnames(mm))
        if (verbose) 
            message("Argument 'coef' not specified; ", 
                sprintf("testing for %s.", dQuote(coef)))
    }
    
    contrast <- getContrast(v, as.formula(formula), cd, coef)
    fit <- dream(v, formula, cd, contrast, ddf = ddf)
    fit <- eBayes(fit, trend = trended, robust = TRUE)
    if (n_threads > 1) stopCluster(cl)
    
    topTable(fit, number = Inf, sort.by = "none") %>% 
        rename(p_val = "P.Value", p_adj.loc = "adj.P.Val")
}

#' @describeIn mmDS
#' 
#' see details.
#' 
#' @param vst method to use as variance-stabilizing transformations. 
#'   \code{"sctransform"} for \code{\link[sctransform]{vst}};
#'   \code{"DESeq2"} for \code{\link[DESeq2]{varianceStabilizingTransformation}}.
#' @param bayesian logical; whether to use bayesian mixed models.
#' @param blind logical; whether to ignore experimental design for the vst.
#' @param REML logical; whether to maximize REML instead of log-likelihood.
#' 
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateDispersions 
#'   sizeFactors varianceStabilizingTransformation
#' @importFrom dplyr last
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom sctransform vst
#' @importFrom SingleCellExperiment counts sizeFactors sizeFactors<-
#' @importFrom SummarizedExperiment assay
#' @importFrom stats p.adjust
.mm_vst <- function(x, 
    coef, covs, n_threads, verbose, 
    vst = c("sctransform", "DESeq2"),
    bayesian = FALSE, blind = TRUE, REML = TRUE,
    ddf = c("Satterthwaite", "Kenward-Roger", "lme4")) {
    
    vst <- match.arg(vst)
    ddf <- match.arg(ddf)
    cd <- .prep_cd(x, covs)
    
    # variance-stabilizing transformation
    y <- switch(vst, 
        sctransform = {
            # assure correct vst() function is used
            fun <- getFromNamespace("vst", vst)
            fun_call <- expression(fun(
                umi = assay(x), 
                cell_attr = data.frame(cd,
                    logUMI = log10(colSums(assay(x)) + 1)),
                #logDRT = log10(colSums(assay(x) > 0))),
                latent_var = 'logUMI',#c('logUMI', 'logDRT'),
                batch_var = 'sample_id', show_progress = verbose)$y)
            if (verbose) eval(fun_call) else suppressMessages(eval(fun_call))
        }, 
        DESeq2 = {
            if (is.null(sizeFactors(x)))
                x <- computeSumFactors(x)
            formula <- paste("~", paste(c(covs, "sample_id"), collapse="+"))
            formula <- as.formula(formula)
            y <- as.matrix(counts(x))
            y <- suppressMessages(DESeqDataSetFromMatrix(y, cd, formula))
            sizeFactors(y) <- sizeFactors(x)
            if (!blind) y <- estimateDispersions(y)
            varianceStabilizingTransformation(y, blind)
        })
    
    # get formula
    formula <- paste(c("~(1|sample_id)", covs, "group_id"), collapse = "+")
    if (verbose) print(formula)
    formula <- as.formula(paste("u", formula))
    
    # get coefficient to test
    if (is.null(coef)) {
        coef <- paste0("group_id", last(levels(x$group_id)))
        if (verbose) 
            message("Argument 'coef' not specified; ", 
                sprintf("testing for %s.", dQuote(coef)))
    }
    
    # fit mixed model for ea. gene
    fits <- bplapply(seq_len(nrow(y)), function(i) {
        .fit_lmer(df = data.frame(u = y[i, ], cd), 
            formula, coef, bayesian, REML, ddf)
    }, BPPARAM = MulticoreParam(n_threads)) %>% 
        set_names(rownames(y))
    
    if (verbose) message("Applying empirical Bayes moderation...")
    res <- .eBayes(fits, coef, trended)
    res$p_adj.loc <- p.adjust(res$p_val)#, method = "BH")
    return(res)
}

# helper to prepare colData for .mm_dream/vst
#' @importFrom dplyr %>% mutate_at
#' @importFrom magrittr set_rownames
#' @importFrom SummarizedExperiment colData
.prep_cd <- function(x, covs) {
    cd <- colData(x)[c("sample_id", "group_id", covs)]
    data.frame(cd, check.names = FALSE) %>% 
        mutate_at(covs, function(u) if (is.numeric(u)) scale(u)) %>% 
        set_rownames(colnames(x))
}

# fits mixed models and returns fit information required for eBayes
#' @importFrom blme blmer
#' @importFrom dplyr last
#' @importFrom lmerTest contest
#' @importFrom lme4 lmer lmerControl
#' @importFrom purrr set_names
#' @importFrom stats residuals
.fit_lmer <- function(df, formula, coef, 
    bayesian = FALSE, REML = TRUE, ddf = "Kenward-Roger") {
    
    # here we should do some handling of convergence/singularity
    fun <- ifelse(bayesian, blmer, getFromNamespace("lmer", "lme4"))
    mod <- tryCatch(fun(formula, df, REML, control = lmerControl(
        check.conv.singular = .makeCC(action = "ignore", tol = 1e-4))),
        error = function(e) NULL)
    if (is.null(mod)) return(mod)
    
    tryCatch({
        coefs <- colnames(coef(mod)[[1]])
        cs <- as.numeric(coefs == coef)
        re <- cbind(coef(summary(mod)), p_val = NA_real_)
        re[, "p_val"] <- ifelse(
            "Pr(>|t|)" %in% colnames(re), 
            re[, "Pr(>|t|)"], NA_real_)
        re[cs == 1, ncol(re)] <- last(contest(mod, cs, ddf = ddf))
        re <- re[, c(seq_len(3), ncol(re))]
        re <- split(re, col(re)) %>% 
            map(set_names, coefs) %>% 
            set_names(c("beta", "SE", "stat", "p_val"))
        c(re, list(
            Amean = mean(df$u),
            sigma = sd(residuals(mod)), 
            df.residual = df.residual(mod)))
    }, error=function(e) NULL)
}

# formats a list of .fit_lmer results into 
# an eBayes compatible list & performs moderation
#' @importFrom dplyr %>% bind_cols
#' @importFrom limma eBayes
#' @importFrom purrr map set_names
.eBayes <- function(fits, coef, trended = TRUE) {
    fits <- fits[!vapply(fits, is.null, logical(1))]
    vars <- names(fits[[1]]) %>% set_names()
    res <- lapply(vars, map, .x = fits) %>% 
        map(data.frame) %>% map(t) %>% 
        map(function(u) {if (ncol(u) == 1) c(u) else u})
    names(res)[seq_len(4)] <- c("coefficients", "stdev.unscaled", "z", "PValue")
    res <- eBayes(res, trend = trended, robust = TRUE)
    res[c("coefficients", "z", "PValue", "p.value")] %>% 
        set_names(c("beta", "stat", "p_val0", "p_val")) %>% 
        map(as.data.frame) %>% map(pull, coef) %>% 
        bind_cols %>% data.frame(row.names = names(fits))
}
