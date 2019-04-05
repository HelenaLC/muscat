#' @describeIn mmDS
#' 
#' A SCE wrapper around the voom-lme4-implementation 
#' \code{\link[variancePartition]{dream}} of mixed models for RNAseq data. 
#' \code{.mm_dream} expects cells from a single cluster, and
#' does not perform filtering or handle incorrect parameters well.
#' Meant to be called by \code{mmDS} with \code{method = "dream"}
#' to be applied across all clusters.
#' 
#' @param dup_corr logical; whether to use 
#'   \code{\link[limma]{duplicateCorrelation}}.
#' @param trended logical; whether to use expression-dependent variance priors
#'  in \code{\link[limma]{eBayes}}.
#' @param ddf character string specifying the method for estimating 
#'  the effective degrees of freedom. Either \code{"Satterthwaite"}
#'  (faster) or \code{"Kenward-Roger"} (more accurate) or . 
#'  See \code{\link[variancePartition]{dream}} for details.
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
    fit <- eBayes(fit, trend=trended, robust=TRUE)
    if (n_threads > 1) stopCluster(cl)
    
    topTable(fit, number = Inf, sort.by = "none") %>% 
        rename(p_val = "P.Value", p_adj.loc = "adj.P.Val")
}

#' @describeIn mmDS
#' 
#' A SCE wrapper applying \code{lmer} or \code{blmer} mixed models on 
#' \code{\link[sctransform]{vst}} transformed counts.
#' \code{.mm_vst} expects cells from a single cluster, and 
#' does not perform filtering or handle incorrect parameters well. 
#' Meant to be called by \code{mmDS} with \code{method = "vst"}
#' to be applied across all clusters.
#' 
#' @param ddf logical; methods for estimating degrees of freedom; either
#'  'Kenward-Roger' (default, more accurate) or 'Satterthwaite' (faster).
#' @param REML logical; whether to maximize REML instead of log-likelihood.
#' @param bayesian logical; whether to use bayesian mixed models.
#' 
#' @importFrom sctransform vst
#' @importFrom dplyr last
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom SingleCellExperiment counts sizeFactors sizeFactors<-
#' @importFrom SummarizedExperiment assay
#' @importFrom stats p.adjust
.mm_vst <- function(x, coef, covs, n_threads, verbose, 
                    ddf = "Kenward-Roger", REML = TRUE, bayesian=FALSE ) {
    y <- vst(assay(x), show_progress=verbose)$y
    
    cd <- .prep_cd(x, covs)
    
    formula <- paste0("~(1|sample_id)+", 
                      paste(c(covs, "group_id"), collapse="+"))
    if (verbose) print(formula)
    formula <- as.formula(paste0("x", formula))
    
    if (is.null(coef)) {
        coef <- paste0("group_id", last(levels(x$group_id)))
        if (verbose) 
            message("Argument 'coef' not specified; ", 
                    sprintf("testing for %s.", dQuote(coef)))
    }
    
    # we fit mixed models on each gene
    fits <- bplapply( y, form=formula, df=cd, testcoef=coef, REML=REML, ddf=ddf, 
                      BPPARAM=MulticoreParam(n_threads), bayesian=bayesian, FUN=.fitlmer)
    
    if(verbose) message("Applying empirical Bayes moderation")
    res <- .mmEBayesWrapper(fits, coef, trended)
    
    res$p_adj.loc <- p.adjust(res$p_val, method = "BH")
    return(res)
}

#' @describeIn mmDS
#' 
#' A SCE wrapper around \code{\link[DESeq2]{varianceStabilizingTransformation}} 
#' followed by \code{lme4} mixed models. 
#' \code{.mm_deseq} expects cells from a single cluster, and 
#' does not perform filtering or handle incorrect parameters well. 
#' Meant to be called by \code{mmDS} with \code{method = "deseq"}
#' to be applied across all clusters.
#' 
#' @param ddf logical; methods for estimating degrees of freedom; either
#'  'Kenward-Roger' (default, more accurate) or 'Satterthwaite' (faster).
#' @param blind logical; whether to ignore experimental design for the vst.
#' @param REML logical; whether to maximize REML instead of log-likelihood.
#' @param bayesian logical; whether to use bayesian mixed models.
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateDispersions 
#'   sizeFactors varianceStabilizingTransformation
#' @importFrom dplyr last
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom SingleCellExperiment counts sizeFactors sizeFactors<-
#' @importFrom SummarizedExperiment assay
#' @importFrom stats p.adjust
.mm_deseq <- function(x, coef, covs, n_threads, verbose, ddf="Satterthwaite", 
                      bayesian=FALSE, blind = TRUE, REML = TRUE) {
    
    if (is.null(sizeFactors(x)))
        x <- computeSumFactors(x)
    
    cd <- .prep_cd(x, covs)
    formula <- paste0("~", paste(c(covs, "sample_id"), collapse="+"))
    formula <- as.formula(formula)
    y <- suppressMessages(DESeqDataSetFromMatrix(
        as.matrix(counts(x)), cd, formula))
    
    sizeFactors(y) <- sizeFactors(x)
    if (!blind) y <- estimateDispersions(y)
    vst <- varianceStabilizingTransformation(y, blind)
    
    formula <- paste(c("~(1|sample_id)", covs, "group_id"), collapse="+")
    if (verbose) print(formula)
    formula <- as.formula(paste("x", formula))
    
    if (is.null(coef)) {
        gids <- levels(x$group_id)
        coef <- paste0("group_id", last(gids))
        if (verbose) 
            message("Argument 'coef' not specified; ", 
                    sprintf("testing for %s.", dQuote(coef)))
    }

    # we fit mixed models on each gene
    fits <- bplapply( y, 1, form=formula, df=cd, testcoef=coef, REML=REML, ddf=ddf, 
                      BPPARAM=MulticoreParam(n_threads), bayesian=bayesian, FUN=.fitlmer)

    if(verbose) message("Applying empirical Bayes moderation")
    res <- .mmEBayesWrapper(fits, coef, trended)
    res$p_adj.loc <- p.adjust(res$p_val, method = "BH")
    return(res)
}

# helper to prepare colData for .mm_dream/vst
#' @importFrom dplyr %>% mutate_at
#' @importFrom methods is
#' @importFrom magrittr set_rownames
#' @importFrom SummarizedExperiment colData
.prep_cd <- function(x, covs) {
    cd <- colData(x)[c("sample_id", "group_id", covs)]
    data.frame(cd, check.names = FALSE) %>% 
        mutate_at(covs, function(u) if (is.numeric(u)) scale(u)) %>% 
        set_rownames(colnames(x))
}

# fits mixed models and returns fit information required for eBayes
#' @import lmerTest Matrix
#' @importFrom blme blmer
#' @importFrom dplyr last
.fitlmer <- function(x, form, df, testcoef, bayesian=FALSE, ddf="Kenward-Roger", REML=TRUE){
    df$x <- x
    mod <- tryCatch({
        # here we should do some handling of convergence/singularity
        if(bayesian){
            blmer(form, df, REML=REML)
        }else{
            lmer(form, df, REML=REML)
        }
    }, error=function(e){ message(e); NULL })
    if(is.null(mod)) return(mod)
    tryCatch({
        cvec <- as.numeric(colnames(coef(mod)[[1]])==testcoef)
        d <- cbind(coef(summary(mod)), p=NA_real_)
        if("Pr(>|t|)" %in% colnames(d)){
            d[,"p"] <- d[,"Pr(>|t|)"]
        }else{
            d[,"p"] <- NA_real_
        }
        d[which(cvec==1),"p"] <- last(contest(mod, cvec, ddf=ddf))
        co <- d
        list( sigma=sd(residuals(mod)),
              beta=co[,1],
              df.residual=df.residual(mod),
              SE=co[,2],
              stat=co[,3],
              pval=co[,"p"]
        )
    }, error=function(e){ message(e); NULL })
}

# formats a list of .fitlmer results into an eBayes compatible list and performs moderation
#' @importFrom limma eBayes
.mmEBayesWrapper <- function( fit.res, testcoef, trended=TRUE ){
    rl <- fit.res[!sapply(fit.res,is.null)]
    res <- list( coefficients=t(sapply(rl,FUN=function(x) x$beta )),
                 stdev.unscaled=t(sapply(rl,FUN=function(x) x$SE) ),
                 df.residual=rep(rl[[1]]$df.residual, length(rl)),
                 sigma=sapply(rl, FUN=function(x) x$sigma),
                 z=t(sapply(rl, FUN=function(x) x$stat)),
                 Amean=rowMeans(assay(sce)[names(rl),]),
                 PValue=t(sapply(rl, FUN=function(x) x$pval)) )
    res <- limma::eBayes(res, trend=trended, robust=TRUE)
    data.frame( row.names=names(rl), 
                beta=res$coefficients[,testcoef], 
                p_val.orig=res$PValue[,testcoef], 
                stat=res$z[,testcoef],
                p_val=res$p.value[,testcoef] )
}