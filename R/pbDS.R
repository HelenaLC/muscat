#' @rdname pbDS
#' @title pseudobulk DS analysis
#'
#' @description \code{pbDS} tests for DS after aggregating single-cell 
#'   measurements to pseudobulk data, by applying bulk RNA-seq DE methods, 
#'   such as \code{edgeR}, \code{DESeq2} and \code{limma}.
#'
#' @param pb a \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   containing pseudobulks as returned by \code{\link{aggregateData}}.
#' @param method a character string.
#' @param design For methods \code{"edegR"} and \code{"limma"}, a design matrix 
#'   with row & column names(!) created with \code{\link[stats]{model.matrix}}; 
#'   For \code{"DESeq2"}, a formula with variables in \code{colData(pb)}.
#'   Defaults to \code{~ group_id} or the corresponding \code{model.matrix}.
#' @param contrast a matrix of contrasts to test for
#'   created with \code{\link[limma]{makeContrasts}}.
#' @param coef passed to \code{\link[edgeR]{glmQLFTest}},
#'   \code{\link[limma]{contrasts.fit}}, \code{\link[DESeq2]{results}}
#'   for \code{method = "edgeR", "limma-x", "DESeq2"}, respectively.
#' @param min_cells a numeric. Specifies the minimum number of cells in a given 
#'   cluster-sample required to consider the sample for differential testing.
#' @param verbose logical. Should information on progress be reported?
#'
#' @return a list containing \itemize{
#' \item a data.frame with differential testing results,
#' \item a \code{\link[edgeR]{DGEList}} object of length nb.-clusters, and
#' \item the \code{design} matrix, and \code{contrast} or \code{coef} used.}
#'
#' @examples
#' # simulate 5 clusters, 20% of DE genes
#' data(sce)
#'     
#' # compute pseudobulk sum-counts & run DS analysis
#' pb <- aggregateData(sce)
#' res <- pbDS(pb, method = "limma-trend")
#'
#' names(res)
#' names(res$table)
#' head(res$table$`stim-ctrl`$`B cells`)
#' 
#' # count nb. of DE genes by cluster
#' vapply(res$table$`stim-ctrl`, function(u) 
#'   sum(u$p_adj.loc < 0.05), numeric(1))
#' 
#' # get top 5 hits for ea. cluster w/ abs(logFC) > 1
#' library(dplyr)
#' lapply(res$table$`stim-ctrl`, function(u)
#'   filter(u, abs(logFC) > 1) %>% 
#'     arrange(p_adj.loc) %>% 
#'     slice(seq_len(5)))
#'
#' @author Helena L Crowell & Mark D Robinson
#'
#' @references 
#' Crowell, HL, Soneson, C, Germain, P-L, Calini, D, 
#' Collin, L, Raposo, C, Malhotra, D & Robinson, MD: 
#' On the discovery of population-specific state transitions from 
#' multi-sample multi-condition single-cell RNA sequencing data. 
#' \emph{bioRxiv} \strong{713412} (2018). 
#' doi: \url{https://doi.org/10.1101/713412}
#'
#' @importFrom edgeR filterByExpr
#' @importFrom dplyr last rename
#' @importFrom limma makeContrasts
#' @importFrom Matrix qr
#' @importFrom methods is
#' @importFrom scater isOutlier
#' @importFrom stats model.matrix
#' @importFrom SummarizedExperiment colData
#' @export

pbDS <- function(pb, 
    method = c("edgeR", "DESeq2", "limma-trend", "limma-voom"),
    design = NULL, coef  = NULL, contrast = NULL, min_cells = 10, 
    filter = c("both", "genes", "samples", "none"), verbose = TRUE) {
    
    # check validity of input arguments
    method <- match.arg(method)
    filter <- match.arg(filter)
    .check_pbs(pb, check_by = TRUE)
    .check_args_pbDS(as.list(environment()))
    
    if (is.null(design)) {
        formula <- ~ group_id
        cd <- as.data.frame(colData(pb))
        design <- model.matrix(formula, cd)
        colnames(design) <- levels(pb$group_id)
    }
    
    if (is.null(coef) & is.null(contrast)) {
        c <- colnames(design)[ncol(design)]
        contrast <- makeContrasts(contrasts = c, levels = design)
    }

    # ct: type of comparison - "contrast" or "coef"
    # cs: named list of 'coef's or 'contrast's
    if (!is.null(contrast)) {
        coef <- NULL
        names(cs) <- cs <- colnames(contrast)
    } else if (!is.null(coef)) {
        cs <- vapply(coef, function(i)
            paste(colnames(design)[i], collapse = "-"),
            character(1))
        names(cs) <- names(coef) <- cs
    }
    ct <- ifelse(is.null(coef), "contrast", "coef")
    
    # compute cluster-sample counts
    n_cells <- metadata(pb)$n_cells
    names(kids) <- kids <- assayNames(pb)
    
    if (!is.function(method)) {
        fun <- switch(method,
            "DESeq2" = .DESeq2,
            "edgeR" = .edgeR, 
            "limma-trend" = .limma_trend, 
            "limma-voom" = .limma_voom)
    } else {
        fun_call <- 1
    }
    fun_args <- names(as.list(args(fun)))
    fun_args <- fun_args[-length(fun_args)]
    
    # for ea. cluster, run DEA
    res <- lapply(kids, function (k) {
        if (verbose) cat(k, "..", sep = "")
        rmv <- n_cells[k, ] < min_cells
        d <- design[colnames(y <- pb[ , !rmv]), , drop = FALSE]
        if (filter %in% c("samples", "both")) {
            ls <- colSums(assay(y, k))
            ol <- isOutlier(ls, log = TRUE, type = "lower", nmads = 3)
            d <- d[colnames(y <- y[, !ol]), , drop = FALSE]
        }
        if (any(tabulate(y$group_id) < 2) 
            || qr(d)$rank == nrow(d) 
            || qr(d)$rank < ncol(d)) 
            return(NULL)
        if (filter %in% c("genes", "both"))
            y <- y[filterByExpr(assay(y, k)), ]
        args <- list(
            x = y, k = k, design = d, coef = coef, 
            contrast = contrast, ct = ct, cs = cs)
        args <- args[intersect(names(args), fun_args)]
        do.call(fun, args)
    })
    
    # remove empty clusters
    rmv <- vapply(res, is.null, logical(1))
    res <- res[!rmv]
    kids <- kids[names(res)]
    
    # reorganize & do global p-value adjustment
    names(i) <- i <- c("table", "data", "fit")
    res <- lapply(i, map, .x = res)
    res$table <- .p_adj_global(res$table)
    return(res)
}
