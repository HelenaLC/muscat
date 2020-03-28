#' @name prepSim
#' 
#' @title SCE preparation for \code{\link{simData}}
#' 
#' @description \code{prepSim} prepares an input SCE for simulation 
#'   with \code{muscat}'s \code{\link{simData}} function by 
#' \enumerate{
#'   \item{basic filtering of genes and cells}
#'   \item{(optional) filtering of subpopulation-sample instances}
#'   \item{estimation of cell (library sizes) and gene parameters 
#'   (dispersions and sample-specific means), respectively.}
#' }
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param min_count,min_cells used for filtering of genes; only genes with 
#'   a count > \code{min_count} in >= \code{min_cells} will be retained.
#' @param min_genes used for filtering cells;
#'   only cells with a count > 0 in >= \code{min_genes} will be retained.
#' @param min_size used for filtering subpopulation-sample combinations;
#'   only instances with >= \code{min_size} cells will be retained.
#'   Specifying \code{min_size = NULL} skips this step.
#' @param group_keep character string; if \code{nlevels(x$group_id) > 1},
#'   specifies which group of samples to keep (see details). The default
#'   NULL retains samples from \code{levels(x$group_id)[1]}; otherwise, 
#'   if `colData(x)$group_id` is not specified, all samples will be kept.
#' @param verbose logical; should information on progress be reported?
#' 
#' @details For each gene \eqn{g}, \code{prepSim} fits a model to estimate 
#'   sample-specific means \eqn{\beta_g^s}, for each sample \eqn{s}, 
#'   and dispersion parameters \eqn{\phi_g} using \code{edgeR}'s 
#'   \code{\link[edgeR]{estimateDisp}} function with default parameters. 
#'   Thus, the reference count data is modeled as NB distributed: 
#'   \deqn{Y_{gc} \sim NB(\mu_{gc}, \phi_g)}
#'   for gene \eqn{g} and cell \eqn{c}, where the mean 
#'   \eqn{\mu_{gc} = \exp(\beta_{g}^{s(c)}) \cdot \lambda_c}. Here, 
#'   \eqn{\beta_{g}^{s(c)}} is the relative abundance of gene \eqn{g} 
#'   in sample \eqn{s(c)}, \eqn{\lambda_c} is the library size 
#'   (total number of counts), and \eqn{\phi_g} is the dispersion.
#'   
#' @return a \code{\link[SingleCellExperiment]{SingleCellExperiment}} 
#'   containing, for each cell, library size (\code{colData(x)$offset})
#'   and, for each gene, dispersion and sample-specific mean estimates 
#'   (\code{rowData(x)$dispersion} and \code{$beta.sample_id}, respectively).
#' 
#' @examples
#' data(sce)
#' library(SingleCellExperiment)
#' 
#' ref <- prepSim(sce)
#' 
#' # nb. of genes/cells before vs. after
#' ns <- cbind(before = dim(sce), after = dim(ref)) 
#' rownames(ns) <- c("#genes", "#cells"); ns
#' 
#' head(rowData(ref)) # gene parameters
#' head(colData(ref)) # cell parameters
#' 
#' @author Helena L Crowell
#' 
#' @references 
#' Crowell, HL, Soneson, C, Germain, P-L, Calini, D, 
#' Collin, L, Raposo, C, Malhotra, D & Robinson, MD: 
#' On the discovery of population-specific state transitions from 
#' multi-sample multi-condition single-cell RNA sequencing data. 
#' \emph{bioRxiv} \strong{713412} (2018). 
#' doi: \url{https://doi.org/10.1101/713412}
#' 
#' @importFrom edgeR DGEList estimateDisp glmFit
#' @importFrom Matrix colSums rowSums
#' @importFrom SingleCellExperiment SingleCellExperiment counts
#' @importFrom SummarizedExperiment colData rowData<-
#' @importFrom stats model.matrix
#' @export

prepSim <- function(x, 
    min_count = 1, min_cells = 10, 
    min_genes = 100, min_size = 100, 
    group_keep = NULL, verbose = TRUE) {
    
    .check_sce(x, req_group = FALSE)
    stopifnot(is.numeric(min_count), 
        is.numeric(min_cells), is.numeric(min_genes), 
        is.null(min_size) || is.numeric(min_size),
        is.logical(verbose), length(verbose) == 1)
    
    n_cells0 <- ncol(x)
    if (is.null(group_keep)) {
        if ("group_id" %in% colnames(colData(x))) {
            group_keep <- levels(x$group_id)[1]
            if (verbose) message(sprintf(paste(
                "Argument `group_keep` unspecified;",
                "defaulting to retaining %s-group samples."),
                dQuote(group_keep)))
            cells_keep <- x$group_id == group_keep
        } else {
            cells_keep <- seq_len(ncol(x))
        }
    } else {
        stopifnot(is.character(group_keep), 
            group_keep %in% levels(x$group_id))
        cells_keep <- x$group_id %in% group_keep
    }
    x <- .update_sce(x[, cells_keep])
    
    # keep genes w/ count > `min_count` in at least `min_cells`;
    # keep cells w/ at least `min_genes` detected genes
    if (verbose) message("Filtering...")
    genes_keep <- rowSums(counts(x) > min_count) >= min_cells
    cells_keep <- colSums(counts(x) > 0) >= min_genes
    if (verbose) message(sprintf(
        "- %s/%s genes and %s/%s cells retained.",
        sum(genes_keep), nrow(x), sum(cells_keep), n_cells0))
    x <- x[genes_keep, cells_keep, drop = FALSE]
    
    # keep cluster-samples w/ at least 'min_size' cells
    if (!is.null(min_size)) {
        n_cells <- table(x$cluster_id, x$sample_id)
        n_cells <- .filter_matrix(n_cells, n = min_size)
        if (ncol(n_cells) == 1)
            stop("Current 'min_size' retains only 1 sample,\nbut",
                " mean-dispersion estimation requires at least 2.")
        if (verbose) message(sprintf(
            "- %s/%s subpopulations and %s/%s samples retained.",
            nrow(n_cells), nlevels(x$cluster_id), 
            ncol(n_cells), nlevels(x$sample_id)))
        x <- .filter_sce(x, rownames(n_cells), colnames(n_cells))
    }
    
    # estimate gene/cell parameters 
    if (verbose) 
        message("Estimating gene and cell parameters...")
    y <- DGEList(counts(x))
    mm <- model.matrix(~ 0 + x$sample_id)
    y <- estimateDisp(y, mm)
    y <- glmFit(y, prior.count = 0)
    
    # update row- & colData
    x$offset <- c(y$offset)
    rowData(x)$dispersion <- y$dispersion
    betas <- paste0("beta.", levels(x$sample_id))
    colnames(y$coefficients) <- betas
    rowData(x)[, betas] <- y$coefficients
    
    # return SCE
    return(x)
}