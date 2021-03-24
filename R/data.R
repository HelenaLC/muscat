#' @rdname data
#' @name data
#' @aliases data sce
#' 
#' @title Example datasets
#' 
#' @description 
#' A \code{\link[SingleCellExperiment]{SingleCellExperiment}} containing 
#' 10x droplet-based scRNA-seq PBCM data from 8 Lupus patients befor and after 
#' 6h-treatment with INF-beta (16 samples in total).
#' 
#' The original data has been filtered to
#' \itemize{
#' \item{remove unassigned cells & cell multiplets}
#' \item{retain only 4 out of 8 samples per experimental group}
#' \item{retain only 5 out of 8 subpopulations (clusters)}
#' \item{retain genes with a count > 1 in > 50 cells}
#' \item{retain cells with > 200 detected genes}
#' \item{retain at most 100 cells per cluster-sample instance}
#' }
#' 
#' Assay \code{logcounts} corresponds to log-normalized values 
#' obtained from \code{\link[scater]{logNormCounts}} with default parameters.
#'   
#' The original measurement data, as well as gene and cell metadata 
#' is available through the NCBI GEO accession number GSE96583;
#' code to reproduce this example dataset from the original data 
#' is provided in the examples section.
#' 
#' @return a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' 
#' @examples
#' \donttest{
#' # set random seed for cell sampling
#' set.seed(2929)
#' 
#' # load data
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' sce <- eh[["EH2259"]]
#' 
#' # drop unassigned cells & multiplets
#' sce <- sce[, !is.na(sce$cell)]
#' sce <- sce[, sce$multiplets == "singlet"]
#' 
#' # keep 4 samples per group
#' sce$id <- paste0(sce$stim, sce$ind)
#' inds <- sample(sce$ind, 4)
#' ids <- paste0(levels(sce$stim), rep(inds, each = 2))
#' sce <- sce[, sce$id %in% ids]
#' 
#' # keep 5 clusters
#' kids <- c("B cells", "CD4 T cells", "CD8 T cells",
#'     "CD14+ Monocytes", "FCGR3A+ Monocytes")
#' sce <- sce[, sce$cell %in% kids]
#' sce$cell <- droplevels(sce$cell)
#' 
#' # basic filtering on  genes & cells
#' gs <- rowSums(counts(sce) > 1) > 50
#' cs <- colSums(counts(sce) > 0) > 200
#' sce <- sce[gs, cs]
#' 
#' # sample max. 100 cells per cluster-sample
#' cs_by_ks <- split(colnames(sce), list(sce$cell, sce$id))
#' cs <- sapply(cs_by_ks, function(u) 
#'     sample(u, min(length(u), 100)))
#' sce <- sce[, unlist(cs)]
#' 
#' # compute logcounts
#' library(scater)
#' sce <- computeLibraryFactors(sce)
#' sce <- logNormCounts(sce)
#' 
#' # re-format for 'muscat'
#' sce <- prepSCE(sce, 
#'     cluster_id = "cell", 
#'     sample_id = "id", 
#'     group_id = "stim", 
#'     drop = TRUE)
#' } 
#' 
#' @references 
#' Kang et al. (2018). Multiplexed droplet single-cell RNA-sequencing 
#' using natural genetic variation. \emph{Nature Biotechnology},
#' \bold{36}(1): 89-94. DOI: 10.1038/nbt.4042.
#'  
#' @author Helena L Crowell

NULL
