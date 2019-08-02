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
#' The original data has been minimally filtered to remove undetected genes, 
#' cell multiplets, and unassigned cells. 
#' Assay \code{logcounts} corresponds to log-normalized values 
#' obtained from \code{\link[scater]{normalize}} with default parameters.
#'   
#' The complete raw, gene, and cell metadata 
#' is available through the NCBI GEO accession number 
#' [GSE96583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583).
#' 
#' @return a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' 
#' @references 
#' Kang et al. (2018). 
#' Multiplexed droplet single-cell RNA-sequencing 
#' using natural genetic variation. 
#' \emph{Nature Biotechnology}, 
#' \bold{36}(1): 89-94. DOI: 10.1038/nbt.4042.
#'  
#' @author Helena L Crowell

NULL