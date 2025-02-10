#' bbhw: Bulk-based hypothesis weighing
#'
#' This is a method based on Independent Hypothesis Weighing (Ignatiadis et al.,
#' 2016) to increase the power of low-sample-size, per-celltype differential 
#' state analysis, using a larger dataset of bulk RNAseq.
#'
#' @param pbDEA A data.frame of pseudo-bulk DEA results, as for instance 
#'   produced by \code{\link{pbDS}} or \code{\link{mmDS}} (specifically this 
#'   should be a data.frame for one contrast, e.g. an element of `res$table` of 
#'   the output). This should contain the columns "gene", "cluster_id", "p_val"
#'   and, optionally "logFC".
#' @param bulkDEA A data.frame of bulk DEA results, with gene names as row.names
#'   and including the columns "p_val", "gene", "cluster_id" and, ideally,
#'   "logFC". Alternatively, a named vector of significance values. Not that 
#'   these samples should be independent form the single-cell samples on which
#'   `pbDEA` is based.
#' @param pb A pseudo-bulk SummarizedExperiment object as produced by 
#'   \code{\link{aggregateData}}.
#' @param method The method(s) for which to compiled the adjusted p-values.
#'   By default all methods are computed, but the recommended method is "PAShw"
#'   or, if `pb` is unavailable, `ihw.global`.
#' @param nbins The number of significance bins to use for the covariate (i.e. 
#'   prior). For `combIHW`, the effective number of bins will be doubled. We 
#'   recommend leaving this to the default values (which varies across methods).
#' @param alpha The nominal level for FDR control for \code{\link[IHW]{ihw}}.
#' @param ... Passed to \code{\link[IHW]{ihw}}.
#' @param verbose Logical; whether to print helpful information.
#' @param BPPARAM An optional BiocParallel BPPARAM object for multithreading.
#'
#' @return The `pbDEA` object including extra columns.
#' 
#' @author Pierre-Luc Germain
#' @references 
#' Ignatiadis, N., Klaus, B., Zaugg, J. et al. Data-driven hypothesis weighting 
#' increases detection power in genome-scale multiple testing. Nat Methods 13, 
#' 577–580 (2016). https://doi.org/10.1038/nmeth.3885
#' 
#' @details
#' This function contains the following methods:
#' * ihw.local and ihw.global simply bin the genes by their significance on the
#'    bulk RNA-seq, and perform standard independent hypothesis weighing (
#'    \code{\link[IHW]{ihw}} ) either locally (i.e. for each cell type 
#'    separately) or globally.
#' * Proportion-adjusted significance hypothesis weighting (PAShw and 
#'   PAShw.local) also rely on ihw, but first adjust the covariate according to
#'   whether the direction of the difference matches (assuming a logFC is 
#'   available), and then according to how much the cell type contributes to the
#'   bulk for each specific gene.
#' * combIHW is the standard global IHW, but the significance bins are further 
#'   split into low- and high- contributions of the celltype to the bulk, and 
#'   the bins are treated as nominal.
#' 
#' @export
#' @importFrom IHW ihw adj_pvalues
#' @importFrom SummarizedExperiment assays
#' @importFrom stats quantile
#' @importFrom dplyr bind_rows
#' @importFrom gtools logit inv.logit
bbhw <- function(pbDEA, bulkDEA, pb=NULL, 
                  method=c("PAShw","PAShw.local","combIHW","ihw.local",
                           "ihw.global"), nbins=NULL, alpha=0.1, ..., 
                 verbose=TRUE, BPPARAM=SerialParam(progress=verbose)){
  
  method <- match.arg(method, several.ok = TRUE)
  stopifnot(length(method)>0)
  stopifnot(is.null(nbins) || nbins>1L)
  stopifnot(alpha>0 && alpha<1)
  
  if(any(c("PAShw","PAShw.local","combIHW") %in% method)){
    if(is.null(pb)){
      stop("The selected adjustment method(s) require the `pb` argument to be ",
           "given. If you don't have this data, use `ihw.global` instead.")
    }
  }else{
    pb <- NULL
  }
  
  stopifnot(is.data.frame(pbDEA) &&
              all(c("p_val","gene","cluster_id") %in% colnames(pbDEA)))
  if(is.data.frame(bulkDEA)){
    if("logFC" %in% colnames(bulkDEA) && "logFC" %in% colnames(pbDEA)){
      bulkDEA <- setNames(sign(bulkDEA$logFC)*bulkDEA$PValue,
                          row.names(bulkDEA))
    }else{
      bulkDEA <- setNames(bulkDEA$PValue, row.names(bulkDEA))
    }
  }else{
    stopifnot(is.numeric(bulkDEA) && !is.null(names(bulkDEA)))
  }
  ig <- intersect(as.character(unique(pbDEA$gene)), names(bulkDEA))
  if(ig < 2000) stop("Too few genes in common for the procedure.")
  if(verbose){
    prop.common.genes <- length(ig)/max(length(bulkDEA),length(unique(pbDEA)))
    if(prop.common.genes<0.5){
      warning("`pbDEA` and `bulkDEA` have a lot of genes (", 
              round(100*(1-prop.common.genes)), "%) that are not in common.")
    }else if(prop.common.genes<1){
      message(round(100*(1-prop.common.genes)), "% genes in common")
    }
  }
    
  if(!is.null(pb)){
    if(inherits(pb, "SummarizedExpriment")){
      stopifnot(length(assays(pb))>1)
      # compute proportion of reads expected to be contributed to the bulk for
      # each gene/celltype
      rs <- vapply(assays(pb), FUN=rowSums, FUN.VALUE=numeric(nrow(pb)))
      props <- (1L+rs)/(rowSums(1L+rs))
    }else if(is.matrix(pb)){
      props <- pb
      rs <- rowSums(props)
      if(sum(rs>1,na.rm=TRUE)/length(rs) > 0.25){
        props <- (1L+props)/(rowSums(1L+props))
      }
    }
    if(!all(ig %in% row.names(props))){
      stop("All the genes that intersect between `pbDEA` and `bulkDEA` should ",
           "be present in `pb`.")
    }
    ct <- intersect(colnames(props), as.character(unique(pbDEA$cluster_id)))
    if(length(ct)==0)
      stop("The `cluster_id` column of `pbDEA` does not match the assays of `pb`.")
    all.cts <- union(colnames(props), unique(pbDEA$cluster_id))
    if(verbose && length(ct)<length(all.cts)){
      message("Only ",length(ct),"/",length(all.cts), " celltypes have data in",
              "both `pbDEA` and `pb`. The results will be restricted to those.")
    }
    props <- data.frame(gene=rep(row.names(props),ncol(props)),
                        cluster_id=rep(factor(colnames(props)),
                                       each=nrow(props)),
                        readProportion=as.numeric(as.matrix(props)))
    pbDEA <- merge(pbDEA, props, by=c("gene","cluster_id"))
  }
  
  pbDEA$bulkUnsigned <- pbDEA$bulk <- bulkDEA[as.character(pbDEA$gene)]
  if(any(pbDEA$bulk<0)){
    pbDEA$bulk[which(sign(pbDEA$logFC) != sign(pbDEA$bulk))] <- 0.6
  }
  pbDEA$bulk[is.na(pbDEA$bulk)] <- 0
  pbDEA$bulkUnsigned[is.na(pbDEA$bulkUnsigned)] <- 0
  pbDEA$bulk <- abs(pbDEA$bulk)
  pbDEA$bulkUnsigned <- abs(pbDEA$bulkUnsigned)
  
  if("ihw.global" %in% methods){
    if(verbose) message("Running ihw.global")
    a <- IHW::ihw(pbDEA$p_val, abs(pbDEA$bulkUnsigned), alpha=alpha)
    pbDEA$padj.ihw_glb <- IHW::adj_pvalues(a)
    if(verbose) .checkIhwRes(a)
  }
  
  if("PAShw" %in% methods){
    if(verbose) message("Running PAShw")
    pbDEA$padj.PAShw <- .PAShw(pbDEA, verbose=verbose, alpha=alpha, nfolds=4L,
                               nbins=nbins)
  }
  
  if("combIHW" %in% methods){
    if(verbose) message("Running combIHW")
    pbDEA$padj.combIHW <- .combIHW(pbDEA, verbose=verbose, alpha=alpha, 
                                   nbins=nbins)
  }
  
  if(any(grepl("local",method))){
    pbDEA <- dplyr::bind_rows(bplapply(split(pbDEA, pbDEA$cluster_id),
                                       BPPARAM=BPPARAM, FUN=function(pbDEA){
      if("ihw.local" %in% methods){
        a <- IHW::ihw(pbDEA$p_val, abs(pbDEA$bulkUnsigned), alpha=alpha)
        pbDEA$padj.ihw_loc <- IHW::adj_pvalues(a)
      }
      
      if("PAShw.local" %in% methods){
       pbDEA$padj.PAShw <- .PAShw(pbDEA, verbose=FALSE, alpha=alpha, nfolds=3L,
                                  nbins=nbins)
      }
      pbDEA
    }))
  }

  pbDEA
}

.getQBreaks <- function(x, nbins){
  c(0,quantile(x, seq_len(nbins-1)/(nbins)),1)
}

.checkIhwRes <- function(a, msg=TRUE){
  cc <- suppressWarnings(cor(a@weights))
  cc <- cc[lower.tri(cc)]
  if(sum(cc>0.5, na.rm=TRUE)/length(cc) > 0.5) return(TRUE)
  if(!msg) return(FALSE)
  message("The prior (i.e. bulk data) leads to unstable weights and might ",
          "not be informative. If this is the case, you'd better rely on ",
          " standard local FDR adjustment.")
}

.PAShw <- function(pbDEA, verbose=TRUE, alpha=0.1, nfolds=3L, nbins=NULL, ...){
  pbDEA$combCov <- 
    gtools::inv.logit(gtools::logit(pbDEA$bulk)*sqrt(pbDEA$readProportion))
  pbDEA$combCov[is.na(pbDEA$combCov)] <- 1
  if(is.null(nbins))
    nbins <- ifelse(nrow(pbDEA)>10000,9,ifelse(nrow(pbDEA)>7000,7,5))-1L
  breaks <- c(0,quantile(pbDEA$combCov,
                       seq_len(nbins)[-ceiling((nbins+1L)/2)]/(nbins+1L),
                       na.rm=TRUE), 1)
  pbDEA$combCov <- cut(pbDEA$combCov, unique(breaks))
  pbDEA$combCov[which(is.na(pbDEA$combCov) | pbDEA$bulk>1)] <- 
    levels(pbDEA$combCov)[ceiling(length(levels(pbDEA$combCov))/2)]
  a <- IHW::ihw(pbDEA$p_val, droplevels(pbDEA$combCov), alpha=alpha,
                nfolds=nfolds, ...)
  .checkIhwRes(a, verbose)
  IHW::adj_pvalues(a)
}

.combIHW <- function(pbDEA, verbose=TRUE, alpha=0.1, nbins=NULL, ...){
  if(is.null(nbins)) nbins <- ifelse(nrow(pbDEA)>10000,5,3)
  breaks <- .getQBreaks(pbDEA$bulk, nbins)
  pbDEA$tmpROWindex <- seq_len(nrow(pbDEA))
  pbDEA$sigBin <- cut(pbDEA$bulk, breaks, include.lowest=TRUE)
  pbDEA <- dplyr::bind_rows(lapply(split(pbDEA, pbDEA$sigBin), FUN=\(x){
    if(min(table(highProp <- x$readProportion>p.th))>1000){
      x$propBin <- highProp
    }else if(nrow(x)>2000){
      x$propBin <- x$readProportion > median(x$readProportion)
    }else{
      x$propBin <- FALSE
    }
    x
  }))
  pbDEA$covBin <- as.factor(paste(pbDEA$sigBin, pbDEA$propBin))
  pbDEA <- pbDEA[order(pbDEA$tmpROWindex),]
  a <- IHW::ihw(pbDEA$p_val, pbDEA$covBin, alpha=alpha,
                covariate_type="nominal", nfolds=3L, ...)
  .checkIhwRes(a, verbose)
  IHW::adj_pvalues(a)
}
