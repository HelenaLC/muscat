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
#'   \code{\link{aggregateData}}. Alternatively, a matrix with cell types as 
#'   columns and gene as rows, giving the read counts or proportion of 
#'   contribution for each gene. If neither is given, the only available methods
#'   are "ihw.local" and "ihw.global".
#' @param method The method(s) for which to compiled the adjusted p-values (see
#'   details below for options). By default the local and global versions of the
#'    top methods are computed. If `pb` is unavailable, `ihw.global` is the 
#'    recommended method. Use `method="all"` to compute adjusted p-values using
#'    all methods.
#' @param nbins The number of significance bins to use for the covariate (i.e. 
#'   prior). For `combIHW`, the effective number of bins will be doubled. We 
#'   recommend leaving this to the default values (which varies across methods).
#' @param nfolds The number of cross-validation folds, passed to 
#'   \code{\link[IHW]{ihw}}. If null, will use appropriate defaults for 
#'   different methods based on the number of hypotheses.
#' @param alpha The nominal level for FDR control for \code{\link[IHW]{ihw}}.
#' @param ... Passed to \code{\link[IHW]{ihw}}.
#' @param BPPARAM An optional BiocParallel BPPARAM object for multithreading.
#' @param verbose Logical; whether to print helpful information.
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
#' \itemize{
#' \item{**ihw** : the IHW procedure is applied with default settings, 
#'  using the bulk p-value as covariate.}
#' \item{**combIHW** : first, whenever the direction of the change is 
#'  different between bulk and pseudobulk datasets for a gene in a given cell 
#'  type, we set the bulk p-value to 0.7 for that cell type. Then, we divide the 
#'  bulk p-values into five quantile bins. We further divide each bin into two 
#'  depending on the proportion of the bulk reads for that gene that is 
#'  contributed by the given cell type. Then we apply IHW (globally) on this 
#'  covariate, using the bins as nominal. Note that this method is slightly 
#'  inferior to PASW and PABW, which are instead recommended.}
#' \item{**PASW** (Proportion-Adjusted Significance Weighing): first, whenever 
#'  the direction of the change is different between bulk and pseudobulk 
#'  datasets for a gene in a given cell type, we increase the bulk p-value to 
#'  0.7 for that cell type (if it was below). Then, we adjust the covariate 
#'  (i.e. bulk p-value) based on the proportion of the bulk reads for that gene 
#'  that is contributed by the given cell type, using 
#'  `inv.logit( logit(p) * sqrt(c) )` where `p` and `c` are respectively 
#'  the bulk p-value and the proportion of bulk reads contributed by the cell 
#'  type). We then split this covariate into quantile bins and apply IHW.
#'  If unspecified (recommended), the number of bins will be determined 
#'  (somewhere between 6 and 10) based on the number of hypotheses.}
#' \item{**PABW** (Proportion-Adjusted Bin-Wise correction): the covariate bins
#'  are prepared in the same fashion as for PASW. However, instead of using IHW,
#'  we simply compute FDR separately for each bin.}
#' }
#' 
#' If you have information about the contribution of each cell type to each 
#' bulk gene (i.e. the `pb` argument), the recommended methods are PASW and 
#' PABW, which give similar results (PABW is much faster and deterministic). If
#' you do not have such information, use the IHW method.
#' 
#' Each method exists in two flavors: a local one, which is applied for each 
#' cell type separately, and a global one, which is applied once across all cell
#' types. Unfortunately, which of the two is preferable seems to depends on the
#' context.
#' 
#' The `method` argument should indicate both the method and whether it should 
#' be applied locally or globally, e.g. `method="PABW.global"`.
#' 
#' @export
#' @importFrom IHW ihw adj_pvalues
#' @importFrom SummarizedExperiment assays
#' @importFrom stats quantile
#' @importFrom dplyr bind_rows
#' @importFrom gtools logit inv.logit
bbhw <- function(pbDEA, bulkDEA, pb=NULL, method=c("PASW.local","PASW.global"),
                 nbins=NULL, alpha=0.1, nfolds=NULL, ...,
                 BPPARAM=SerialParam(progress=verbose), verbose=TRUE){
  allmethods <- c("PASW.global","PASW.local","PABW.local","PABW.global",
                  "combIHW.global","combIHW.local","ihw.local", "ihw.global")
  method <- match.arg(method, choices=c("all",allmethods), several.ok=TRUE)
  if(any("all" %in% method)) method <- allmethods
  stopifnot(is.null(nbins) || nbins>1L)
  stopifnot(alpha>0 && alpha<1)
  
  if(any(grepl("PASW|combIHW|PABW",method))){
    if(is.null(pb)){
      stop("The selected adjustment method(s) require the `pb` argument to be ",
           "given. If you don't have this data, use `ihw.global` instead.")
    }
  }else{
    pb <- NULL
  }
  
  stopifnot(is.data.frame(pbDEA))
  pbDEA <- .homogenizeDEA4bbhw(pbDEA)
  stopifnot(all(c("p_val","gene","cluster_id") %in% colnames(pbDEA)))
  if(is.data.frame(bulkDEA)){
    bulkDEA <- .homogenizeDEA4bbhw(bulkDEA)
    stopifnot("p_val" %in% colnames(bulkDEA))
    if("logFC" %in% colnames(bulkDEA) && "logFC" %in% colnames(pbDEA)){
      bulkDEA <- setNames(sign(bulkDEA$logFC)*bulkDEA$p_val,
                          row.names(bulkDEA))
    }else{
      bulkDEA <- setNames(bulkDEA$p_val, row.names(bulkDEA))
    }
  }else{
    stopifnot(is.numeric(bulkDEA) && !is.null(names(bulkDEA)))
  }
  ig <- intersect(as.character(unique(pbDEA$gene)), names(bulkDEA))
  if(length(ig) < 2000) stop("Too few genes in common for the procedure.")
  if(verbose){
    prop.common.genes <- length(ig)/max(length(bulkDEA),
                                        length(unique(pbDEA$gene)))
    if(prop.common.genes<0.5){
      warning("`pbDEA` and `bulkDEA` have a lot of genes (", 
              round(100*(1-prop.common.genes)), "%) that are not in common.")
    }else if(prop.common.genes<1){
      message(round(100*(prop.common.genes)), "% genes in common")
    }
  }
    
  if(!is.null(pb)){
    if(inherits(pb, "SummarizedExpriment") || 
       inherits(pb, "SingleCellExperiment")){
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
    }else{
      stop("Invalid `pb` input.")
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
    pbDEA <- merge(pbDEA[,setdiff(colnames(pbDEA), "readProportion")], props,
                   by=c("gene","cluster_id"))
  }
  
  pbDEA$bulkUnsigned <- pbDEA$bulk <- bulkDEA[as.character(pbDEA$gene)]
  if(any(pbDEA$bulk<0)){
    w <- which(sign(pbDEA$logFC) != sign(pbDEA$bulk))
    pbDEA$bulk[w] <- pmax(pbDEA$bulk[w],0.7)
  }
  pbDEA$bulk[is.na(pbDEA$bulk)] <- 0
  pbDEA$bulkUnsigned[is.na(pbDEA$bulkUnsigned)] <- 0
  pbDEA$bulk <- abs(pbDEA$bulk)
  pbDEA$bulkUnsigned <- abs(pbDEA$bulkUnsigned)
  
  if("ihw.global" %in% method){
    if(verbose && length(method)>1) message("Running ihw.global")
    a <- IHW::ihw(pbDEA$p_val, abs(pbDEA$bulkUnsigned), alpha=alpha,
                  nfolds=ifelse(is.null(nfolds),5L,nfolds))
    pbDEA$padj.ihw_glb <- IHW::adj_pvalues(a)
    if(verbose) .checkIhwRes(a)
  }
  
  if("PASW.global" %in% method){
    if(verbose && length(method)>1) message("Running PASW")
    pbDEA$padj.PASW_glb <- .PASW(pbDEA, verbose=verbose, alpha=alpha, 
                                 nfolds=nfolds, nbins=nbins)
  }
  
  if("combIHW.global" %in% method){
    if(verbose && length(method)>1) message("Running combIHW")
    pbDEA$padj.combIHW_glb <- .combIHW(pbDEA, verbose=verbose, alpha=alpha, 
                                       nfolds=nfolds, nbins=nbins)
  }
  
  if("PABW.global" %in% method){
    if(verbose && length(method)>1) message("Running PASW.global")
    pbDEA$padj.PABW_glb <- .PABW(pbDEA, verbose=FALSE, nbins=nbins)
  }
  
  
  
  if(any(grepl("local",method))){
    if(verbose && length(method)>1) message("Running local methods")
    pbDEA <- dplyr::bind_rows(bplapply(split(pbDEA, pbDEA$cluster_id),
                                       BPPARAM=BPPARAM, FUN=function(pbDEA){
      if("ihw.local" %in% method){
        a <- IHW::ihw(pbDEA$p_val, abs(pbDEA$bulkUnsigned), alpha=alpha, 
                      nfolds=ifelse(is.null(nfolds),4L,nfolds))
        pbDEA$padj.ihw_loc <- IHW::adj_pvalues(a)
      }
      
      if("PASW.local" %in% method){
       pbDEA$padj.PASW_loc <- .PASW(pbDEA, verbose=FALSE, alpha=alpha,
                                    nfolds=nfolds, nbins=nbins)
      }
                                         
      if("PABW.local" %in% method){
        pbDEA$padj.PABW_loc <- .PABW(pbDEA, verbose=FALSE, nbins=nbins)
      }
                                         
      if("combIHW.local" %in% method){
        pbDEA$padj.combIHW_loc <- .combIHW(pbDEA, verbose=FALSE, alpha=alpha, 
                                           nfolds=nfolds, nbins=nbins)
      }
      pbDEA
    }))
  }

  pbDEA
}

.getQBreaks <- function(x, nbins){
  c(0,quantile(x, seq_len(nbins-1)/(nbins)),1)
}

.getNfolds <- function(x, nfolds){
  if(!is.null(nfolds)) return(nfolds)
  ifelse(x>1e5, 5L, ifelse(x>1e4, 4L, 3L))
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

.PASW <- function(pbDEA, verbose=TRUE, alpha=0.1, nfolds=NULL, nbins=NULL, ...){
  nfolds <- .getNfolds(nrow(pbDEA), nfolds)
  pbDEA <- .getAdjustedBins(pbDEA, nbins=nbins)
  a <- IHW::ihw(pbDEA$p_val, droplevels(pbDEA$combCov), alpha=alpha,
                nfolds=nfolds, ...)
  .checkIhwRes(a, verbose)
  IHW::adj_pvalues(a)
}

.PABW <- function(pbDEA, verbose=TRUE, nbins=NULL){
  pbDEA <- .getAdjustedBins(pbDEA, nbins=nbins)
  pbDEA$tmpROWindex <- seq_len(nrow(pbDEA))
  pbDEA <- dplyr::bind_rows(lapply(split(pbDEA, pbDEA$combCov), FUN=\(x){
    x$FDR <- p.adjust(x$p_val, method="fdr")
    x
  }))
  pbDEA[order(pbDEA$tmpROWindex),"FDR"]
}

.getAdjustedBins <- function(pbDEA, nbins){
  pbDEA$combCov <- 
    gtools::inv.logit(gtools::logit(pbDEA$bulk)*sqrt(pbDEA$readProportion))
  pbDEA$combCov[is.na(pbDEA$combCov)] <- 1
  if(is.null(nbins))
    nbins <- ifelse(nrow(pbDEA)>10000,11,ifelse(nrow(pbDEA)>7000,9,7))-1L
  breaks <- c(0,quantile(pbDEA$combCov,
                         seq_len(nbins)[-ceiling((nbins+1L)/2)]/(nbins+1L),
                         na.rm=TRUE), 1)
  pbDEA$combCov <- cut(pbDEA$combCov, unique(breaks))
  pbDEA$combCov[which(is.na(pbDEA$combCov) | pbDEA$bulk>1)] <- 
    levels(pbDEA$combCov)[ceiling(length(levels(pbDEA$combCov))/2)]
  pbDEA  
}

.combIHW <- function(pbDEA, verbose=TRUE, alpha=0.1, nbins=NULL, nfolds=NULL, 
                     ...){
  nfolds <- .getNfolds(nrow(pbDEA), nfolds)
  if(is.null(nbins)) nbins <- ifelse(nrow(pbDEA)>10000,5,4)
  breaks <- .getQBreaks(pbDEA$bulk, nbins)
  pbDEA$tmpROWindex <- seq_len(nrow(pbDEA))
  pbDEA$sigBin <- cut(pbDEA$bulk, unique(breaks), include.lowest=TRUE)
  pbDEA <- dplyr::bind_rows(lapply(split(pbDEA, pbDEA$sigBin), FUN=\(x){
    if(nrow(x)<=2000){
      x$propBin <- FALSE
    }else{
      q <- quantile(x$readProportion, ifelse(nrow(x)>4000,0.75,0.5))
      x$propBin <- x$readProportion >= as.numeric(q)[1]
    }
    x
  }))
  pbDEA$covBin <- as.factor(paste(pbDEA$sigBin, pbDEA$propBin))
  pbDEA <- pbDEA[order(pbDEA$tmpROWindex),]
  a <- IHW::ihw(pbDEA$p_val, pbDEA$covBin, alpha=alpha,
                covariate_type="nominal", nfolds=nfolds, ...)
  .checkIhwRes(a, verbose)
  IHW::adj_pvalues(a)
}

.homogenizeDEA4bbhw <- function(x){
  if(is.null(x$p_val) && !is.null(x$PValue)) x$p_val <- x$PValue
  if(is.null(x$cluster_id) && !is.null(x$celltype)) x$cluster_id <- x$celltype
  x
}