#' bbhw: Bulk-based hypothesis weighing
#'
#' This is a method to increase the power of low-sample-size, per-celltype 
#' differential state analysis by using a larger dataset of bulk RNAseq. In at 
#' nutshell, it uses bulk data to create a covariate according to which the 
#' hypotheses are grouped, and then uses this grouping either either on via 
#' independent hypothesis weighing or grouped Benjamini-Hochberg correction to
#' increase power.
#'
#' @param pbDEA A data.frame of pseudo-bulk DEA results, as for instance 
#'   produced by \code{\link{pbDS}} or \code{\link{mmDS}} (specifically this 
#'   should be a data.frame for one contrast, e.g. an element of `res$table` of 
#'   the output). This should contain the columns "gene", "cluster_id", "p_val"
#'   and, optionally "logFC".
#' @param bulkDEA A data.frame of bulk DEA results, with gene names as row.names
#'   and including the column "p_val" and, ideally, "logFC". Alternatively, a 
#'   named vector of significance values. Not that these samples should be 
#'   independent form the single-cell samples on which `pbDEA` is based.
#' @param pb A pseudo-bulk SummarizedExperiment object as produced by 
#'   \code{\link{aggregateData}}. Alternatively, a matrix with cell types as 
#'   columns and gene as rows, giving the read counts or proportion of 
#'   contribution for each gene. If neither is given, the only available methods
#'   are "ihw.local" and "ihw.global".
#' @param bin.method The method for creating the bulk-based bins. Either "PAS" 
#'   (recommended and default), "combined", or "sig". Note that only 
#'   method="sig" is available if `pb` is not provided.
#' @param correction.method Determines with which method the bins are used. 
#'   See details for the different options. We recommend "gBH.LSL". 
#' @param local Logical; whether to apply the adjustment locally, i.e. 
#'   separately for each cell type (default TRUE).
#' @param useSign Logical; whether to discount bulk p-values for which the 
#'   change is in the opposite direction as in the given celltype.
#' @param nbins The number of significance bins to use for the covariate (i.e. 
#'   prior). For `combIHW`, the effective number of bins will be doubled (to 
#'   accommodate proportion-based bins). If omitted, a decent number of bins 
#'   will be set based on the method and number of hypotheses.
#' @param nfolds The number of cross-validation folds, passed to 
#'   \code{\link[IHW]{ihw}}. If null, will use appropriate defaults based on 
#'   the number of hypotheses per bin.
#' @param NAsep Logical; whether to put NA bulk p-values into their own bin 
#'   (assuming there is a sufficient number of them). Otherwise, NA values will 
#'   be set to 0.5. In practice there is often an enrichment for small p-values
#'   in genes that are undetectable at the bulk level, so we recommend setting 
#'   this to TRUE (default).
#' @param alpha The nominal level for FDR control for \code{\link[IHW]{ihw}}.
#' @param ... Passed to \code{\link[IHW]{ihw}}.
#' @param BPPARAM An optional BiocParallel BPPARAM object for multithreading.
#' @param verbose Logical; whether to print helpful information.
#'
#' @return The `pbDEA` object including extra columns, in particular the `padj`
#'   column.
#' 
#' @author Pierre-Luc Germain
#' @references 
#' Hu, J. X. and Zhao, H. and Zhou, H. H. (2010). False Discovery Rate 
#'   Control With Groups. J Am Stat Assoc, 105(491):1215-1227.
#' Ignatiadis, N., Klaus, B., Zaugg, J. et al. (2016) Data-driven hypothesis 
#' weighting increases detection power in genome-scale multiple testing. Nat 
#' Methods 13, 577–580.
#' 
#' @details
#' This function contains different methods to create the bulk-based evidence 
#' bins (defined by the `bin.method` argument), as well as different methods to 
#' use this grouping for multiple testing correction `correction.method`.
#' 
#' Here we call a 'hypothesis' a differential expression test on one gene in one
#' cell type. We define the 'contribution' of the cell type to the bulk 
#' expression of the gene as the proportion of the total pseudobulk reads for 
#' that gene that is contributed by the cell type (across all samples).
#' 
#' The following `bin.method` options are available (if `pb` is missing, only 
#' `sig` is available):
#' \itemize{
#' \item{**sig** : the bulk significance values are used as is, eventually 
#'   taking the direction of the logFC into account if provided in `bulkDEA`.
#'   `nbins` are created using quantiles.}
#' \item{**combined** : first, significance-based bins are created in the same
#'   fashion as in `bin.method="sig"`. Each significance bin is then further 
#'   split into genes for which the cell type contributes much to the bulk, and
#'   genes for which the cell type contributes little.}
#' \item{**PAS** (Proportion-Adjusted Significance): for each hypothesis, the 
#'  bulk significance is adjusted based on the cell type contribution to the 
#'  bulk of that gene using `inv.logit( logit(p) * sqrt(c) )` where `p` and `c` 
#'  are respectively the bulk p-value and the proportion of bulk reads 
#'  contributed by the cell type). We then split this covariate into quantile 
#'  bins as is done for the "sig" method. *This is the recommended method.*}
#' \item{**asNA** : for hypotheses for which the cell type contributes little 
#'  to the bulk profile, the covariate (i.e. bulk p-value) is set to NA, 
#'  resulting it in making up its own bin.}
#' }
#' 
#' In all cases, if `useSign=TRUE` (default) and `bulkDEA` contains logFC 
#' information, then whenever the direction of the change is different between 
#' bulk and pseudobulk datasets for a gene in a given cell type, we increase the
#' bulk p-value to 0.7 (if it was below) for that cell type.
#' 
#' Once the bins are created, the following `correction.method` options are 
#' available:
#' 
#' \itemize{
#' \item{**binwise** : the Benjamini-Hochberg (BH) procedure is applied 
#'  separately for each bin. Doing this can lead to an increase in false 
#'  positives if the number of bins is large, and to correct for this the 
#'  resulting adjusted p-values are multiplied by `pmin(1,nbins/rank(p)`. This
#'  results is proper FDR control even across a large number of bins, but the
#'  method is more conservative than others.}
#' \item{**IHW** : The Independent Hypothesis Weighing (IHW) method of 
#'  Ignatiadis et al. (2016) is applied. See \code{\link[IHW]{ihw}}.}
#' \item{**gBH.LSL** and **gBH.TST**: the Grouped BH method of Hu, Zhao and Zhou
#'  (2010) is applied. The method has two options to compute the groups' rate 
#'  of true null hypotheses, LSL and TST, which make the corresponding 
#'  `correction.method` options (see \code{\link{gBH}} for more detail).
#'  *We recommend using `gBH.LSL`*.}
#' }
#' 
#' Each method exists in two flavors: a local one, which is applied for each 
#' cell type separately, and a global one, which is applied once across all cell
#' types (see the `local` argument). We recommend using the local one.
#' 
#' @export
#' @importFrom IHW ihw adj_pvalues
#' @importFrom SummarizedExperiment assays
#' @importFrom stats quantile
#' @importFrom dplyr bind_rows
bbhw <- function(pbDEA, bulkDEA, pb=NULL, local=TRUE, useSign=TRUE, nbins=NULL,
                 bin.method=c("PAS","combined","asNA","sig"), NAsep=TRUE, 
                 correction.method=c("gBH.LSL","IHW","binwise","gBH.TST"),
                 alpha=0.1, nfolds=NULL, BPPARAM=SerialParam(progress=verbose),
                 verbose=TRUE, ...){
  correction.method <- match.arg(correction.method)
  bin.method <- match.arg(bin.method)
  if(bin.method!="sig"){
    if(is.null(pb)){
      stop("The selected bin.method require the `pb` argument to be given. If ",
           "you don't have this data, use `bin.method='sig'` instead.")
    }
  }else{
    pb <- NULL
  }
  
  stopifnot(is.null(nbins) || nbins>1L)
  stopifnot(alpha>0 && alpha<1)
  
  stopifnot(is.data.frame(pbDEA))
  pbDEA <- .homogenizeDEA4bbhw(pbDEA)
  stopifnot(all(c("p_val","gene","cluster_id") %in% colnames(pbDEA)))
  if(is.data.frame(bulkDEA)){
    bulkDEA <- .homogenizeDEA4bbhw(bulkDEA)
    stopifnot("p_val" %in% colnames(bulkDEA))
    if(useSign && "logFC" %in% colnames(bulkDEA) && 
       "logFC" %in% colnames(pbDEA)){
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
  prop.common.genes <- length(ig)/max(length(bulkDEA),
                                      length(unique(pbDEA$gene)))
  if(verbose) .checkCommonGenes(prop.common.genes)

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
  
  pbDEA$bulk <- bulkDEA[as.character(pbDEA$gene)]
  if(useSign && any(pbDEA$bulk<0)){
    w <- which(sign(pbDEA$logFC) != sign(pbDEA$bulk))
    pbDEA$bulk[w] <- pmax(pbDEA$bulk[w],0.7)
  }
  pbDEA$bulk <- abs(pbDEA$bulk)
  if(!NAsep) pbDEA$bulk[which(is.na(pbDEA$bulk) | pbDEA$bulk>1)] <- 1

  pbDEA$tmpROWindex <- seq_len(nrow(pbDEA))
  
  if(local){
    pbDEA <- dplyr::bind_rows(bplapply(split(pbDEA, pbDEA$cluster_id),
                                       BPPARAM=BPPARAM, FUN=function(pbDEA){
      pbDEA <- .bbhwGetBins(pbDEA, bin.method=bin.method, nbins=nbins)
      .bbhw(pbDEA, correction.method=correction.method, nfolds=nfolds, 
            alpha=alpha, ...)
    }))
  }else{
    pbDEA <- .bbhwGetBins(pbDEA, bin.method=bin.method, nbins=nbins)
    pbDEA <- .bbhw(pbDEA, correction.method=correction.method, nfolds=nfolds, 
                   alpha=alpha, ...)
  }
  pbDEA <- pbDEA[order(pbDEA$tmpROWindex),]
  pbDEA$tmpROWindex <- NULL
  pbDEA
}

.bbhwGetBins <- function(pbDEA, bin.method, nbins, binBy=NULL){
  if(!is.null(binBy)){
    pbDEA <- dplyr::bind_rows(lapply(split(pbDEA, binBy), \(x){
      .bbhwGetBins(x, bin.method=bin.method, nbins=nbins)
    }))
    pbDEA$hbin <- factor(paste(pbDEA$binBy, pbDEA$hbin))
    return(pbDEA)
  }
  if(length(w <- which(is.na(pbDEA$bulk)))>=50){
    pbDEA$bulk[w] <- 2
  }else{
    pbDEA$bulk[w] <- 0.5
  }
  if(bin.method=="combined"){
    pbDEA <- .bbhwGetCombinedBins(pbDEA, nbins)
  }else if(bin.method=="asNA"){
    w <- which(pbDEA$readProportion < mean(pbDEA$readProportion, na.rm=TRUE))
    b2 <- pbDEA$bulk
    b2[w] <- 2L
    if(is.null(nbins)) nbins <- max(1, min(15, floor(nrow(pbDEA)/1000))-1L)
    pbDEA$hbin <- cut(b2, .getQBreaks(b2, nbins), include.lowest=TRUE)
  }else if(bin.method=="sig"){
    if(is.null(nbins)) nbins <- max(1, min(12, floor(nrow(pbDEA)/1000))-1L)
    breaks <- .getQBreaks(pbDEA$bulk, nbins)
    pbDEA$hbin <- cut(pbDEA$bulk, breaks, include.lowest=TRUE)
  }else{
    pbDEA <- .getAdjustedBins(pbDEA, nbins)
  }
  pbDEA$hbin[which(is.na(pbDEA$hbin) | pbDEA$bulk>1)] <- 
    rev(levels(pbDEA$hbin))[1]
  pbDEA$hbin <- droplevels(pbDEA$hbin)
  pbDEA
}

.bbhw <- function(pbDEA, correction.method, nfolds, alpha, verbose=FALSE, ...){
  stopifnot(!is.null(pbDEA$hbin))
  
  if(length(unique(pbDEA$hbin))==1){
    pbDEA$padj <- p.adjust(pbDEA$p_val, method="fdr")
    return(pbDEA)
  }
  
  if(correction.method=="binwise"){
    pbDEA <- dplyr::bind_rows(lapply(split(pbDEA, pbDEA$hbin), FUN=\(x){
      x$padj <- p.adjust(x$p_val, method="fdr")
      x
    }))
    pr <- rank(pbDEA$FDR, ties.method = "min")
    pbDEA$padj <- pbDEA$padj*pmax(1,pr/length(unique(pbDEA$hbin)))
  }else if(correction.method=="gBH.LSL"){
    pbDEA$padj <- gBH(pbDEA$p_val, pbDEA$hbin, pi0="LSL")
  }else if(correction.method=="gBH.TST"){
    pbDEA$padj <- gBH(pbDEA$p_val, pbDEA$hbin, pi0="TST", alpha=alpha)
  }else if(correction.method=="IHW"){
    a <- IHW::ihw(pbDEA$p_val, pbDEA$hbin, alpha=alpha,
                  nfolds=.getNfolds(pbDEA$p_val, nfolds), ...)
    if(verbose) .checkIhwRes(a)
    pbDEA$padj <- IHW::adj_pvalues(a)
  }else{
    stop("Unknown method")
  }
  pbDEA
}

.getQBreaks <- function(x, nbins){
  unique(c(0,quantile(x, seq_len(nbins-1)/(nbins)), 1,
             ifelse(sum(x>1)>=50,3,1)))
}

.getNfolds <- function(x, nfolds){
  if(!is.null(nfolds)) return(nfolds)
  if(length(x)==1) return(ifelse(x>1e5, 5L, ifelse(x>1e4, 4L, 3L)))
  x <- 
  q <- quantile(table(x))
  if(q[[2]]>1000) return(5L)
  if(q[[2]]>400) return(4L)
  return(3L)
}

.checkCommonGenes <- function(prop.common.genes){
  if(prop.common.genes<0.5){
    warning("`pbDEA` and `bulkDEA` have a lot of genes (", 
            round(100*(1-prop.common.genes)), "%) that are not in common.")
  }else if(prop.common.genes<1){
    message(round(100*(prop.common.genes)), "% genes in common")
  }
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

.getAdjustedBins <- function(pbDEA, nbins){
  p <- pmin(pbDEA$bulk,1)
  lip <- log(p/(1 - p))*sqrt(pbDEA$readProportion)
  p <- exp(lip)/(1 + exp(lip))
  p[which(pbDEA$bulk>1 | is.na(pbDEA$bulk))] <- 2
  pbDEA$combCov <- ifelse(is.na(p) & !is.na(pbDEA$bulk), 1, p)
  pbDEA$combCov[is.na(pbDEA$combCov)] <- 1
  w <- which(pbDEA$bulk>1)
  pbDEA$combCov[w] <- ifelse(length(w)>50,2,0.5)
  
  if(is.null(nbins))
    nbins <- max(1, min(12, floor(nrow(pbDEA)/1000))-1L)
    breaks <- c(0,quantile(pbDEA$combCov[pbDEA$combCov<=1],
                         seq_len(nbins)[-ceiling((nbins+1L)/2)]/(nbins+1L),
                         na.rm=TRUE), 1, ifelse(sum(pbDEA$combCov>1)>50,3,1))
  pbDEA$hbin <- cut(pbDEA$combCov, unique(breaks), include.lowest=TRUE)
  pbDEA$hbin[which(is.na(pbDEA$hbin) | pbDEA$bulk>1)] <- 
    rev(levels(pbDEA$hbin))[1]
  pbDEA$hbin <- droplevels(pbDEA$hbin)
  pbDEA
}

.bbhwGetCombinedBins <- function(pbDEA, nbins=NULL){
  if(is.null(nbins)) nbins <- max(1, min(6, floor(nrow(pbDEA)/2000)))
  hasNAbin <- length(w <- which(pbDEA$bulk>1))>=50
  if(!hasNAbin) pbDEA$bulk[w] <- 0.5
  breaks <- .getQBreaks(pbDEA$bulk, nbins)
  pbDEA$sigBin <- cut(pbDEA$bulk, unique(breaks), include.lowest=TRUE)
  pbDEA <- dplyr::bind_rows(lapply(split(pbDEA, pbDEA$sigBin), FUN=\(x){
    if(nrow(x)<=2000 || 
       (hasNAbin & x$sigBin[[1]]==rev(levels(pbDEA$sigBin))[1])){
      x$propBin <- FALSE
    }else{
      q <- quantile(x$readProportion, ifelse(nrow(x)>4000,0.75,0.5))
      x$propBin <- x$readProportion >= as.numeric(q)[1]
    }
    x
  }))
  pbDEA$hbin <- as.factor(paste(pbDEA$sigBin, pbDEA$propBin))
  pbDEA$hbin[is.na(pbDEA$hbin)] <- rev(levels(pbDEA$hbin))[1]
  pbDEA
}

.homogenizeDEA4bbhw <- function(x){
  if(is.null(x$p_val) && !is.null(x$PValue)) x$p_val <- x$PValue
  if(is.null(x$cluster_id) && !is.null(x$celltype)) x$cluster_id <- x$celltype
  x
}


#' gBH - Grouped Benjamini-Hochberg procedure
#'
#' This computes adjusted p-values using a user-defined grouping of the 
#' hypotheses, following the method from Hu, Zhao and Zhou (2010).
#'
#' @param p A vector of p-values
#' @param bins A factor of same length as `p` indicating to which bin the 
#'  p-value belongs
#' @param pi0 The pi0 estimation method, either LSL (default) or TST.
#' @param alpha The desired FDR control (ignored for the LSL method)
#' 
#' @return A vector of same length as `p` with the adjusted p-values.
#' 
#' @author Pierre-Luc Germain
#' @details
#' This is partly inspired from code in the c212 package by Raymond Carragher,
#' which followed the implementation described in Hu, Zhao and Zhou (2010). The 
#' implementation was adapted to use vector operations to be (a lot) 
#' faster, and to produce adjusted p-values (rather than a rejection rule).
#' 
#' The method has two variants which differ in the way the rate of true null
#' hypotheses (pi0) in each group/bin is estimated. The Two-Stage (TST) method 
#' uses Benjamini & Hochberg's FDR method to estimate the proportion of 
#' rejections in each group, and bases the pi0 on this. The LSL method uses
#' the Least-Slope estimator proposed by Benjamini and Hochberg (2000). We 
#' recommend using the LSL method (default), which was more robust in our hands
#' and has the virtue of not being dependent on an input alpha.
#'
#' @references 
#' Hu, J. X. and Zhao, H. and Zhou, H. H. (2010). False Discovery Rate 
#'   Control With Groups. J Am Stat Assoc, 105(491):1215-1227.
#' Benjamini Y, Hochberg Y. (2000). On the Adaptive Control of the False 
#'   Discovery Rate in Multiple Testing With Independent Statistics. Journal of 
#'   Educational and Behavioral Statistics, 25(1):60–83.
#'
#' @return A vector of adjusted p-values
#' @export
#' @importFrom stats p.adjust
#' @examples
#' # generate data with fake p-values and bins that are somewhat informative:
#' d <- data.frame(
#'   p=c(abs(rnorm(500,sd=0.01)), runif(4500)),
#'   truth=rep(c(TRUE,FALSE), c(500,4500)),
#'   bins=c(sample(LETTERS[1:10], 500, c(0.55,rep(0.05,9)), replace=TRUE),
#'          sample(LETTERS[1:10], 4500, replace=TRUE)))
#' # compute grouped adjusted p-values:
#' d$padj_grouped <- gBH(d$p, d$bins)
#' table(truth=d$truth, grouped=d$padj_grouped<0.05)
#' # compare with the normal BH:
#' d$padj_normal <- p.adjust(d$p, method="fdr")
#' table(truth=d$truth, normal=d$padj_normal<0.05)
gBH <- function(p, bins, pi0=c("LSL","TST"), alpha=0.05){
  if(is.character(bins)) bins <- as.factor(bins)
  stopifnot(is.factor(bins))
  bins <- droplevels(bins)
  stopifnot(length(bins) == length(p))
  stopifnot(!any(is.na(bins)) && !any(is.na(p)))
  pi0 <- match.arg(pi0)
  if(pi0!="LSL")
    stopifnot(length(alpha)==1 && alpha < 1 && alpha > 0)
  
  # estimate the proportion of true null hypotheses in each bin
  if(pi0=="TST"){
    pi0_est <- vapply(split(p, bins), FUN.VALUE=numeric(1), FUN=\(x){
      1-sum(p.adjust(x, method="fdr")<=(alpha/(1 + alpha)))/length(x)
    })
  }else{
    pi0_est <- vapply(split(p, bins), FUN.VALUE=numeric(1), FUN=\(x){
      x <- sort(x)
      l <- rev(seq_along(x))/(1L-x)
      j <- 1L
      if(length(x)>1) j <- head(which(l[-1] > l[-length(l)]),1)
      min((floor(l[j]) + 1L)/length(x), 1L)
    })
  }
  if(all(pi0_est==1)) return(rep(1L, length(p)))
  # weigh p-values
  w <- as.numeric(pi0_est/(1-pi0_est))
  pw <- p * w[as.integer(bins)]
  # compute BH on pooled weighted p-values, using the weighted alpha, i.e.
  # alpha/(1-pi_global)
  o <- order(pw, decreasing=TRUE)
  global_pi0 <- sum(as.numeric(pi0_est)*as.integer(table(bins)))/length(p)
  padj <- (1-global_pi0)*length(p)*pw[o]/rev(seq_along(pw))
  # report the smallest alpha for which each hypothesis is rejected
  pmin(1, cummin(padj)[order(o)])
}
