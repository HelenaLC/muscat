#' resDS
#' Formatting of DS analysis results
#' 
#' \code{resDS} provides a simple wrapper to format cluster-level
#' differential testing results into an easily filterable table, and
#' to optionally append gene expression frequencies by cluster-sample
#' & -group, as well as cluster-sample-wise CPM.
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param y a list of DS testing results as returned 
#'   by \code{\link{pbDS}} or \code{\link{mmDS}}.
#' @param bind character string specifying the output format (see details).
#' @param frq logical or a pre-computed list of expression frequencies 
#'   as returned by \code{\link{calcExprFreqs}}.
#' @param cpm logical specifying whether CPM by cluster-sample 
#'   should be appendeded to the output result table(s).
#' @param digits integer value specifying the 
#'   number of significant digits to maintain.
#' @param sep character string to use as separator 
#'   when constructing new column names.
#' @param ... optional arguments passed to 
#'   \code{\link{calcExprFreqs}} if \code{frq = TRUE}.
#' 
#' @details When \code{bind = "col"}, the list of DS testing results at 
#' \code{y$table} will be merge vertically (by column) into a single table
#' in tidy format with column \code{contrast/coef} specifying the comparison.
#'   
#' Otherwise, when \code{bind = "row"}, an identifier of the respective
#' contrast or coefficient will be appended to the column names,
#' and all tables will be merge horizontally (by row).
#'   
#' Expression frequencies pre-computed with \code{\link{calcExprFreqs}} 
#' may be provided with \code{frq}. Alternatively, when \code{frq = TRUE}, 
#' expression frequencies can be computed directly, and additional arguments 
#' may be passed to \code{\link{calcExprFreqs}} (see examples below).
#' 
#' @return returns a `data.frame`.
#' 
#' @examples
#' data(sce)
#' 
#' # compute pseudobulks (sum of counts)
#' pb <- aggregateData(sce, assay = "counts", fun = "sum")
#' 
#' # run DS analysis (edgeR on pseudobulks)
#' res <- pbDS(pb, method = "edgeR")
#' 
#' head(resDS(sce, res, bind = "row")) # tidy format
#' head(resDS(sce, res, bind = "col", digits = Inf))
#' 
#' # append CPMs & expression frequencies
#' head(resDS(sce, res, cpm = TRUE))
#' head(resDS(sce, res, frq = TRUE))
#' 
#' # pre-computed expression frequencies & append
#' frq <- calcExprFreqs(sce, assay = "counts", th = 0)
#' head(resDS(sce, res, frq = frq))
#' 
#' @author Helena L Crowell & Mark D Robinson
#' 
#' @importFrom dplyr %>% bind_rows inner_join full_join mutate mutate_if select
#' @importFrom edgeR cpm
#' @importFrom methods is
#' @importFrom purrr reduce
#' @importFrom reshape2 melt
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#' @export  

resDS <- function(x, y, bind = c("col", "row"),
    frq = FALSE, cpm = FALSE, digits = 3, sep = "__", ...) {
    
    # check validity of input arguments
    .check_sce(x, req_group = TRUE)
    #.check_res(x, y)
    bind <- match.arg(bind)
    if (!is.logical(frq)) 
        .check_frq(x, frq)
    stopifnot(is.infinite(digits) || is.numeric(digits) &
        digits > 0 & as.integer(digits) == digits)

    ei <- metadata(x)$experiment_info
    kids <- levels(x$cluster_id)
    
    res <- switch(bind,
        row = {
            ct <- ifelse(!is.null(y$contrast), "contrast", "coef")
            cs <- names(y$table)
            res <- lapply(cs, function(c) {
                df <- bind_rows(y$table[[c]])
                df <- select(df, -ct)
                i <- !colnames(df) %in% c("gene", "cluster_id")
                colnames(df)[i] <- paste(colnames(df)[i], c, sep = sep)
                return(df)
            })
            reduce(res, full_join, by = c("gene", "cluster_id"))
        },
        col = {
            bind_rows(lapply(y$table, bind_rows))
        })

    .tidy <- function(u, ei, append = "") {
        m1 <- match(ei$sample_id, colnames(u))
        m2 <- match(levels(ei$group_id), colnames(u))
        if (all(is.na(m2))) m2 <- 0
        colnames(u)[m1] <- paste0(ei$sample_id, append)
        colnames(u)[m2] <- paste0(colnames(u)[m2], append)
        k <- seq_len(ncol(u))[-c(m1, m2)]
        u[, c(k, m1[order(ei$group)], m2)]
    }
    
    # append expression frequencies
    if (is.logical(frq))
        if (frq) frq <- calcExprFreqs(x, ...) else frq <- NULL
    if (!is.null(frq)) {
        frq <- data.frame(
            gene = rep(rownames(x), length(assays(frq))),
            cluster_id = rep(assayNames(frq), each = nrow(x)), 
            do.call("rbind", as.list(assays(frq))),
            row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
        frq <- .tidy(frq, ei, append = ".frq")
        res <- inner_join(frq, res, by = c("gene", "cluster_id"))
    }

    # append CPMs
    if (cpm) {
        cpm <- lapply(kids, function(k) {
            if (is.null(y$data[[k]])) return(NULL)
            cpm <- cpm(y$data[[k]])
            data.frame(gene = rownames(cpm), cluster_id = k, cpm, 
                row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
        })
        cpm <- bind_rows(cpm)
        cpm <- .tidy(cpm, ei, append = ".cpm")
        res <- inner_join(cpm, res, by = c("gene", "cluster_id"))
    }
    mutate_if(res, is.numeric, signif, digits)
}
