#' resDS
#' Formatting of DS analysis results
#' 
#' \code{resDS} provides a simple wrapper to format cluster-level
#' differential testing results into an easily filterable table, and
#' to optionally append gene expression frequencies by cluster-sample
#' & -group, as well as cluster-sample-wise CPM.
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param y a list of cluster-level differential testing results 
#'   as returned by \code{\link{runDS}}.
#' @param bind character string specifying the output format (see details).
#' @param frq logical or a pre-computed list of expression frequencies 
#'   as returned by \code{\link{clacExprFreqs}}.
#' @param cpm logical specifying whether CPM by cluster-sample 
#'   should be appendeded to the output result table(s).
#' @param digits integer value specifying the 
#'   number of significant digits to maintain.
#' @param ... optional arguments passed to 
#'   \code{\link{clacExprFreqs}} if \code{frq = TRUE}.
#' 
#' @details When \code{bind = "col"}, the list of DS testing results at 
#'   \code{y$table} will be merge vertically (by column) into a single table
#'   in tidy format with column \code{contrast/coef} specifying the comparison.
#'   
#'   Otherwise, when \code{bind = "row"}, an identifier of the respective
#'   contrast or coefficient will be appended to the column names,
#'   and all tables will be merge horizontally (by row).
#'   
#'   Expression frequencies pre-computed with \code{\link{calcExprFreqs}} 
#'   may be provided with \code{frq}. Alternatively, when \code{frq = TRUE}, 
#'   expression frequencies can be computed directly, and additional arguments 
#'   may be passed to \code{\link{calcExprFreqs}} (see examples below).
#' 
#' @return returns a data.frame.
#' 
#' @examples
#' # simulate 5 clusters, 20% of DE genes
#' data(kang)
#' sim <- simData(kang, n_genes = 100, n_cells = 200, 
#'   p_dd = c(0.8, 0, 0.2, 0, 0, 0), fc = 4)
#' 
#' # compute pseudo-bulk counts
#' pb <- aggregateData(sim, data = "counts", fun = "sum")
#' 
#' # specify design & contrast matrix
#' ei <- S4Vectors::metadata(sim)$experiment_info
#' design <- stats::model.matrix(~ 0 + ei$group_id)
#' dimnames(design) <- list(ei$sample_id, levels(ei$group_id))
#' contrast <- limma::makeContrasts("B-A", levels = design)
#' 
#' # test for cluster-specific DE 
#' res <- runDS(sim, pb, design, contrast, method = "edgeR")
#' 
#' head(resDS(sim, res, bind = "row"))
#' head(resDS(sim, res, bind = "col", digits = Inf))
#' 
#' head(resDS(sim, res, cpm = TRUE)) # append CPM by sample
#' head(resDS(sim, res, frq = TRUE)) # append expr. freqs. by sample
#' 
#' # pre-computed expr. freqs. & append
#' frq <- calcExprFreqs(sim, th = 10)
#' head(resDS(sim, res, frq = frq))
#' 
#' @author Helena L. Crowell \email{helena.crowell@uzh.ch} and Mark D. Robinson.
#' 
#' @importFrom dplyr %>% bind_rows inner_join full_join mutate mutate_if select
#' @importFrom edgeR cpm
#' @importFrom methods is
#' @importFrom purrr reduce
#' @importFrom reshape2 melt
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#'  
#' @export  

resDS <- function(x, y, bind = c("col", "row"),
    frq = FALSE, cpm = FALSE, digits = 3, ...) {
    
    stopifnot(is(x, "SingleCellExperiment"))
    cluster_ids <- levels(colData(x)$cluster_id)
    ei <- metadata(x)$experiment_info
    
    # check_res(x, y)
    # if (!is.logical(frq)) 
    #     check_frq(x, frq)
    bind <- match.arg(bind)
    stopifnot(is.infinite(digits) || is.numeric(digits) &
        digits > 0 & as.integer(digits) == digits)

    res <- switch(bind,
        row = {
            ct <- ifelse(!is.null(y$contrast), "contrast", "coef")
            cs <- names(y$table)
            res <- lapply(cs, function(c) {
                df <- bind_rows(y$table[[c]])
                df <- select(df, -ct)
                i <- !colnames(df) %in% c("gene", "cluster_id")
                colnames(df)[i] <- paste(colnames(df)[i], c, sep = "__")
                return(df)
            })
            res %>% reduce(full_join, 
                by = c("gene", "cluster_id"))
        },
        col = {
            bind_rows(lapply(y$table, bind_rows))
        })

    reorder_summary <- function(u, ei, append = "") {
        m1 <- match(ei$sample_id, colnames(u))
        m2 <- match(levels(ei$group_id), colnames(u))
        if (all(is.na(m2))) m2 <- 0
        colnames(u)[m1] <- paste0(ei$sample_id, "__", ei$group_id, append)
        colnames(u)[m2] <- paste0(colnames(u)[m2], append)
        k <- seq_len(ncol(u))[-c(m1, m2)]
        u[, c(k, m1[order(ei$group)], m2)]
    }
    
    # append expression frequencies
    if (is.logical(frq))
        if (frq) frq <- calcExprFreqs(x, ...) else frq <- NULL
    if (!is.null(frq)) {
        frq <- data.frame(gene = rep(rownames(x), length(frq)),
            cluster_id = rep(names(frq), each = nrow(x)), bind_rows(frq), 
            row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
        frq <- reorder_summary(frq, ei, append = ".frq")
        res <- inner_join(frq, res, by = c("gene", "cluster_id"))
    }

    # append CPMs if available
    if (cpm) {
        cpm <- lapply(cluster_ids, function(k) {
            if (is.null(y$data[[k]])) return(NULL)
            cpm <- cpm(y$data[[k]])
            data.frame(gene = rownames(cpm), cluster_id = k,
                cpm, row.names = NULL, stringsAsFactors = FALSE)
        })
        cpm <- bind_rows(cpm)
        cpm <- reorder_summary(cpm, ei, append = ".cpm")
        res <- inner_join(frq, res, by = c("gene", "cluster_id"))
    }
    res %>% mutate_if(is.numeric, signif, digits)
}