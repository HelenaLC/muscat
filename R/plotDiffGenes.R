#' @rdname plotDiffGenes
#' @title Heatmap of mean-marker expression by cluster-sample
#' 
#' @description ...
#' 
#' @param x a \code{[SingleCellExperiment]{SingleCellExperiment}}.
#' @param y a table containing differential testing results
#'   as returned by, for example, \code{run_edgeR}.
#' @param top_n single numeric specifying the number of genes to include.
#' @param sort_by a character string specifying a \code{y} column to sort by.
#' @param clusters a character string specifying which cluster(s) results
#'   should be included for. If NULL (the default), the \code{top_n} hits
#'   for each cluster will be included.
#' @param fdr single numeric specifying the threshold on adjusted p-values
#'   below which results should be considered significant.
#' @param lfc single numeric specifying the threshold on absolute logFCs
#'   above which results should be included.
#' 
#' @return a \code{\link{HeatmapList-class}} object.
#' 
#' @author Helena L. Crowell \email{helena@crowells.eu}
#' 
#' @import ComplexHeatmap
#' @import SingleCellExperiment
#' @importFrom dplyr %>% group_by summarise_all
#' @importFrom grid gpar
#' @importFrom methods is
#' @importFrom reshape2 melt
#' @importFrom tidyr complete
#' @importFrom viridis viridis
#' @export

plotDiffGenes <- function(x, y, 
    top_n = 20, sort_by = c("FDR", "logFC"), 
    clusters = NULL, fdr = 0.05, lfc = 1) {
    
    # validity checks of input arguments
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot(is.numeric(top_n), length(top_n) == 1, top_n > 1)
    sort_by <- match.arg(sort_by)
    
    cluster_ids <- colData(x)$cluster_id
    contrasts <- levels(y$contrast)
    
    # filter results
    if (is.null(clusters)) clusters <- levels(cluster_ids)
    y <- y %>% filter(FDR < fdr, abs(logFC) > lfc, cluster_id %in% clusters)
    
    # order results
    if (sort_by == "FDR") {
        o <- order(y$FDR)
    } else if (sort_by == "logFC") {
        o <- order(abs(y$logFC), decreasing = TRUE)
    }
    y <- y[o, ]
    
    # split results by contrast & cluster
    y <- split(y, y$contrast)
    y <- lapply(y, function(x) {
        rownames(x) <- sprintf("%s(%s)", x$gene, x$cluster_id)
        return(x)
    })
    y <- lapply(y, function(x) split(x, x$cluster_id))
    
    # subset top_n hits
    y <- lapply(y, lapply, "[", i = seq_len(top_n), j = TRUE)
    y <- lapply(y, lapply, function(x) x[!is.na(x$gene), ])
    
    gs <- lapply(y, lapply, "[[", "gene")
    gs <- unlist(Reduce(c, gs))
    gs <- make.names(gs)
    gs <- unique(gs)
    
    rownames(x) <- gsub("[-]|[(]|[)]", ".", rownames(x))
    es <- logcounts(x[gs, ])
    es <- as.matrix(t(es))
    es <- CATALYST:::scale_exprs(es)
    df <- data.frame(
        check.names = FALSE, es,
        cluster_id = colData(x)$cluster_id,
        sample_id = colData(x)$sample_id)
    
    zeros <- setNames(as.list(numeric(nrow(x))), rownames(x))
    ms <- df %>% 
        group_by_(~cluster_id, ~sample_id) %>% 
        summarise_at(gs, mean) %>% 
        complete(sample_id, fill = zeros) %>% 
        melt(id.var = c("cluster_id", "sample_id"),
            variable.name = "gene", value.name = "expr")
    
    # reformat (rows = genes, columns = samples)
    ms <- split(ms, ms$gene)
    ms <- lapply(ms, acast, cluster_id ~ sample_id, value.var = "expr")
    ms <- do.call(rbind, ms)
    rownames(ms) <- sprintf("%s(%s)",
        rep(gs, each = nlevels(cluster_ids)),
        rep(levels(cluster_ids), length(gs)))
    
    # split by contrast
    ms <- setNames(lapply(contrasts, function(i)
        ms[unlist(sapply(y[[i]], rownames)), ]),
        contrasts)
    
    # retain only samples of compaired groups & reorder
    ei <- metadata(x)$experiment_info
    ei$group <- gsub("[-]|[(]|[)]|[/]", ".", ei$group)
    g1 <- setNames(gsub(".*-", "", contrasts), contrasts)
    g2 <- setNames(gsub("-.*", "", contrasts), contrasts)
    groups <- sapply(contrasts, function(c)
        unlist(sapply(list(g1 = g1[c], g2 = g2[c]), function(x) 
            as.character(ei$sample_id[ei$group %in% x]))))
    ms <- setNames(lapply(contrasts, function(i) ms[[i]][, groups[, i]]), contrasts)
    
    lapply(contrasts, function(i) {
        Heatmap(ms[[i]],
            column_title = i,
            col = viridis(10),
            name = "avg. scaled\nexpression",
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            row_names_gp = gpar(fontsize = 6))
    })
}