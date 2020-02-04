# validity checks for objects & function arguments
# ==============================================================================

# check input SCE
#' @importFrom methods is
#' @importFrom SummarizedExperiment colData
.check_sce <- function(x, req_group = TRUE) {
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot(c("cluster_id", "sample_id") %in% colnames(colData(x)))
    if (req_group)
        stopifnot("group_id" %in% colnames(colData(x)))
}

# check of 'assay' argument
#' @importFrom SummarizedExperiment assayNames
.check_arg_assay <- function(x, y) {
    stopifnot(is.character(y), length(y) == 1, y %in% assayNames(x))
    if (sum(assayNames(x) == y) > 1)
        stop("Argument 'assay' was matched to multiple times.\n ", 
            " Please assure that the input SCE has unique 'assayNames'.")
}

# check pseudo-bulks for DS analysis
# (must have be aggregated by cluster-sample)
#   x = SCE used for aggregation
#   y = SCE containing pseudo-bulks as returned by 
#`      aggregateData(x, by = c("cluster_id", "sample_id"))
#' @importFrom methods is
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment assayNames
.check_pbs <- function(pbs, sce = NULL, check_by = TRUE) {
    stopifnot(is(pbs, "SingleCellExperiment"),
        !is.null(ei <- metadata(pbs)$experiment_info),
        !is.null(agg_pars <- metadata(pbs)$agg_pars),
        !is.null(n_cells <- metadata(pbs)$n_cells),
        identical(assayNames(pbs), rownames(n_cells)),
        identical(colnames(pbs), colnames(n_cells)))
    if (!is.null(sce)) {
        stopifnot(identical(ei, metadata(sce)$experiment_info),
            identical(assayNames(pbs), levels(sce[[agg_pars$by[1]]])),
            identical(rownames(pbs), rownames(sce)))
        if (length(agg_pars$by == 2))
            stopifnot(identical(colnames(pbs), levels(sce[[agg_pars$by[2]]])))
    }
    if (check_by)
        stopifnot(!is.null(pbs[["group_id"]]),
            identical(agg_pars$by, c("cluster_id", "sample_id")))
}

# check validity of runDS() output
#' @importFrom methods is
#' @importFrom S4Vectors metadata
.check_res <- function(x, y) {
    ei <- metadata(x)$experiment_info
    nk <- length(kids <- levels(x$cluster_id))
    nms <- c("table", "data", "method", "design", "contrast", "coef")
    stopifnot(is(y, "list"), all.equal(names(y), nms))
    # table
    stopifnot(is(y$table, "list"),
        vapply(y$table, is, class = "list", logical(1)),
        identical(names(y$table), colnames(y$contrast))
        | identical(names(y$table), names(y$coef)),
        apply(vapply(y$table, names, character(nk)), 2, identical, kids))
    # data
    stopifnot(is(y$data, "list"), names(y$data) %in% kids,
        vapply(y$data, is, class = "DGEList", logical(1)))
    # design
    stopifnot(is(y$design, "matrix"),
        colnames(y$design) %in% ei$group_id,
        rownames(y$design) %in% ei$sample_id)
    # contrast & coef
    stopifnot(is.null(y$contrast) | is(y$contrast, "matrix"))
    stopifnot(is.null(y$coef) | is(y$coef, "numeric") | is(y$coef, "list"))
}

# check validity of calcExprFreqs() output
#' @importFrom methods is
#' @importFrom SummarizedExperiment assays colData
.check_frq <- function(x, y) {
    stopifnot(
        is(x, "SingleCellExperiment"), 
        is(y, "SingleCellExperiment"))
    kids <- levels(x$cluster_id)
    
    ids <- levels(x$sample_id)
    if ("group_id" %in% colnames(colData(x)))
        ids <- c(ids, levels(x$group_id))
    stopifnot(identical(ids, colnames(y)))
    
    vals <- unlist(assays(y))
    stopifnot(all(vals <= 1), all(vals >= 0))
}

.check_args_simData <- function(u) {
    stopifnot(
        is.numeric(u$nc), length(u$nc) == 1, u$nc > 0, as.integer(u$nc) == u$nc,
        is.numeric(u$nk), length(u$nk) == 1, u$nk > 0, as.integer(u$nk) == u$nk,
        is.numeric(u$ns), length(u$ns) %in% c(1, 2), u$ns > 0, as.integer(u$ns) == u$ns,
        is.numeric(u$p_dd), length(u$p_dd) == 6, sum(u$p_dd) == 1, u$p_dd >= 0, u$p_dd <= 1,
        is.logical(u$paired), length(u$paired) == 1,
        is.numeric(u$p_ep), length(u$p_ep) == 1, u$p_ep > 0, u$p_ep < 1,
        is.numeric(u$p_dp), length(u$p_dp) == 1, u$p_dp > 0, u$p_dp < 1,
        is.numeric(u$p_dm), length(u$p_dm) == 1, u$p_dm > 0, u$p_dm < 1,
        is.numeric(u$p_type), length(u$p_type) == 1, u$p_type >= 0, u$p_type <= 1,
        is.numeric(u$lfc), is.numeric(u$lfc), length(u$lfc) == 1, u$lfc >= 1,
        is.numeric(u$ng), length(u$ng) == 1, u$ng > 0, as.integer(u$ng) == u$ng,
        is.logical(u$force), length(u$force) == 1)
    if (!is.null(u$rel_lfc))
        stopifnot(is.numeric(u$rel_lfc), 
            length(u$rel_lfc) == u$nk, u$rel_lfc >= 0)
}

#' @importFrom SummarizedExperiment colData
.check_args_aggData <- function(u) {
    stopifnot(is.character(u$by), length(u$by) <= 2, 
        u$by %in% colnames(colData(u$x)))
    stopifnot(is.logical(u$scale), length(u$scale) == 1)
    if (u$scale & (!u$assay %in% c("cpm", "CPM") | u$fun != "sum"))
        stop("Option 'scale = TRUE' only valid for", 
            " 'assay = \"cpm/CPM\"' and 'fun = \"sum\"'.")
}

.check_args_pbDS <- function(u) {
    if (!is.null(u$design))
        stopifnot(is.matrix(u$design),
            !is.null(rownames(u$design)),
            !is.null(colnames(u$design)))
    stopifnot(
        is.null(u$contrast) | is.matrix(u$contrast),
        is.null(u$coef) | is.numeric(u$coef),
        is.numeric(u$min_cells), length(u$min_cells) == 1,
        is.logical(u$verbose), length(u$verbose) == 1,
        is.logical(u$treat), length(u$treat) == 1)
}

.check_args_pbHeatmap <- function(u) {
    if (!is.null(u$k))
        stopifnot(is.character(u$k), u$k %in% levels(u$x$cluster_id))
    if (!is.null(u$g))
        stopifnot(is.character(u$g), u$g %in% rownames(u$x))
    if (!is.null(u$c))
        stopifnot(is.character(u$c), u$c %in% names(u$y$table))
    stopifnot(
        is.numeric(u$top_n), length(u$top_n) == 1, u$top_n > 1,
        is.numeric(u$fdr), length(u$fdr) == 1, u$fdr > 0,
        is.numeric(u$lfc), length(u$lfc) == 1)
    stopifnot(
        is.character(u$sort_by), 
        length(u$sort_by) == 1,
        u$sort_by %in% names(u$y$table[[1]][[1]]),
        is.numeric(u$y$table[[1]][[1]][[u$sort_by]]))
    stopifnot(is.function(u$fun))
    stopifnot(is.logical(u$normalize), length(u$normalize) == 1)
    stopifnot(
        is.logical(u$row_anno), length(u$row_anno) == 1,
        is.logical(u$col_anno), length(u$col_anno) == 1)
}
