# validity checks for objects & function arguments
# ==============================================================================

# check input SCE
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

# check validity of runDS() output
.check_res <- function(x, y) {
    ei <- metadata(x)$experiment_info
    kids <- levels(x$cluster_id)
    nk <- length(kids)

    stopifnot(is(y, "list"), all.equal(names(y),
        c("table", "data", "design", "contrast", "coef")))
    # table
    stopifnot(is(y$table, "list"))
    stopifnot(vapply(y$table, function(u) is(u, "list"), logical(1)))
    stopifnot(identical(names(y$table), colnames(y$contrast))
        | identical(names(y$table), names(y$coef)))
    stopifnot(apply(vapply(y$table, names, character(nk)), 2, identical, kids))
    # data
    stopifnot(is(y$data, "list"))
    stopifnot(names(y$data) %in% kids)
    stopifnot(vapply(y$data, function(u) is(u, "DGEList"), logical(1)))
    # design
    stopifnot(is(y$design, "matrix"))
    stopifnot(colnames(y$design) %in% ei$group_id)
    stopifnot(rownames(y$design) %in% ei$sample_id)
    # contrast & coef
    stopifnot(is.null(y$contrast) | is(y$contrast, "matrix"))
    stopifnot(is.null(y$coef) | is(y$coef, "numeric") | is(y$coef, "list"))
}

# check validity of calcExprFreqs() output
.check_frq <- function(x, y) {
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot(is(y, "SummarizedExperiment"))
    stopifnot(identical(levels(x$cluster_id), assayNames(y)))
    
    ids <- levels(x$sample_id)
    if ("group_id" %in% colnames(colData(x)))
        ids <- c(ids, levels(x$group_id))
    stopifnot(identical(ids, colnames(y)))
    
    vals <- unlist(assays(y))
    stopifnot(all(vals <= 1), all(vals >= 0))
}
