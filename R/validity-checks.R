# validity checks for objects & function arguments

# check validity of runDS() output
check_res <- function(sce, res) {
    ei <- metadata(sce)$experiment_info
    cluster_ids <- levels(factor(colData(sce)$cluster_id))
    n_clusters <- length(cluster_ids)
    nms <- c("table", "data", "design", "contrast", "coef")
    stopifnot(
        is(res, "list"), all.equal(names(res), nms),
        # table
        apply(vapply(res$table, names, character(n_clusters)), 
            2, function(u) all(u %in% cluster_ids)),
        # data
        is(res$data, "list"),
        names(res$data) %in% cluster_ids,
        vapply(res$data, function(u) is(u, "DGEList"), logical(1)),
        # design
        is(res$design, "matrix"), 
        colnames(res$design) %in% ei$group_id,
        rownames(res$design) %in% ei$sample_id,
        # contrast & coef
        is.null(res$contrast) | is(res$contrast, "matrix"),
        is.null(res$coef) | is(res$coef, "numeric") | is(res$coef, "list"))
}

# check validity of calcExprFreqs() output
check_frq <- function(sce, frq) {
    cluster_ids <- levels(factor(colData(sce)$cluster_id))
    group_ids <- levels(colData(sce)$group_id)
    sample_ids <- levels(colData(sce)$sample_id)
    nms <- vapply(frq, function(u) sort(names(u)), 
        character(length(group_ids) + length(sample_ids)))
    stopifnot(
        length(frq) == length(cluster_ids),
        names(frq) %in% cluster_ids,
        apply(nms, 1, function(u) length(unique(u)) == 1),
        apply(nms, 2, function(u) u %in% c(group_ids, sample_ids)))
}
