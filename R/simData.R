#' simData
#' 
#' Simulation of complex scRNA-seq data 
#' 
#' \code{simData} simulates multiple clusters and samples 
#' across 2 experimental conditions from a real scRNA-seq data set.
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param n_genes # of genes to simulate. 
#' @param n_cells # of cells to simulate. 
#'   Either a single numeric or a range to sample from.
#' @param ns nb. of genes common to 1, 2, ..., all clusters.
#' @param p_dd numeric vector of length 6.
#'   Specifies the probability of a gene being
#'   EE, EP, DE, DP, DM, or DB, respectively.
#' @param fc numeric value to use as mean logFC
#'   for DE, DP, DM, and DB type of genes.
#' 
#' @return a \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   containing multiple clusters & samples across 2 groups.
#' 
#' @examples
#' data(kang)
#' simData(kang,
#'     n_genes = 10, n_cells = 10,
#'     p_dd = c(1,0,0,0,0,0))
#' 
#' @importFrom data.table data.table
#' @importFrom dplyr mutate_all mutate_at
#' @importFrom edgeR DGEList estimateDisp glmFit
#' @importFrom purrr modify_at
#' @importFrom stats model.matrix rgamma setNames
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors split
#' @importFrom tibble column_to_rownames
#' @importFrom zeallot %<-%
#' 
#' @export

simData <- function(x, n_genes = 500, n_cells = 300, probs = NULL, p_dd = diag(6)[1, ], fc = 2) {
    
    # throughout this code...
    # k: cluster ID
    # s: sample ID
    # c: gene category
    
    # check validity of input arguments
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot(is.numeric(n_genes), length(n_genes) == 1)
    stopifnot(is.numeric(n_cells), length(n_cells) == 1 | length(n_cells) == 2)
    stopifnot(is.numeric(p_dd), length(p_dd) == 6, sum(p_dd) == 1)
    stopifnot(is.numeric(fc), is.numeric(fc), fc > 1)
    
    kids <- levels(colData(x)$cluster_id)
    sids <- levels(colData(x)$sample_id)
    gids <- c("A", "B")
    names(kids) <- kids
    names(sids) <- sids
    names(gids) <- gids
    nk <- length(kids)
    
    # initialize count matrix
    gs <- paste0("gene", seq_len(n_genes))
    cs <- paste0("cell", seq_len(n_cells))
    y <- matrix(0, n_genes, n_cells, dimnames = list(gs, cs))
    
    # sample cell metadata
    cd <- .sample_cell_md(
        n = n_cells, probs = NULL,
        ids = list(kids, sids, gids)) %>% set_rownames(cs)
    cs_idx <- .split_cells(cd, by = colnames(cd))
    n_cs <- modify_depth(cs_idx, -1, length)
    
    # split input cells by cluster-sample
    cs_by_ks <- .split_cells(x)
    
    # sample nb. of genes to simulate per category & gene indices
    n_dd <- replicate(nk, 
        table(sample(factor(cats, levels = cats), n_genes, TRUE, p_dd))) %>% 
        set_colnames(kids)
    gs_idx <- .sample_gene_inds(gs, n_dd)
    
    # for ea. cluster, sample unique set of genes to simulate from
    gs_by_k <- replicate(nk, 
        setNames(sample(rownames(x), n_genes, TRUE), gs),
        simplify = FALSE) %>% set_names(kids)
    gs_by_kc <- lapply(kids, function(k) 
        lapply(cats, function(c)
            gs_by_k[[k]][gs_idx[[c, k]]]) %>% 
            set_names(cats))
    
    # sample logFCs
    lfc <- vapply(kids, function(k) 
        lapply(cats, function(c) { 
            n <- n_dd[c, k]
            if (c == "ee") return(rep(NA, n))
            signs <- sample(c(-1, 1), n, TRUE)
            lfc <- rgamma(n, 4, 4/fc) * signs
            names(lfc) <- gs_by_kc[[k]][[c]]
            return(lfc)
        }), vector("list", length(cats))) %>% 
        set_rownames(cats)
    
    # compute NB parameters
    b <- exp(rowData(x)$beta)
    o <- exp(colData(x)$offset)
    m <- vapply(o, function(l) b*l, numeric(nrow(x)))
    dimnames(m) <- dimnames(x)
    d <- rowData(x)$dispersion
    names(d) <- rownames(x)
    
    for (k in kids) {
        for (s in sids) {
            for (c in cats[n_dd[, k] != 0]) {
                gs_kc <- gs_by_kc[[k]][[c]]
                cs_ks <- cs_by_ks[[k]][[s]]
                
                g1 <- cs_idx[[k]][[s]]$A
                g2 <- cs_idx[[k]][[s]]$B
                
                ng1 <- length(g1)
                ng2 <- length(g2) 
                
                cs_g1 <- sample(cs_ks, ng1, replace = TRUE)
                cs_g2 <- sample(cs_ks, ng2, replace = TRUE)
                
                m_g1 <- m[gs_kc, cs_g1, drop = FALSE]
                m_g2 <- m[gs_kc, cs_g2, drop = FALSE]
                d_kc <- d[gs_kc]
                lfc_kc <- lfc[[c, k]]
                
                counts <- .sim(c, cs_g1, cs_g2, m_g1, m_g2, d = d_kc, lfc = lfc_kc)
                y[gs_idx[[c, k]], c(g1, g2)] <- counts
            }
        }
    }
    
    # construct gene metadata table storing
    # gene | cluster_id | category | logFC
    gi <- data.frame(
        gene = unlist(gs_idx),
        cluster_id = rep.int(rep(kids, each = length(cats)), c(n_dd)),
        category = rep.int(rep(cats, nk), c(n_dd)),
        logFC = unlist(lfc)) %>% 
        mutate_at("gene", as.character)
    o <- order(as.numeric(gsub("[a-z]", "", gi$gene)))
    gi <- gi[o, ] %>% set_rownames(NULL)
    
    # construct SCE
    cd$sample_id <- factor(paste(cd$sample_id, cd$group_id, sep = "."))
    m <- match(levels(cd$sample_id), cd$sample_id)
    gids <- cd$group_id[m]
    o <- order(gids)
    sids <- levels(cd$sample_id)[o]
    ei <- data.frame(sample_id = sids, group_id = gids[o])
    cd <- cd %>% mutate_at("sample_id", factor, levels = sids)
    
    md <- list(
        experiment_info = ei,
        n_cells = table(cd$sample_id),
        gene_info = gi)
    
    SingleCellExperiment(
        assays = list(counts = as.matrix(y)),
        colData = cd, 
        metadata = md)
}


















