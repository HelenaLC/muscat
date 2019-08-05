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
#' @param probs a list of length 3 containing probabilities of a cell belonging
#'   to each cluster, sample, and group, respectively. List elements must be 
#'   NULL (equal probabilities) or numeric values in [0, 1] that sum to 1.
#' @param p_dd numeric vector of length 6.
#'   Specifies the probability of a gene being
#'   EE, EP, DE, DP, DM, or DB, respectively.
#' @param p_type numeric. Probaility of EE/EP gene being a type-gene.
#'   If a gene is of class "type" in a given cluster, a unique mean 
#'   will be used for that gene in the respective cluster.
#' @param lfc numeric value to use as mean logFC
#'   for DE, DP, DM, and DB type of genes.
#' @param rel_lfc numeric vector of relative logFCs for each cluster. 
#'   Should be of length \code{nlevels(x$cluster_id)} with 
#'   \code{levels(x$cluster_id)} as names. 
#'   Defaults to factor of 1 for all clusters.
#'   
#' @return a \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   containing multiple clusters & samples across 2 groups.
#' 
#' @examples
#' data(sce)
#' 
#' # prep. SCE for simulation
#' sce <- prepSim(sce)
#' 
#' # simulate data
#' (sim <- simData(sce,
#'   n_genes = 100, n_cells = 10,
#'   p_dd = c(0.9, 0, 0.1, 0, 0, 0)))
#' 
#' # simulation metadata
#' head(gi <- metadata(sim)$gene_info)
#' 
#' # should be ~10% DE  
#' table(gi$category)
#' 
#' # unbalanced sample sizes
#' sim <- simData(sce,
#'   n_genes = 10, n_cells = 100,
#'   probs = list(NULL, c(0.1, 0.3, 0.6), NULL))
#' table(sim$sample_id)
#' 
#' # one group only
#' sim <- simData(sce,
#'   n_genes = 10, n_cells = 100,
#'   probs = list(NULL, NULL, c(1, 0)))
#' levels(sim$group_id)
#'     
#' @author Helena L Crowell
#' 
#' @references 
#' Crowell, HL, Soneson, C, Germain, P-L, Calini, D, 
#' Collin, L, Raposo, C, Malhotra, D & Robinson, MD: 
#' On the discovery of population-specific state transitions from 
#' multi-sample multi-condition single-cell RNA sequencing data. 
#' \emph{bioRxiv} \strong{713412} (2018). 
#' doi: \url{https://doi.org/10.1101/713412}
#' 
#' @importFrom data.table data.table
#' @importFrom dplyr mutate_all mutate_at
#' @importFrom edgeR DGEList estimateDisp glmFit
#' @importFrom magrittr set_colnames set_rownames
#' @importFrom purrr modify_at set_names
#' @importFrom stats model.matrix rgamma setNames
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors split
#' @export

simData <- function(x, n_genes = 500, n_cells = 300, 
    probs = NULL, p_dd = diag(6)[1, ], p_type = 0,
    lfc = 2, rel_lfc = NULL) {
    
    # throughout this code...
    # k: cluster ID
    # s: sample ID
    # g: group ID
    # c: DD category
    
    # check validity of input arguments
    .check_sce(x, req_group = FALSE)
    .check_args_simData(as.list(environment()))
    
    kids <- set_names(levels(x$cluster_id))
    sids <- set_names(levels(x$sample_id))
    gids <- set_names(c("A", "B"))
    nk <- length(kids)
    
    if (is.null(rel_lfc)) 
        rel_lfc <- rep(1, nk)
    if (is.null(names(rel_lfc))) 
        names(rel_lfc) <- kids
    
    # initialize count matrix
    gs <- paste0("gene", seq_len(n_genes))
    cs <- paste0("cell", seq_len(n_cells))
    y <- matrix(0, n_genes, n_cells, dimnames = list(gs, cs))
    
    # sample cell metadata
    cd <- .sample_cell_md(
        n = n_cells, probs = probs,
        ids = list(kids, sids, gids)) %>% 
        set_rownames(cs)
    cs_idx <- .split_cells(cd, by = colnames(cd))
    n_cs <- modify_depth(cs_idx, -1, length)
    
    # split input cells by cluster-sample
    cs_by_ks <- .split_cells(x)
    
    # sample nb. of genes to simulate per category & gene indices
    n_dd <- replicate(nk, 
        table(sample(factor(cats, levels = cats), n_genes, TRUE, p_dd))) %>% 
        set_colnames(kids)
    gs_idx <- .sample_gene_inds(gs, n_dd)
    
    # for ea. cluster, sample set of genes to simulate from
    gs_by_k <- setNames(sample(rownames(x), n_genes, TRUE), gs)
    gs_by_k <- replicate(nk, gs_by_k) %>% set_colnames(kids)

    # impute type-genes
    if (p_type != 0)
        gs_by_k <- .impute_type_genes(x, gs_by_k, gs_idx, p_type)

    # split by cluster & categroy
    gs_by_k <- split(gs_by_k, col(gs_by_k))
    gs_by_k <- map(gs_by_k, set_names, gs)
    names(gs_by_k) <- kids
    
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
            lfcs <- rgamma(n, 4, 4/lfc) * signs
            names(lfcs) <- gs_by_kc[[k]][[c]]
            lfcs * rel_lfc[k]
        }), vector("list", length(cats))) %>% 
        set_rownames(cats)
    
    # compute NB parameters
    o <- exp(colData(x)$offset)
    m <- lapply(sids, function(s) {
        cn <- paste("beta", s, sep = ".")
        k <- grep(cn, names(rowData(x)))
        b <- exp(rowData(x)[[k]])
        vapply(o, "*", b, FUN.VALUE = numeric(nrow(x))) %>% 
            set_rownames(rownames(x)) %>% 
            set_colnames(colnames(x)) %>% 
            round
    })
    d <- rowData(x)$dispersion %>% 
        set_names(rownames(x))
    
    sim_mean <- lapply(kids, function(k) 
        lapply(gids, function(g)
            setNames(numeric(n_genes), rownames(y))))
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
                
                m_g1 <- m[[s]][gs_kc, cs_g1, drop = FALSE]
                m_g2 <- m[[s]][gs_kc, cs_g2, drop = FALSE]
                d_kc <- d[gs_kc]
                lfc_kc <- lfc[[c, k]]
                
                gidx <- gs_idx[[c, k]]
                cidx <- c(g1, g2)
                
                re <- .sim(c, cs_g1, cs_g2, m_g1, m_g2, d_kc, lfc_kc)
                y[gidx, cidx] <- re$cs
                
                for (g in c("A", "B")) sim_mean[[k]][[g]][gidx] <- 
                    ifelse(is.null(re$ms[[g]]), NA, list(re$ms[[g]]))[[1]]
            }
        }
    }
    
    sim_mean <- sim_mean %>%
        map(bind_cols) %>% 
        bind_rows(.id = "cluster_id") %>% 
        mutate_at("cluster_id", factor) %>% 
        mutate(gene = rep(gs, nk))
    
    # construct gene metadata table storing ------------------------------------
    # gene | cluster_id | category | logFC, gene, disp, mean used for sim.
    gi <- data.frame(
        gene = unlist(gs_idx),
        cluster_id = rep.int(rep(kids, each = length(cats)), c(n_dd)),
        category = rep.int(rep(cats, nk), c(n_dd)),
        logFC = unlist(lfc),
        sim_gene = unlist(gs_by_kc),
        sim_disp = d[unlist(gs_by_kc)]) %>% 
        mutate_at("gene", as.character)
    # add true simulation means
    gi <- full_join(gi, sim_mean, by = c("gene", "cluster_id")) %>% 
        rename("sim_mean.A" = "A", "sim_mean.B" = "B")
    # reorder
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
