#' simData
#' 
#' Simulation of complex scRNA-seq data 
#' 
#' \code{simData} simulates multiple clusters and samples 
#' across 2 experimental conditions from a real scRNA-seq data set.
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param ng,nc,ns,nk # of genes, cells, samples and clusters to simulate. 
#' @param probs a list of length 3 containing probabilities of a cell belonging
#'   to each cluster, sample, and group, respectively. List elements must be 
#'   NULL (equal probabilities) or numeric values in [0, 1] that sum to 1.
#' @param p_dd numeric vector of length 6.
#'   Specifies the probability of a gene being
#'   EE, EP, DE, DP, DM, or DB, respectively.
#' @param paired logial specifying whether a paired design should 
#'   be simulated (both groups use the same set of reference samples) 
#'   or not (reference samples are drawn at random).
#' @param p_ep,p_dp,p_dm numeric specifying the proportion of cells
#'   to be shifted to a different expression state in one group (see details).
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
#' library(SingleCellExperiment)
#' 
#' # prep. SCE for simulation
#' sce <- prepSim(sce)
#' 
#' # simulate data
#' (sim <- simData(sce, ng = 100, nc = 10,
#'   p_dd = c(0.9, 0, 0.1, 0, 0, 0)))
#' 
#' # simulation metadata
#' head(gi <- metadata(sim)$gene_info)
#' 
#' # should be ~10% DE  
#' table(gi$category)
#' 
#' # unbalanced sample sizes
#' sim <- simData(sce, ng = 10, nc = 100,
#'   probs = list(NULL, c(0.25, 0.75), NULL))
#' table(sim$sample_id)
#' 
#' # one group only
#' sim <- simData(sce, ng = 10, nc = 100,
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
#' @importFrom purrr modify_depth set_names
#' @importFrom stats model.matrix rgamma setNames
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors split
#' @export

simData <- function(x, ng = nrow(x), nc = 2e3, ns = 3, nk = 3,
    probs = NULL, p_dd = diag(6)[1, ], paired = FALSE,
    p_ep = 0.5, p_dp = 0.3, p_dm = 0.5,
    p_type = 0, lfc = 2, rel_lfc = NULL) {
    
    # throughout this code...
    # k: cluster ID
    # s: sample ID
    # g: group ID
    # c: DD category
    # 0: reference
    
    # check validity of input arguments
    .check_sce(x, req_group = FALSE)
    .check_args_simData(as.list(environment()))
    
    # reference IDs
    nk0 <- length(kids0 <- set_names(levels(x$cluster_id)))
    ns0 <- length(sids0 <- set_names(levels(x$sample_id)))
    
    # simulation IDs
    nk <- length(kids <- set_names(paste0("cluster", seq_len(nk))))
    sids <- set_names(paste0("sample", seq_len(ns)))
    gids <- set_names(c("A", "B"))
    
    # sample reference clusters & samples
    ref_kids <- setNames(sample(kids0, nk, nk > nk0), kids)
    if (paired) { 
        # use same set of reference samples for both groups
        ref_sids <- sample(sids0, ns, ns > ns0)
        ref_sids <- replicate(length(gids), ref_sids)
    } else {
        # draw reference samples at random for each group
        ref_sids <- replicate(length(gids), 
            sample(sids0, ns, ns > ns0))
    }
    dimnames(ref_sids) <- list(sids, gids)
    
    if (is.null(rel_lfc)) 
        rel_lfc <- rep(1, nk)
    if (is.null(names(rel_lfc))) {
        names(rel_lfc) <- kids
    } else {
        stopifnot(names(rel_lfc) %in% kids0)
    }
    
    # initialize count matrix
    gs <- paste0("gene", seq_len(ng))
    cs <- paste0("cell", seq_len(nc))
    y <- matrix(0, ng, nc, dimnames = list(gs, cs))
    
    # sample cell metadata
    cd <- .sample_cell_md(
        n = nc, probs = probs,
        ids = list(kids, sids, gids))
    rownames(cd) <- cs
    cs_idx <- .split_cells(cd, by = colnames(cd))
    n_cs <- modify_depth(cs_idx, -1, length)
    
    # split input cells by cluster-sample
    cs_by_ks <- .split_cells(x)
    
    # sample nb. of genes to simulate per category & gene indices
    n_dd <- table(sample(cats, ng, TRUE, p_dd))
    n_dd <- replicate(nk, n_dd)
    colnames(n_dd) <- kids
    gs_idx <- .sample_gene_inds(gs, n_dd)
    
    # for ea. cluster, sample set of genes to simulate from
    gs_by_k <- setNames(sample(rownames(x), ng, TRUE), gs)
    gs_by_k <- set_colnames(replicate(nk, gs_by_k), kids)

    # impute type-genes
    if (p_type != 0)
        gs_by_k <- .impute_type_genes(x, gs_by_k, gs_idx, p_type)

    # split by cluster & categroy
    gs_by_k <- split(gs_by_k, col(gs_by_k))
    gs_by_k <- setNames(map(gs_by_k, set_names, gs), kids)
    
    gs_by_kc <- lapply(kids, function(k) 
        lapply(unfactor(cats), function(c) 
            gs_by_k[[k]][gs_idx[[c, k]]])) 
    
    # sample logFCs
    lfc <- vapply(kids, function(k) 
        lapply(unfactor(cats), function(c) { 
            n <- n_dd[c, k]
            if (c == "ee") return(rep(NA, n))
            signs <- sample(c(-1, 1), n, TRUE)
            lfcs <- rgamma(n, 4, 4/lfc) * signs
            names(lfcs) <- gs_by_kc[[k]][[c]]
            lfcs * rel_lfc[k]
        }), vector("list", length(cats)))

    # compute NB parameters
    m <- lapply(sids0, function(s) {
        b <- paste0("beta.", s)
        b <- exp(rowData(x)[[b]])
        m <- outer(b, exp(x$offset), "*")
        dimnames(m) <- dimnames(x); m
    })
    d <- rowData(x)$dispersion 
    names(d) <- rownames(x)
    
    sim_mean <- lapply(kids, function(k) 
        lapply(gids, function(g) 
            setNames(numeric(ng), gs)))
    
    for (k in kids) {
        for (s in sids) {
            # get reference samples, clusters & cells
            s0 <- ref_sids[s, ]
            k0 <- ref_kids[k]
            cs0 <- cs_by_ks[[k0]][s0]
            
            # get output cell indices
            ci <- unlist(cs_idx[[k]][[s]])
            
            for (c in cats[n_dd[, k] != 0]) {
                # sample cells to simulate from
                cs_g1 <- sample(cs0[[1]], n_cs[[k]][[s]][[1]], TRUE)
                cs_g2 <- sample(cs0[[2]], n_cs[[k]][[s]][[2]], TRUE)
                
                # get reference genes & output gene indices
                gs0 <- gs_by_kc[[k]][[c]] 
                gi <- gs_idx[[c, k]]
                
                # get NB parameters
                m_g1 <- m[[s0[[1]]]][gs0, cs_g1, drop = FALSE]
                m_g2 <- m[[s0[[2]]]][gs0, cs_g2, drop = FALSE]
                d_kc <- d[gs0]
                lfc_kc <- lfc[[c, k]]
                
                re <- .sim(c, cs_g1, cs_g2, m_g1, m_g2, d_kc, lfc_kc, p_ep, p_dp, p_dm)
                y[gi, ci] <- re$cs
                
                for (g in gids) sim_mean[[k]][[g]][gi] <- ifelse(
                    is.null(re$ms[[g]]), NA, list(re$ms[[g]]))[[1]]
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
    gi <- set_rownames(gi[o, ], NULL)
    
    # construct SCE
    cd$group_id <- droplevels(cd$group_id)
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
        gene_info = gi,
        ref_sids = ref_sids,
        ref_kids = ref_kids)
    
    SingleCellExperiment(
        assays = list(counts = as.matrix(y)),
        colData = cd, metadata = md)
}
