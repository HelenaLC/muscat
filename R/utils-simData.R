# ------------------------------------------------------------------------------
# gene categories for simulation of differential distributions
# ------------------------------------------------------------------------------
#   ee = equivalently expressed
#   ep = equivalent proportions
#   de = 'classical' differential expression (shift in means)
#   dp = differential proportions
#   dm = differential modality
#   db = shift in mean & dm
# ------------------------------------------------------------------------------
names(cats) <- cats <- c("ee", "ep", "de", "dp", "dm", "db")
cats <- factor(cats, levels = cats)

# ------------------------------------------------------------------------------
# for ea. cluster-sample, samples the number of cells for 2 groups
# ------------------------------------------------------------------------------
#   n: single numeric or range to sample from
#   k: character vector of cluster IDs
#   s: sample IDs
# > array of dim. #s x #k x 2 with row names = sample IDs & 
#   column names = cluster IDs. Ea. entry is a list of length 2 
#   with numberic values equal to or in the range of n.
# ------------------------------------------------------------------------------
.sample_n_cells <- function(n, k, s) {
    nk <- length(k)
    ns <- length(s)
    if (length(n) == 1) {
        n <- list(rep(n, 2))
    } else {
        n <- replicate(nk * ns, 
            list(sample(seq(n[1], n[2]), 2)))
    }
    matrix(n, 
        nrow = ns, ncol = nk, 
        dimnames = list(s, k))
}

# ------------------------------------------------------------------------------
# generate a randomized data.frame/colData of cell metadata
# (cluster IDs, sample IDs, and group IDs)
# ------------------------------------------------------------------------------
#   n     = nb. of cells
#   ids   = list of IDs to sample from
#   probs = list of probabilities for ea. set of IDs
#           (in order of cluster, sample, group)
# ------------------------------------------------------------------------------
.sample_cell_md <- function(n, ids, probs = NULL) {
    ns <- vapply(ids, length, numeric(1))
    if (is.null(probs)) 
        probs <- vector("list", 3)
    probs <- lapply(seq_along(probs), function(i) {
        if (!is.null(probs[[i]])) {
            return(probs[[i]])
        } else {
            rep(1 / ns[i], ns[i])
        }
    })
    cd <- vapply(seq_along(probs), function(i) 
        sample(ids[[i]], n, TRUE, probs[[i]]), 
        character(n))
    cd <- data.frame(cd, row.names = NULL)
    colnames(cd) <- c("cluster_id", "sample_id", "group_id")
    cd$group_id <- factor(cd$group_id, levels = ids[[3]])
    return(cd)
}

# ------------------------------------------------------------------------------
# get gene indeces for by gene category & cluster
#   gs = character vector of gene names
#   ns = nb. of genes in ea. category
#   > array of dim. #(categories) x #(clusters);
#     ea. entry is a character vector of genes
# ------------------------------------------------------------------------------
.sample_gene_inds <- function(gs, ns) {
    cluster_ids <- colnames(ns)
    vapply(cluster_ids, function(k)
        split(sample(gs), rep.int(cats, ns[, k])),
        vector("list", length(cats)))
}

# helper to extract cluster IDs as "clusterN" from phylogeny
.get_clusters_from_phylo <- function(u) {
    pat <- "(?<=').*?(?=')"
    v <- gregexpr(pat, u, perl = TRUE)
    u <- unlist(regmatches(u, v))
    grep("cluster[0-9]+", u, value = TRUE)
}

# ------------------------------------------------------------------------------
#  Read the nodes of a phylogram to create a table of state/ type 
#  that corresponds to the relations between clusters. 
#  Recursively calls itself at each node to update the class_tbl with shared 
#  category among the relevant clusters. 
# ------------------------------------------------------------------------------
#   phylo_tree   = input cell phylogram (see '?simData')
#   class_tbl    = class of gene category (type or state)
#   used_tg      = genes already used as 'type' in previous recursions
#   phylo_pars   = distance parameters for the number of 
#                  genes shared b/w branches (see '?simData')
#   > i) updated 'class_tbl' for the current node (if in recursion) or 
#     updated 'class_tbl' for the whole tree (if all nodes were read);
#     ii) genes already used as 'shared' in previous recursions
# ------------------------------------------------------------------------------
#' @importFrom dplyr %>%
.read_branch <- function(phylo_tree, class_tbl, used, phylo_pars) {
    # assure there's no linebreaks
    phylo <- gsub("\n", "", phylo_tree)
    phygs <- gsub("^\\(|\\);$", "", phylo) 
    # decompose groups from phylogram
    phygs <- strsplit(phygs, "\\([^)]+,(*SKIP)(*FAIL)|,\\s*", perl = TRUE)[[1]]
    # grep distances & remove them from groups
    ds <- lapply(phygs, function(u) as.numeric(gsub(".*:", "", u)))
    phygs <- lapply(phygs, function(u) gsub("\\:[0-9]*\\.[0-9]*$", ";", u))
    # identify clusters to assign type genes
    k_shared <- lapply(phygs, .get_clusters_from_phylo)
    # compute number of shared genes b/w these clusters as 
    # Exp w/ intercept nb. genes x theta1 & rate distance x theta2
    ng <- nrow(class_tbl)
    n_shared <- lapply(ds, function(d) 
        ceiling(phylo_pars[1]*ng*exp(-phylo_pars[2]*d)))
    if (length(used) != 0) {
        gs <- rownames(class_tbl)
        gs <- gs[!gs %in% used]
    } else gs <- rownames(class_tbl)
    if (sum(unlist((n_shared))) > length(gs)) stop(
        "Ran out of genes to sample from;\n  ",
        "please simulate more genes or adjust 'phylo_pars'.")
    type_gs <- list()
    for (i in seq_along(n_shared)) {
        type_gs[[i]] <- sample(gs, n_shared[[i]])
        used <- c(used, type_gs[[i]])
        gs <- setdiff(gs, used)
    }
    # update class table
    for (i in seq_along(type_gs))
        class_tbl[type_gs[[i]], k_shared[[i]]] <- "type" 
    # remove phylogeny groups that reached leaf (recognized by missing ",")
    is_not_leaf <- vapply(phygs, function(u) 
        length(grep("\\,", u)) > 0, logical(1))
    phygs <- phygs[is_not_leaf]
    # stop if no nodes left, otherwise recursion on further nodes
    if (length(phygs) != 0) 
        for (node in phygs) {
            res <- .read_branch(node, class_tbl, used, phylo_pars)
            class_tbl <- res$class_tbl; used <- res$used
        }
    list(class_tbl = class_tbl, used = used)
}

# ------------------------------------------------------------------------------
# sample marker classes based on a cell phylogram by calling .read_branch
# ------------------------------------------------------------------------------
#   x           = input SCE
#   gs_by_k     = n_genes x n_clusters matrix of 'x' genes to use for sim.
#   gs_idx      = n_category x n_clusters matrix of output gene indices
#   phylo_tree  = input cell phylogram (see also '?simData')
#   phylo_pars  = parameters to define the number of shared genes, 
#                 based on branch distance (see also '?simData')
#   positive_type = force type genes to be selected from ref genes that have 
#   > returns a list of 
#   1 an updated marker class matrix
#   2 a vector of genes already used as 'type' 
#     to avoid re-use by '.impute_type_genes()'
#   3 a gene information matrix to be returned as part of the
#     simulation metadata stored in `metadata(x)$gene_info`
# ------------------------------------------------------------------------------
.impute_shared_type_genes <- function(x, gs_by_k, gs_idx, phylo_tree, 
                                      phylo_pars, positive_type) {
    # sample gene-classes for genes of categroy EE & EP
    ex_cats <- c("ee", "ep")
    n_not_de <- sum(vapply(gs_idx[ex_cats, 1], length, numeric(1)))
    not_de <- sample(rownames(gs_by_k), n_not_de)
    names(kids) <- kids <- colnames(gs_idx)
    # initialize temporary copy of 'gs_by_k' 
    # such that all genes are of class "state"
    class_tbl <- gs_by_k
    ij <- lapply(dim(gs_by_k), seq_len)
    class_tbl[ij[[1]], ij[[2]]] <- "state"
    # track type-genes already used while looping
    res <- .read_branch(phylo_tree, class_tbl[not_de, ], used = c(), phylo_pars)
    class_tbl[not_de, ] <- res$class_tbl[not_de, ]
    used <- res$used
    # sample genes that were defined as type 
    # (same across related clusters)
    # for (g in used) {
    #     k <- kids[class_tbl[g, ] == "type"]
    #     gs <- setdiff(rownames(x), gs_by_k[g, !kids %in% k])
    #     gs_by_k[g, k] <- sample(gs, 1)
    # }
    
    mean_beta <- apply(rowData(x)[, grep("beta\\.[0-9]*", colnames(rowData(x)))], 
                       1, mean)
    for (g in used) {
        k <- kids[class_tbl[g, ] == "type"]
        g_beta <- mean_beta[names(mean_beta) %in% gs_by_k[g,]]
        if (length(positive_type) == 0) {
          gs_sel <- rownames(x)
        } else if (positive_type == TRUE) {
          gs_sel <- rownames(x)[mean_beta > max(g_beta)]
        } else if (positive_type == FALSE) {
          gs_sel <- rownames(x)[mean_beta < min(g_beta)]
        } 
        if (length(gs_sel) == 0) {
          disp <- rowData(x)$dispersion
          g_disp <- unique(rowData(x)[gs_by_k[g,], "dispersion"])
          if (positive_type == TRUE) {
            gs_sel <- rownames(x)[g_disp < disp]
          } else {
            gs_sel <- rownames(x)[g_disp > disp]
          }
          if (length(length(gs_sel)) == 0) {
            # take the second max/min if a smaller beta can not be found
            gs_sel <- names(sort(mean_beta, decreasing = positive_type)[2])
          }
        }
        gs <- setdiff(gs_sel, gs_by_k[g, !kids %in% k])
        gs_by_k[g, k] <- sample(gs, 1)
    }
    
    # split classes by gene
    cs_by_g <- split(class_tbl, row(class_tbl))
    names(cs_by_g) <- rownames(gs_by_k)
    # get gene specificities
    specs <- lapply(cs_by_g, function(u) {
        if (all(u == "state")) return(NA)
        unname(kids[u == "type"])
    })
    # get gene classes
    class <- vapply(specs, function(u)
        ifelse(isTRUE(is.na(u)), "state", "shared"),
        character(1))
    # get single-type genes
    is_type <- unlist(lapply(cs_by_g, function(u) sum(u == "type") == 1))
    if (!all(!is_type)) class[is_type] <- "type"
    list(gs_by_k = gs_by_k, used = used, class = class, specs = specs)
}

# ------------------------------------------------------------------------------
# for ea. cluster, sample marker classes
#   x       = input 'SingleCellExperiment'
#   gs_by_k = n_genes x n_clusters matrix of 'x' genes to use for sim
#   gs_idx  = n_category x n_clusters matrix of output gene indices
#   p_type  = prob. of EE/EP gene being of class "type"
#   positive_type = force type genes to be selected from ref genes that have 
#   a higher beta than the genes from other clusters. 
#   > type-genes may only be of categroy EE & EP type of genes,
#     and use a cluster-specific mean in the NB count simulation
# ------------------------------------------------------------------------------
#' @importFrom data.table data.table
#' @importFrom dplyr %>%
#' @importFrom purrr map
.impute_type_genes <- function(x, gs_by_k, gs_idx, p_type, positive_type) {
    names(kids) <- kids <- colnames(gs_idx)
    if (length(p_type) == 1) {
        p_type <- rep(p_type, ncol(gs_by_k)) 
        names(p_type) <- colnames(gs_by_k)
    }
    # sample gene-classes for genes of categroy EE & EP
    non_de <- c("ee", "ep")
    
    class_tbl <- list()
    used_gs <- c()
    for(k in kids){
      gs <- unlist(gs_idx[non_de, k])
      gs0 <- gs[gs %in% used_gs]
      gs <- gs[!gs %in% used_gs]
      n <- length(gs)
      # avoid re-sampling of the same genes that would be specific to multiple 
      # clusters
      if (k == kids[1]){
        p_type_k <- p_type[k]
      } else {
        p_type_k <- p_type[k]/(1-(length(used_gs)/nrow(x)))
        if (p_type_k > 1) stop("Not enough genes to simulate. Please reduce p_type.")
      }
      class_tbl[[k]] <- data.table(
        stringsAsFactors = FALSE,
        gene = gs, cluster_id = k,
        class = sample(factor(c("state", "type")), n,
                       prob = c(1 - p_type_k, p_type_k), replace = TRUE))
      if (k != kids[1]){
        class_tbl[[k]] <- rbind(class_tbl[[k]], 
                                data.table(
                                   stringsAsFactors = FALSE,
                                   gene = gs0, cluster_id = k,
                                   class = "state")) 
        class_tbl[[k]] <- class_tbl[[k]][match(unlist(gs_idx[non_de, k]), 
                                               class_tbl[[k]]$gene), ]
      }
      used_gs <- c(used_gs, class_tbl[[k]]$gene[class_tbl[[k]]$class == "type"])
    }
    class_tbl <- class_tbl %>% map(split, by = "class", flatten = FALSE)
    
    # class_tbl <- lapply(kids, function(k) {
    #     gs <- unlist(gs_idx[non_de, k])
    #     n <- length(gs)
    #     data.table(
    #         stringsAsFactors = FALSE,
    #         gene = gs, cluster_id = k,
    #         class = sample(factor(c("state", "type")), n,
    #             prob = c(1 - p_type[k], p_type[k]), replace = TRUE))
    # }) %>% map(split, by = "class", flatten = FALSE)
    mean_beta <- apply(rowData(x)[, grep("beta\\.[0-9]*", colnames(rowData(x)))], 
                       1, mean)
    for (k in kids) {
      type_gs <- class_tbl[[k]]$type$gene
      gs_by_k[type_gs, k] <- apply(
        gs_by_k[type_gs, kids != k, drop = FALSE], 1,
        function(ex){
          # only sample from genes that have lower/ greater beta to control
          # for up/down marker genes
          if (length(positive_type) == 0) {
            gs_sel <- rownames(x)
          } else if (positive_type == TRUE) {
            gs_sel <- rownames(x)[mean_beta > max(mean_beta[ex])]
          } else if (positive_type == FALSE) {
            gs_sel <- rownames(x)[mean_beta < min(mean_beta[ex])]
          } 
          if (length(gs_sel) == 0) {
            disp <- rowData(x)$dispersion
            g_disp <- unique(rowData(x)[ex, "dispersion"])
            if (positive_type == TRUE) {
              gs_sel <- rownames(x)[g_disp < disp]
            } else {
              gs_sel <- rownames(x)[g_disp > disp]
            }
            if (length(length(gs_sel)) == 0) {
              # take the second max/min if a smaller beta can not be found
              gs_sel <- names(sort(mean_beta, decreasing = positive_type)[2])
            }
          }
          sample(setdiff(gs_sel, ex), 1)
        } )
    }
    # # sample cluster-specific genes for ea. cluster & type-gene
    # for (k in kids) {
    #     type_gs <- class_tbl[[k]]$type$gene
    #     gs_by_k[type_gs, k] <- apply(
    #         gs_by_k[type_gs, kids != k, drop = FALSE], 1,
    #         function(ex) sample(setdiff(rownames(x), ex), 1))
    # }
    ng <- nrow(gs_by_k)
    gs <- rownames(gs_by_k)
    is_type <- map(class_tbl, "type")
    type_gs <- unlist(map(is_type, "gene"))
    shared_gs <- setdiff(gs, unlist(gs_idx))
    stopifnot(!any(shared_gs %in% type_gs))
    # get gene classes
    class <- rep("state", ng)
    names(class) <- gs
    class[shared_gs] <- "shared"
    class[unique(type_gs)] <- "type" 
    # get gene specificities
    specs <- rep(NA, ng)
    names(specs) <- gs
    ns <- vapply(is_type, nrow, numeric(1))
    specs[sort(unique(type_gs))] <- split(rep.int(kids, ns), type_gs)
    return(list(gs_by_k = gs_by_k, class = class, specs = specs))
}

# ------------------------------------------------------------------------------
# helper to sample from a NB across a grid 
# of dispersions ('size') and means ('mu')
# ------------------------------------------------------------------------------
#' @importFrom stats rnbinom
.nb <- function(cs, d, m, lfc = NULL, f = 1) {
    n_gs <- length(d)
    n_cs <- length(cs)
    if (is.null(lfc)) {
        lfc <- rep(0, n_gs)
    } else {
        lfc[lfc < 0] <- 0
    }
    fc <- f * (2 ^ lfc)
    fc <- rep(fc, each = n_cs)
    ds <- rep(1/d, each = n_cs)
    ms <- c(t(m[, cs])) * fc 
    y <- rnbinom(n_gs * n_cs, size = ds, mu = ms)
    y <- matrix(y, byrow = TRUE, 
        nrow = n_gs, ncol = n_cs, 
        dimnames = list(names(d), cs))
    ms <- split(ms, rep(seq_len(nrow(m)), each = n_cs))
    list(counts = y, means = ms)
}

# ------------------------------------------------------------------------------
# helper to simulate differential distributions 
# ------------------------------------------------------------------------------
#   gs:  character vector of genes to simulate from
#   cs:  character vector of cells to simulate from
#   ng1: nb. of cells in 1st group
#   ng2: nb. of cells in 2nd group
#   m:   G x C matrix of NB means
#   d:   numeric vector of dispersions
#   lfc: numeric vector of logFCs
# ------------------------------------------------------------------------------
.sim <- function(
    cat = c("ee", "ep", "de", "dp", "dm", "db"),
    cs_g1, cs_g2, m_g1, m_g2, d, lfc, ep, dp, dm) {
    
    cat <- match.arg(cat)
    ng1 <- length(cs_g1)
    ng2 <- length(cs_g2)
    
    re <- switch(cat,
        "ee" = {
            list(
                .nb(cs_g1, d, m_g1),
                .nb(cs_g2, d, m_g2))
        },
        "ep" = {
            g1_hi <- sample(ng1, round(ng1 * ep))
            g2_hi <- sample(ng2, round(ng2 * ep))
            list(
                .nb(cs_g1[-g1_hi], d, m_g1),
                .nb(cs_g1[ g1_hi], d, m_g1, lfc), # 50% g1 hi
                .nb(cs_g2[-g2_hi], d, m_g2),
                .nb(cs_g2[ g2_hi], d, m_g2, lfc)) # 50% g2 hi
        },
        "de" = {
            list(
                .nb(cs_g1, d, m_g1, -lfc), # lfc < 0 => all g1 hi
                .nb(cs_g2, d, m_g2,  lfc)) # lfc > 0 => all g2 hi
        },
        "dp" = {
            props <- sample(c(dp, 1 - dp), 2)
            g1_hi <- sample(ng1, round(ng1 * props[1]))
            g2_hi <- sample(ng2, round(ng2 * props[2]))
            list(                           
                .nb(cs_g1[-g1_hi], d, m_g1), 
                .nb(cs_g1[ g1_hi], d, m_g1,  lfc), # lfc > 0 => dp/(1-dp)% up
                .nb(cs_g2[-g2_hi], d, m_g2), 
                .nb(cs_g2[ g2_hi], d, m_g2, -lfc)) # lfc < 0 => (1-dp)/dp% up
        },
        "dm" = {
            g1_hi <- sample(ng1, round(ng1 * dm))
            g2_hi <- sample(ng2, round(ng2 * dm))
            list(
                .nb(cs_g1[-g1_hi], d, m_g1),
                .nb(cs_g1[ g1_hi], d, m_g1, -lfc), # lfc < 0 => 50% g1 hi
                .nb(cs_g2[-g2_hi], d, m_g2),
                .nb(cs_g2[ g2_hi], d, m_g2,  lfc)) # lfc > 0 => 50% g2 hi
        }, 
        "db" = {
            if (sample(c(TRUE, FALSE), 1)) {
                # all g1 mi, 50% g2 hi
                g2_hi <- sample(ng2, round(ng2 * 0.5))
                list(
                    .nb(cs_g1, d, m_g1, abs(lfc), 0.5),
                    .nb(cs_g2[-g2_hi], d, m_g2, -lfc), 
                    .nb(cs_g2[ g2_hi], d, m_g2,  lfc)) 
            } else {
                # all g2 mi, 50% g1 hi
                g1_hi <- sample(ng1, round(ng1 * 0.5))
                list(
                    .nb(cs_g2, d, m_g2, abs(lfc), 0.5), 
                    .nb(cs_g1[-g1_hi], d, m_g1, -lfc),  
                    .nb(cs_g1[ g1_hi], d, m_g1,  lfc))  
            }
        }
    )
    cs <- map(re, "counts")
    cs <- do.call("cbind", cs)
    ms <- map(re, "means")
    rmv <- vapply(ms, is.null, logical(1))
    ms <- ms[!rmv] %>% 
        map_depth(2, mean) %>% 
        map_depth(1, unlist) %>% 
        data.frame %>% 
        as.matrix
    ms <- switch(cat, 
        ee = ms,
        de = ms,
        db = if (ng2 == 0) {
            as.matrix(
                ms[, 1])
        } else {
            cbind(
                ms[, 1],
                rowMeans(ms[, c(2, 3)]))
        }, if (ng2 == 0) {
            as.matrix(
                rowMeans(ms[, c(1, 2)]))
        } else {
            cbind(
                rowMeans(ms[, c(1, 2)]),
                rowMeans(ms[, c(3, 4)]))
        })
    ms <- split(ms, col(ms))
    names(ms) <- c("A", "B")[c(ng1, ng2) != 0]
    list(cs = cs, ms = ms)
}
