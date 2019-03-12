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
cats <- c("ee", "ep", "de", "dp", "dm", "db")
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
#   n:     nb. of cells
#   ids:   list of IDs to sample from
#   probs: list of probabilities for ea. set of IDs
# ------------------------------------------------------------------------------
#' @importFrom magrittr set_colnames
.sample_cell_md <- function(n, ids, probs = NULL) {
    ns <- vapply(ids, length, numeric(1))
    if (is.null(probs)) 
        probs <- vector("list", 3)
    probs <- lapply(1:3, function(i) {
        if (!is.null(probs[[i]])) {
            return(probs[[i]])
        } else {
            rep(1 / ns[i], ns[i])
        }
    })
    vapply(1:3, function(i) 
        sample(ids[[i]], n, TRUE, probs[[i]]), 
        character(n)) %>% data.frame(row.names = NULL) %>% 
        set_colnames(c("cluster_id", "sample_id", "group_id"))
}

# ------------------------------------------------------------------------------
# for ea. gene category, sample a non-intersecting set of gene indices
# ------------------------------------------------------------------------------
.sample_gene_inds <- function(gs, ns) {
    cluster_ids <- colnames(ns)
    vapply(cluster_ids, function(k)
        split(sample(gs), rep.int(cats, ns[, k])),
        vector("list", length(cats)))
}

# ------------------------------------------------------------------------------
# helper to sample from a NB across a grid 
# of dispersions ('size') and means ('mu')
# ------------------------------------------------------------------------------
#' @importFrom stats rnbinom
.nb <- function(cs, d, m, lfc = NULL) {
    n_gs <- length(d)
    n_cs <- length(cs)
    if (is.null(lfc))
        lfc <- rep(0, n_gs)
    lfc[lfc < 0] <- 0
    fc <- 2 ^ lfc
    fc <- rep(fc, each = n_cs)
    nb <- rnbinom(n_gs * n_cs, 
        size = rep(1/d, each = n_cs), 
        mu = c(t(m[, cs])) * fc)
    matrix(nb, byrow = TRUE,
        nrow = n_gs, ncol = n_cs, 
        dimnames = list(names(d), cs))
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
    cs_g1, cs_g2, m_g1, m_g2, d, lfc) {
    
    ng1 <- length(cs_g1)
    ng2 <- length(cs_g2)
    
    switch(match.arg(cat),
        ee = {
            cbind(
                .nb(cs_g1, d, m_g1),
                .nb(cs_g2, d, m_g2))
        },
        ep = {
            g1_hi <- sample(ng1, round(ng1 * 0.5))
            g2_hi <- sample(ng2, round(ng2 * 0.5))
            cbind(
                .nb(cs_g1[-g1_hi], d, m_g1),
                .nb(cs_g2[-g2_hi], d, m_g2),
                .nb(cs_g1[ g1_hi], d, m_g1, -lfc), # lfc < 0 => 50% g2 hi
                .nb(cs_g2[ g2_hi], d, m_g2,  lfc)) # lfc > 0 => 50% g2 hi
        },
        de = {
            cbind(
                .nb(cs_g1, d, m_g1, -lfc), # lfc < 0 => all g1 hi
                .nb(cs_g2, d, m_g2,  lfc)) # lfc > 0 => all g2 hi
        },
        dp = {
            g1_hi <- sample(ng1, round(ng1 * 0.3))
            g2_hi <- sample(ng2, round(ng2 * 0.7))
            cbind(                                 # g %  >0 <0
                .nb(cs_g1[-g1_hi], d, m_g1, -lfc), # 1 70 -- up
                .nb(cs_g1[ g1_hi], d, m_g1,  lfc), # 1 30 up --
                .nb(cs_g2[-g2_hi], d, m_g2, -lfc), # 2 30 -- up
                .nb(cs_g2[ g2_hi], d, m_g2,  lfc)) # 2 70 up --
        },
        dm = {
            g1_hi <- sample(ng1, round(ng1 * 0.5))
            g2_hi <- sample(ng2, round(ng2 * 0.5))
            cbind(
                .nb(cs_g1[-g1_hi], d, m_g1),
                .nb(cs_g2[-g2_hi], d, m_g2),
                .nb(cs_g1[ g1_hi], d, m_g1, -lfc), # lfc < 0 => 50% g1 hi
                .nb(cs_g2[ g2_hi], d, m_g2,  lfc)) # lfc > 0 => 50% g2 hi
        }, 
        db = {
            g2_hi <- sample(ng2, round(ng2 * 0.5))
            cbind(
                .nb(cs_g1, d, m_g1, lfc/2),       # all g1 mi
                .nb(cs_g2[-g2_hi], d, m_g2),      # 50% g2 lo
                .nb(cs_g2[ g2_hi], d, m_g2, lfc)) # 50% g2 hi
        }
    )
}