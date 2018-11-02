# ------------------------------------------------------------------------------
# helper to sample from a NB across a grid 
# of dispersions ('size') and means ('mu')
# ------------------------------------------------------------------------------

nb <- function(gs, cs, d, m, fc = 1) {
    t(sapply(gs, function(g) sapply(cs, function(c)
        rnbinom(1, size = d[g], mu = m[g, c] * fc))))
}

# ------------------------------------------------------------------------------
# helper to simulate differential distributions 
# ------------------------------------------------------------------------------
# gs = character vector of genes to simulate from
# cs = character vector of cells to simulate from
# ng1 = nb. of cells in group 1
# ng2 = nb. of cells in group 2
# m = matrix of mu's by genes x cells
# d = named numeric vector of dispersions for each gene
# ------------------------------------------------------------------------------

simdd <- function(
    category = c("ee", "ep", "de", "dp", "dm", "db"),
    gs, cs, ng1, ng2, m, d) {

    jg1 <- sample(cs, ng1, replace = TRUE)
    jg2 <- sample(cs, ng2, replace = TRUE)
    
    switch (match.arg(category),
        ee = {
            nb(gs, c(jg1, jg2), d, m)
        },
        ep = {
            g1_lo <- sample(ng1, round(ng1 * 0.6))
            g2_lo <- sample(ng2, round(ng2 * 0.6))
            cbind(
                nb(gs, jg1[ g1_lo], d, m),
                nb(gs, jg1[-g1_lo], d, m, 2),
                nb(gs, jg2[ g2_lo], d, m),
                nb(gs, jg2[-g2_lo], d, m, 2))
        },
        de = {
            cbind(
                nb(gs, jg1, d, m), 
                nb(gs, jg2, d, m, 2))
        },
        dp = {
            g1_lo <- sample(ng1, round(ng1 * 0.4))
            g2_lo <- sample(ng2, round(ng2 * 0.6))
            cbind(
                nb(gs, jg1[ g1_lo], d, m),
                nb(gs, jg1[-g1_lo], d, m, 2),
                nb(gs, jg2[ g2_lo], d, m),
                nb(gs, jg2[-g2_lo], d, m, 2))
        },
        dm = {
            g2_lo <- sample(ng2, round(ng2 * 0.6))
            cbind(
                nb(gs, jg1, d, m),
                nb(gs, jg2[ g2_lo], d, m),
                nb(gs, jg2[-g2_lo], d, m, 2))
        }, 
        db = {
            g2_lo <- sample(ng2, round(ng2 * 0.5))
            cbind(
                nb(gs, jg1, d, m, 2),
                nb(gs, jg2[ g2_lo], d, m),
                nb(gs, jg2[-g2_lo], d, m, 4))
        }
    )
}