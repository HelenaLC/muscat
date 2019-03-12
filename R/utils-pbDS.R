# wrapper to create output tables
#   k:  cluster ID
#   tt: topTable data.frame
#   ct: comparison type; "contrast" or "coef"
#   c:  character string specifying the comparison
res_df <- function(k, tt, ct, c) {
    df <- data.frame(
        gene = rownames(tt), cluster_id = k, tt,
        row.names = NULL, stringsAsFactors = FALSE)
    df[[ct]] <- c
    return(df)
}
