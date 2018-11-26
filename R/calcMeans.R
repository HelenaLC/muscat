#' Calculate cluster medians
#' 
#' Calculate cluster medians (median expression for each cluster-sample-marker
#' combination)
#' 
#' Calculate median marker expression for each cluster and sample (i.e. medians for each 
#' cluster-sample-marker combination).
#' 
#' The data object is assumed to contain a factor \code{marker_class} in the column
#' meta-data (see \code{\link{prepareData}}), which indicates the protein marker class for
#' each column of data (\code{"type"}, \code{"state"}, or \code{"none"}).
#' 
#' The cluster medians are required for testing for differential states within cell
#' populations, and for plotting purposes.
#' 
#' Variables \code{id_type_markers} and \code{id_state_markers} are saved in the
#' \code{metadata} slot of the output object. These can be used to identify the 'cell
#' type' and 'cell state' markers in the list of \code{assays} in the output
#' \code{\link{SummarizedExperiment}} object, which is useful in later steps of the
#' 'diffcyt' pipeline.
#' 
#' Results are returned as a new \code{\link{SummarizedExperiment}} object, where rows =
#' clusters, columns = samples, sheets (\code{assays} slot) = markers. Note that there is
#' a separate table of values (\code{assay}) for each marker. The \code{metadata} slot
#' also contains variables \code{id_type_markers} and \code{id_state_markers}, which can
#' be used to identify the sets of cell type and cell state markers in the list of
#' \code{assays}.
#' 
#' 
#' @param d_se Data object from previous steps, in \code{\link{SummarizedExperiment}}
#'   format, containing cluster labels as a column in the row meta-data (from
#'   \code{\link{generateClusters}}). Column meta-data is assumed to contain a factor
#'   \code{marker_class}.
#' 
#' 
#' @return \code{d_medians}: \code{\link{SummarizedExperiment}} object, where rows =
#'   clusters, columns = samples, sheets (\code{assays} slot) = markers. The
#'   \code{metadata} slot contains variables \code{id_type_markers} and
#'   \code{id_state_markers}, which can be accessed with
#'   \code{metadata(d_medians)$id_type_markers} and
#'   \code{metadata(d_medians)$id_state_markers}.
#' 
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment assay rowData colData
#' @importFrom dplyr group_by tally summarize
#' @importFrom tidyr complete
#' @importFrom reshape2 acast
#' @importFrom magrittr %>%
#' @importFrom methods is
#' 
#' @export
#' 
#' @examples
#' # For a complete workflow example demonstrating each step in the 'diffcyt' pipeline, 
#' # see the package vignette.
#' 
#' # Function to create random data (one sample)
#' d_random <- function(n = 20000, mean = 0, sd = 1, ncol = 20, cofactor = 5) {
#'   d <- sinh(matrix(rnorm(n, mean, sd), ncol = ncol)) * cofactor
#'   colnames(d) <- paste0("marker", sprintf("%02d", 1:ncol))
#'   d
#' }
#' 
#' # Create random data (without differential signal)
#' set.seed(123)
#' d_input <- list(
#'   sample1 = d_random(), 
#'   sample2 = d_random(), 
#'   sample3 = d_random(), 
#'   sample4 = d_random()
#' )
#' 
#' experiment_info <- data.frame(
#'   sample_id = factor(paste0("sample", 1:4)), 
#'   group_id = factor(c("group1", "group1", "group2", "group2")), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' marker_info <- data.frame(
#'   channel_name = paste0("channel", sprintf("%03d", 1:20)), 
#'   marker_name = paste0("marker", sprintf("%02d", 1:20)), 
#'   marker_class = factor(c(rep("type", 10), rep("state", 10)), 
#'                         levels = c("type", "state", "none")), 
#'   stringsAsFactors = FALSE
#' )
#' 
#' # Prepare data
#' d_se <- prepareData(d_input, experiment_info, marker_info)
#' 
#' # Transform data
#' d_se <- transformData(d_se)
#' 
#' # Generate clusters
#' d_se <- generateClusters(d_se)
#' 
#' # Calculate medians
#' d_medians <- calcMedians(d_se)
#' 
calcMeans <- function(d_se) {
  
  if (!is(d_se, "SummarizedExperiment")) {
    stop("Data object must be a 'SummarizedExperiment'")
  }
  
  if (!("cluster_id" %in% (colnames(rowData(d_se))))) {
    stop("Data object does not contain cluster labels. Run 'generateClusters' to generate cluster labels.")
  }
  
  is_marker <- colData(d_se)$marker_class != "none"
  
  # identify 'cell type' and 'cell state' markers in final list of assays
  id_type_markers <- (colData(d_se)$marker_class == "type")[is_marker]
  id_state_markers <- (colData(d_se)$marker_class == "state")[is_marker]
  
  # calculate cluster medians for each marker
  assaydata_mx <- assay(d_se)
  
  medians <- vector("list", sum(is_marker))
  marker_names_sub <- as.character(colData(d_se)$marker_name[is_marker])
  names(medians) <- marker_names_sub
  
  clus <- rowData(d_se)$cluster_id
  smp <- rowData(d_se)$sample_id
  
  for (i in seq_along(medians)) {
    assaydata_i <- assaydata_mx[, marker_names_sub[i], drop = FALSE]
    assaydata_i <- as.data.frame(assaydata_i)
    assaydata_i <- cbind(assaydata_i, sample_id = smp, cluster_id = clus)
    colnames(assaydata_i)[1] <- "value"
    
    assaydata_i %>% 
      group_by(cluster_id, sample_id) %>% 
      summarize(mean = mean(value)) -> 
      med
    
    med <- acast(med, cluster_id ~ sample_id, value.var = "mean", fill = NA)
    
    # fill in any missing clusters
    if (nrow(med) < nlevels(rowData(d_se)$cluster_id)) {
      ix_missing <- which(!(levels(rowData(d_se)$cluster_id) %in% rownames(med)))
      med_tmp <- matrix(NA, nrow = length(ix_missing), ncol = ncol(med))
      rownames(med_tmp) <- ix_missing
      med <- rbind(med, med_tmp)
      # re-order rows
      med <- med[order(as.numeric(rownames(med))), , drop = FALSE]
    }
    
    medians[[i]] <- med
  }
  
  # check cluster IDs and sample IDs are identical
  for (i in seq_along(medians)) {
    if (!all(rownames(medians[[i]]) == rownames(medians[[1]]))) {
      stop("Cluster IDs do not match")
    }
    if (!all(colnames(medians[[i]]) == colnames(medians[[1]]))) {
      stop("Sample IDs do not match")
    }
  }
  
  # create new SummarizedExperiment (rows = clusters, columns = samples)
  
  row_data <- data.frame(
    cluster_id = factor(rownames(medians[[1]]), levels = levels(rowData(d_se)$cluster_id)), 
    stringsAsFactors = FALSE
  )
  
  col_data <- metadata(d_se)$experiment_info
  
  # rearrange sample order to match 'experiment_info'
  medians <- lapply(medians, function(m) {
    m[, match(col_data$sample_id, colnames(m)), drop = FALSE]
  })
  stopifnot(all(sapply(medians, function(m) {
    col_data$sample_id == colnames(m)
  })))
  
  metadata <- list(id_type_markers = id_type_markers, 
                   id_state_markers = id_state_markers,
                   clustering_name = metadata(d_se)$clustering_name)
  
  d_medians <- SummarizedExperiment(
    assays = list(medians = medians), 
    rowData = row_data, 
    colData = col_data, 
    metadata = metadata
  )
  
  d_medians
}


