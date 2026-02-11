#' Filter out low-expressed probes from expression matrix
#'
#' Removes probes (rows) whose log2 expression is below a given threshold in more than half the samples.
#' Typically used to eliminate background noise prior to downstream analysis.
#'
#' @param exprs_matrix A numeric expression matrix (rows = probes, columns = samples)
#' @param gse_id GEO dataset ID (for logging messages)
#' @param exprs_threshold Minimum log2 expression value to be considered expressed (default = 4)
#' @param min_ratio Minimum proportion of samples that must express the gene (default = 0.5)
#'
#' @return A filtered expression matrix with low-expression probes removed
#' @export
filter_low <- function(exprs_matrix,
                       gse_id,
                       exprs_threshold = 4,
                       min_ratio = 0.5) {
  probe_count <- nrow(exprs_matrix)
  if (!is.null(gse_id)) {
    cat(paste0("Initial probe count for ", gse_id, ": ", probe_count, "\n"))
  }
  
  min_samples_exprs <- ceiling(ncol(exprs_matrix) * min_ratio)
  
  keep_probes <- apply(exprs_matrix, 1, function(x) {
    sum(x >= exprs_threshold, na.rm = TRUE) >= min_samples_exprs
  })
  
  filtered_matrix <- exprs_matrix[keep_probes, ]
  
  if (!is.null(gse_id)) {
    removed <- probe_count - nrow(filtered_matrix)
    cat(paste0("Removed low expression probes for ", gse_id, ": ", removed, "\n"))
    cat(paste0("Remaining probes after filtering for ", gse_id, ": ", nrow(filtered_matrix), "\n"))
  }
  
  return(filtered_matrix)
}
