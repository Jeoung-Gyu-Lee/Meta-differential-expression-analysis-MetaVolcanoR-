#' Correct batch effect using ComBat if batch column is provided or auto-detected
#'
#' Applies ComBat to expression data if batch information is available.
#' Returns batch-corrected expression matrix or original matrix if no batch info is found.
#'
#' @param exprs_matrix Normalized expression matrix (genes x samples)
#' @param pheno_data Phenotype metadata (rows = samples, columns = annotations)
#' @param gse_id GEO Series ID (for plot title and output)
#' @param colour Column name used for PCA color grouping (default = "group")
#' @param batch_colname (Optional) Column name in `pheno_data` containing batch info.
#'                       If NULL, will auto-detect using 'batch' keyword.
#' @param output_dir Directory to save PCA plots (default = result_dir)
#'
#' @return A list containing:
#'   \item{exprs_corrected}{Batch-corrected expression matrix}
#'   \item{pheno_data}{Updated phenotype data with batch factor column}
#' @export
correct_batch <- function(exprs_matrix,
                          pheno_data,
                          gse_id,
                          colour = "group",
                          batch_colname = NULL,
                          output_dir = result_dir) {
  
  # Detect or validate batch column
  if (is.null(batch_colname)) {
    batch_colname <- grep("batch", colnames(pheno_data), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(batch_colname)) {
      message("Auto-detected batch column: ", batch_colname)
    }
  }
  
  # If batch column missing, skip ComBat and still do PCA plot
  if (is.null(batch_colname) || !(batch_colname %in% colnames(pheno_data))) {
    message("No valid batch column found — skipping ComBat. Returning original matrix.")
    return(list(exprs_corrected = exprs_matrix, pheno_data = pheno_data))
  }
  
  # Prepare batch variable
  pheno_data$batch <- factor(pheno_data[[batch_colname]])
  message("Batch distribution:")
  print(table(pheno_data$batch, useNA = "always"))
  
  # Design matrix for ComBat (based on group)
  mod_batch <- model.matrix(~ group, data = pheno_data)
  
  # Apply ComBat
  exprs_corrected <- sva::ComBat(
    dat         = exprs_matrix,
    batch       = pheno_data$batch,
    mod         = mod_batch,
    par.prior   = TRUE,
    prior.plots = FALSE
  )
  
  # Plot after batch correction
  plot_pca(exprs_corrected, pheno_data, paste0(gse_id, "_After_ComBat"),
           colour = colour, batch = "batch", output_dir = output_dir)
  
  return(list(exprs_corrected = exprs_corrected, pheno_data = pheno_data))
}
