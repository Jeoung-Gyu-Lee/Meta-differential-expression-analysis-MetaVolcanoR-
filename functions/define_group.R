#' Assign group labels ("Tumor" vs "Normal") to samples
#'
#' This function applies a dataset-specific rule to assign group labels based on phenotype metadata.
#'
#' @param gse_id GEO Series ID
#' @param pheno_data A data.frame containing phenotype metadata (from pData())
#'
#' @return A factor vector indicating group assignment ("Normal", "Tumor")
#' @export
define_group <- function(gse_id, pheno_data) {
  group <- rep(NA, nrow(pheno_data))
  
  # Dataset-specific group classification rules
  rule_list <- list(
    GSE19804   = list(column = "title", tumor = "cancer", normal = "normal"),
    GSE18842   = list(column = "title", tumor = "Tumor", normal = "control"),
    GSE19188   = list(column = "tissue type:ch1", tumor = "Tumor", normal = "healthy"),
    GSE101929  = list(column = "tumor_normal status:ch1", tumor = "^T$", normal = "^N$"),
    GSE27262   = list(column = "title", tumor = "[0-9]+T$", normal = "[0-9]+N$"),
    GSE118370  = list(column = "title", tumor = "adenocarcinoma", normal = "normal"),
    GSE33532   = list(column = "title", tumor = "tumor", normal = "normal"),
    GSE138682  = list(column = "title", tumor = "Tumor", normal = "control"),
    GSE75037   = list(column = "title", tumor = "+T$", normal = "+N$"),
    GSE32863   = list(column = "title", tumor = "adenocarcinoma", normal = "non-tumor"),
    GSE130779  = list(column = "tissue:ch1", tumor = "adenocarcinoma", normal = "non-tumor"),
    GSE136043  = list(column = "tissue:ch1", tumor = "adenocarcinoma", normal = "non-tumor")
  )
  
  rule <- rule_list[[gse_id]]
  
  if (!is.null(rule)) {
    col <- rule$column
    if (!col %in% colnames(pheno_data)) {
      stop(paste("Column", col, "not found in pheno_data for", gse_id))
    }
    
    group[grepl(rule$tumor, pheno_data[[col]], ignore.case = TRUE)] <- "Tumor"
    group[grepl(rule$normal, pheno_data[[col]], ignore.case = TRUE)] <- "Normal"
    
  } else {
    stop(paste("No group rule defined for", gse_id))
  }
  return(factor(group, levels = c("Normal", "Tumor")))
}
