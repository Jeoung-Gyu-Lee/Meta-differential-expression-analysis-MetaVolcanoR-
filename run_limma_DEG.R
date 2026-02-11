#' Perform differential expression analysis using limma
#'
#' Fits a linear model to the expression matrix using the group column in pheno_data,
#' applies empirical Bayes moderation, and returns the result.
#'
#' @param exprs_matrix A numeric matrix of expression data (genes x samples)
#' @param pheno_data A data.frame containing sample metadata with a 'group' column (e.g., Tumor vs Normal)
#' @param gse_id GEO dataset ID (for logging and file naming)
#' @param output_dir Directory to save output CSV files (default = result_dir)
#' @param pval_cutoff Adjusted p-value cutoff for DEG selection (default = 0.05)
#' @param logfc_cutoff Absolute log2 fold-change cutoff for DEG selection (default = 1.0)
#'
#' @return A data.frame with limma differential expression results and Regulation label
#' @export
run_limma_DEG <- function(exprs_matrix,
                          pheno_data,
                          gse_id,
                          output_dir = result_dir,
                          pval_cutoff = 0.05,
                          logfc_cutoff = 1.0) {
  
  # Design matrix
  design <- with(pheno_data, model.matrix(~ 0 + group))
  colnames(design) <- make.names(colnames(design))
  colnames(design) <- gsub("group", "", colnames(design))  # e.g., groupTumor → Tumor
  
  # Check rank
  if (qr(design)$rank < ncol(design)) {
    stop(paste0("Design matrix is not full rank for ", gse_id, ". Check sample group sizes."))
  }
  
  # Fit linear model and contrast
  fit <- limma::lmFit(exprs_matrix, design)
  contrast.matrix <- limma::makeContrasts(
    contrasts = paste0(levels(pheno_data$group)[2], "-", levels(pheno_data$group)[1]),
    levels = design
  )
  fit2 <- limma::contrasts.fit(fit, contrast.matrix)
  fit2 <- limma::eBayes(fit2)
  
  # Extract results
  deg_limma <- limma::topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "P")
  cat(paste0("Dimension of topTable for ", gse_id, ": ", nrow(deg_limma), " x ", ncol(deg_limma), "\n"))
  
  # Add regulation column using cut off
  deg_limma$Regulation <- ifelse(deg_limma$adj.P.Val < pval_cutoff & deg_limma$logFC > logfc_cutoff, "Up",
                                 ifelse(deg_limma$adj.P.Val < pval_cutoff & deg_limma$logFC < -logfc_cutoff, "Down", "NS"))
  print(table(deg_limma$Regulation))
  
  # Save metadata and DEG results
  write.csv(pheno_data, file.path(output_dir, paste0(gse_id, "_Metadata.csv")), row.names = TRUE)
  write.csv(deg_limma, file.path(output_dir, paste0(gse_id, "_Results.csv")), row.names = TRUE)
  
  return(deg_limma)
}
