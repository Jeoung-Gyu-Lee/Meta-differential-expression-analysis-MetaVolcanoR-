#' Perform PCA and generate a 2D PCA plot
#'
#' Computes PCA from an expression matrix (samples in columns) and visualizes the first two principal components.
#' Optionally includes sample groupings by color and batch effects by shape.
#'
#' @param exprs_matrix A numeric expression matrix (rows = genes, columns = samples)
#' @param pheno_data A phenotype data.frame with sample metadata (row names must match sample names)
#' @param gse_id GEO dataset ID (used in plot title and filename)
#' @param colour Column name in `pheno_data` to use for color grouping (default: "group")
#' @param batch Optional column name in `pheno_data` to use for shape (batch effect); set to NULL to disable
#' @param output_dir Directory to save the PCA plot (default = result_dir)
#'
#' @return A list containing:
#'   \item{pca_result}{Output from `prcomp`}
#'   \item{pca_data}{PCA data.frame used for plotting}
#' @export
plot_pca <- function(exprs_matrix,
                     pheno_data,
                     gse_id,
                     colour = "group",
                     batch = NULL,
                     output_dir = result_dir) {
  
  # Compute PCA (samples in rows)
  pca_result <- prcomp(t(exprs_matrix), center = TRUE, scale. = FALSE)
  
  # Prepare PCA data frame
  pca_df <- as.data.frame(pca_result$x)
  pca_df$Sample <- rownames(pca_df)
  pca_df[[colour]] <- pheno_data[[colour]][match(pca_df$Sample, rownames(pheno_data))]
  
  batch_exists <- !is.null(batch) && batch %in% colnames(pheno_data)
  if (batch_exists) {
    pca_df[[batch]] <- pheno_data[[batch]][match(pca_df$Sample, rownames(pheno_data))]
  }
  
  # Variance explained
  var_explained <- summary(pca_result)$importance[2, 1:2] * 100
  subtitle_text <- paste0("Variance Explained: PC1 = ",
                          round(var_explained[1], 2), "%, PC2 = ",
                          round(var_explained[2], 2), "%")
  
  # File output path
  output_file <- file.path(output_dir, paste0(gse_id, "_PCA.pdf"))
  pdf(output_file, width = 10, height = 7)
  
  # Generate PCA plot
  pca_plot <- autoplot(pca_result,
                        data = pca_df,
                        colour = colour,
                        size = 3,
                        frame = TRUE,
                        frame.type = "norm")
  
  if (batch_exists) {
    pca_plot <- pca_plot + aes_string(shape = batch)
  }
  
  pca_plot <- pca_plot +
    labs(title = paste0(gse_id, " - PCA Plot",
                        if (batch_exists) " (With Batch Shape)" else ""),
         subtitle = subtitle_text) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  print(pca_plot)
  dev.off()
  
  invisible(list(pca_result = pca_result, pca_data = pca_df))
}
