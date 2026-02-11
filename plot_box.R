#' Generate boxplot for expression matrix (before/after normalization)
#'
#' This function generates a sample-wise boxplot of log2 expression values using ggplot2
#' and saves it to a PNG or PDF file.
#'
#' @param exprs_matrix A numeric matrix of expression values (rows = probes, columns = samples)
#' @param gse_id GEO dataset ID used for naming the output file
#' @param stage A string indicating stage of normalization ("Before" or "After")
#' @param output_format Either "png" or "pdf"
#' @param output_dir Directory to save the plot file (default = result_dir)
#'
#' @return A long-format data.frame used for plotting (columns: ProbeID, Sample, Expression)
#' @export
plot_box <- function(exprs_matrix,
                     gse_id,
                     stage = c("Before", "After"),
                     output_format = c("png", "pdf"),
                     output_dir = result_dir) {
  
  stage <- match.arg(stage)
  output_format <- match.arg(output_format)
  
  # Reshape to long format
  exprs_long <- as.data.frame(exprs_matrix) %>%
    tibble::rownames_to_column("ProbeID") %>%
    tidyr::pivot_longer(cols = -ProbeID, names_to = "Sample", values_to = "Expression")
  
  # Define output file name
  file_name <- file.path(output_dir, paste0(gse_id, "_BoxPlot_", stage, ".", output_format))
  
  # Open graphics device
  if (output_format == "png") {
    png(file_name, width = 800, height = 600, res = 100)
  } else if (output_format == "pdf") {
    pdf(file_name, width = 10, height = 7)
  }
  
  # Generate plot
  plot_obj <- ggplot(exprs_long, aes(x = Sample, y = Expression)) +
    geom_boxplot() +
    labs(
      title = paste0(gse_id, " - Log2 Expression Distribution (", stage, " Normalization)"),
      x = "Sample ID", y = "Log2 Expression"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5))
  
  print(plot_obj)
  dev.off()
  
  return(exprs_long)
}
