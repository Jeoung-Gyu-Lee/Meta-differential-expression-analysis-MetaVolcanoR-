#' Generate density plot for expression matrix (per-sample or combined)
#'
#' This function creates:
#' - individual sample-wise density plots (for a single dataset)
#' - combined dataset density plot (after normalization)
#'
#' @param exprs_matrix Numeric expression matrix (rows = probes, columns = samples); used for individual plot
#' @param gse_id GEO Series ID (used in plot title and output file name)
#' @param stage Normalization stage ("Before" or "After")
#' @param all_results Optional list of all GSE results to draw combined density (requires `exprs_nor_after`)
#' @param combined_density Logical; if TRUE and stage == "After", draw combined plot across datasets
#' @param output_dir Directory to save plots (default = result_dir)
#'
#' @return No return value. Saves plot files to disk.
#' @export
plot_density <- function(exprs_matrix = NULL,
                         gse_id = NULL,
                         stage = c("Before", "After"),
                         all_results = NULL,
                         combined_density = TRUE,
                         output_dir = result_dir) {
  
  stage <- match.arg(stage)
  
  # 1. Individual Density Plot
  if (!is.null(exprs_matrix) && !is.null(gse_id)) {
    exprs_long <- as.data.frame(exprs_matrix)
    
    if (!"ProbeID" %in% colnames(exprs_long)) {
      exprs_long <- tibble::rownames_to_column(exprs_long, "ProbeID")
    }
    
    exprs_long <- exprs_long %>%
      tidyr::pivot_longer(cols = -ProbeID, names_to = "Sample", values_to = "Expression")
    
    output_file <- file.path(output_dir, paste0(gse_id, "_DensityPlot_", stage, ".pdf"))
    pdf(output_file, width = 10, height = 7)
    
    p <- ggplot(exprs_long, aes(x = Expression, color = Sample)) +
      geom_density() +
      labs(
        title = paste0(gse_id, " - Log2 Expression Density Distribution (", stage, " Normalization)"),
        x = "Log2 Expression", y = "Density"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)
      )
    
    print(p)
    dev.off()
  }
 
  # 2. Combined Density Plot (Across All Datasets)
  if (!is.null(all_results) && stage == "After" && combined_density) {
    density_all <- do.call(rbind, lapply(all_results, function(x) {
      df <- as.data.frame(x$exprs_final) %>%
        tibble::rownames_to_column("ProbeID") %>%
        tidyr::pivot_longer(cols = -ProbeID, names_to = "Sample", values_to = "Expression")
      df$Dataset <- x$gse_id
      return(df)
    }))
    
    output_file <- file.path(output_dir, paste0(gse_id, "_DensityPlot_Combined.pdf"))
    pdf(output_file, width = 10, height = 7)
    
    p <- ggplot(density_all, aes(x = Expression, color = Dataset, fill = Dataset)) +
      geom_density(alpha = 0.5) +
      labs(
        title = "Expression Density Across Datasets (After Normalization)",
        x = "Log2 Expression", y = "Density"
      ) +
      theme_minimal()
    
    print(p)
    dev.off()
  }
}
