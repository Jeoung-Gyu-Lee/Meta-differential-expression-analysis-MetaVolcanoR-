#' Generate a volcano plot using EnhancedVolcano
#'
#' Creates a volcano plot from a limma result table and saves it as a PDF.
#' Automatically handles 0 p-values and annotates DEG counts.
#'
#' @param deg_limma A data.frame from limma (must contain logFC and adj.P.Val columns)
#' @param gse_id GEO dataset ID (used for plot title and output filename)
#' @param output_dir Directory to save PDF file (default: result_dir)
#' @param p_cutoff Adjusted p-value threshold (default = 0.05)
#' @param fc_cutoff Log2 fold-change threshold (default = 1.0)
#'
#' @return None (side effect: saves PDF file)
#' @export
plot_volcano <- function(deg_limma,
                         gse_id,
                         output_dir = result_dir,
                         p_cutoff = 0.05,
                         fc_cutoff = 1.0) {
  
  # Filter valid entries
  volcano_df <- deg_limma %>% dplyr::filter(!is.na(logFC), !is.na(adj.P.Val))
  
  # Handle adj.P.Val == 0 (prevent plotting error)
  if (any(volcano_df$adj.P.Val == 0)) {
    min_p <- min(volcano_df$adj.P.Val[volcano_df$adj.P.Val > 0], na.rm = TRUE)
    if (is.infinite(min_p)) min_p <- 1e-300
    volcano_df$adj.P.Val[volcano_df$adj.P.Val == 0] <- min_p / 10
    cat("Adjusted 0 adj.P.Val values for plotting.\n")
  }
  
  # Define color groups
  keyvals <- rep("grey75", nrow(volcano_df))
  names(keyvals) <- rep("NS", nrow(volcano_df))
  
  up <- volcano_df$logFC > fc_cutoff & volcano_df$adj.P.Val < p_cutoff
  down <- volcano_df$logFC < -fc_cutoff & volcano_df$adj.P.Val < p_cutoff
  fc_only <- abs(volcano_df$logFC) > fc_cutoff & volcano_df$adj.P.Val >= p_cutoff
  
  keyvals[up] <- "tomato"
  names(keyvals)[up] <- "Up"
  keyvals[down] <- "skyblue"
  names(keyvals)[down] <- "Down"
  keyvals[fc_only] <- "lightgreen"
  names(keyvals)[fc_only] <- "Log2FC only"
  
  # Count DEGs
  up_count <- sum(up)
  down_count <- sum(down)
  
  # Output file path
  pdf_path <- file.path(output_dir, paste0(gse_id, "_VOLCANO.pdf"))
  pdf(pdf_path, width = 10, height = 7)
  
  # Generate plot
  volcano_plot <- EnhancedVolcano(volcano_df,
                                  x = "logFC",
                                  y = "adj.P.Val",
                                  pCutoff = p_cutoff,
                                  FCcutoff = fc_cutoff,
                                  colCustom = keyvals,
                                  title = paste0("Volcano Plot (", gse_id, ")"),
                                  xlab = bquote(~Log[2]~ 'fold change'),
                                  ylab = bquote(~-Log[10] ~ adjusted~italic(P)),
                                  legendLabels = c("NS", "Log2FC only", "Down", "Up"),
                                  colAlpha = 0.8,
                                  pointSize = 2.5,
                                  lab = NA)
  
  print(volcano_plot)
  
  # Annotate DEG count
  grid.text(paste0("Up: ", up_count),
            x = unit(0.80, "npc"), y = unit(0.75, "npc"),
            gp = gpar(col = "black", fontsize = 12))
  
  grid.text(paste0("Down: ", down_count),
            x = unit(0.30, "npc"), y = unit(0.75, "npc"),
            gp = gpar(col = "black", fontsize = 12))
  
  dev.off()
}
