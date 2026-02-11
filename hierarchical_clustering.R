#' Generate hierarchical clustering heatmap for top genes
#'
#' Combines normalized expression matrices across datasets, selects top N variable genes,
#' and generates a clustered heatmap with sample annotations.
#'
#' @param gse_results A named list of results from multiple GSE datasets (must contain exprs_final, pheno_data, gse_id)
#' @param top_n Number of genes to include in heatmap (default = 100)
#' @param output_dir Directory where the heatmap PDF will be saved (default = meta_dir)
#'
#' @return None (side effect: saves heatmap PDF)
#' @export
hierarchical_clustering <- function(gse_results,
                                     top_n = 100,
                                     output_dir = meta_dir) {
      
  # 1. Find intersecting genes across datasets
  common_genes <- Reduce(intersect, lapply(gse_results, function(x) rownames(x$exprs_final)))
  if (length(common_genes) == 0) {
    stop("No overlapping genes found across datasets.")
  }
  
  # 2. Combine expression matrices by common genes
  exprs_combined <- do.call(cbind, lapply(gse_results, function(x) x$exprs_final[common_genes, ]))
  
  # 3. Select top variable genes
  row_var <- apply(exprs_combined, 1, var)
  top_genes <- names(sort(row_var, decreasing = TRUE))[1:min(top_n, length(row_var))]
  exprs_top <- exprs_combined[top_genes, ]
  
  # 4. Build sample annotation dataframe
  sample <- do.call(rbind, lapply(gse_results, function(x) {
    data.frame(Sample = rownames(x$pheno_data),
               Group = x$pheno_data$group,
               Dataset = x$gse_id,
               stringsAsFactors = FALSE)
  }))
  rownames(sample) <- sample$Sample
  
  # 5. z-score row-wise + clipping
  exprs_z <- t(scale(t(exprs_top)))
  exprs_z[exprs_z > 3] <- 3
  exprs_z[exprs_z < -3] <- -3
  
  # 6. dataset color
  dataset_levels <- unique(sample$Dataset)
  dataset_colors <- setNames(RColorBrewer::brewer.pal(n = length(dataset_levels), "Set3"), dataset_levels)
  
  # 7. Create hierarchical_clustering
  output_file <- file.path(output_dir, paste0("hierarchical_clustering_Top", top_n, ".pdf"))
  pdf(output_file, width = 11, height = 10)
  pheatmap::pheatmap(
    exprs_z,
    color = colorRampPalette(c("green", "black", "red"))(100),
    annotation_col = sample[, c("Group", "Dataset")],
    annotation_colors = list(
      Group = c(Normal = "skyblue", Tumor = "tomato"),
      Dataset = dataset_colors
    ),
    treeheight_row = 120,
    treeheight_col = 40,
    show_colnames = FALSE,
    show_rownames = TRUE,
    clustering_distance_cols = "euclidean",
    clustering_distance_rows = "correlation",
    clustering_method = "ward.D2",
    fontsize = 6,
    main = paste("Hierarchical Clustering (Top", length(top_genes), "genes)")
  )
  dev.off()
}
