# main.R - DEG analysis and meta-analysis
# Load Functions
gc()
source("functions/install_packages.R")
source("functions/load_packages.R")
source("functions/download_data.R")
source("functions/define_group.R")
source("functions/normalize_platform.R")
source("functions/plot_box.R")
source("functions/plot_density.R")
source("functions/plot_pca.R")
source("functions/correct_batch.R")
source("functions/filter_low.R")
source("functions/gene_mapping.R")
source("functions/collapse_probes.R")
source("functions/filter_coding.R")
source("functions/run_limma_DEG.R")
source("functions/plot_volcano.R")
source("functions/hierarchical_clustering.R")
source("functions/meta_rem_volcano.R")

# Install and load required packages
# Custorm CDF and 'MetaVolcanoR' source packages are must exist in 'project' directory (eg. hgu133plus2hsentrezg.db..)
install_packages()
load_packages()

# Define GSE IDs
gse_ids <- c(
  "GSE19804", "GSE18842", "GSE19188", "GSE101929", 
  "GSE27262", "GSE118370", "GSE130779", "GSE136043",
  "GSE33532", "GSE138682"
)

# Container for all results
gse_results <- list()

# Log start
sink("analysis_log.txt", split = TRUE)
error_con <- file("analysis_error.txt", open = "wt")
sink(error_con, type = "message")

# Loop through each dataset
for (gse_id in gse_ids) {
  cat("\n========== Processing", gse_id, "==========\n")
  tryCatch({
    
  # 1. Download and Load
  data_list <- download_data(gse_id)
  pheno_data <- data_list$pheno_data
  raw_files <- data_list$raw_files
  gse_dir <- data_list$gse_dir
  result_dir <- data_list$result_dir
  
  # 2. Assign Group
  pheno_data$group <- define_group(gse_id, pheno_data)
  print(table(pheno_data$group, useNA = "ifany"))
  
  # 3. Normalize Platform
  norm_list <- normalize_platform(gse_id, pheno_data, raw_files)
  exprs_before <- norm_list$exprs_before
  exprs_after <- norm_list$exprs_after
  
  # 4. Quality Control Plots
  plot_box(exprs_before, gse_id, stage = "Before", output_format = "png", output_dir = result_dir)
  plot_box(exprs_after, gse_id, stage = "After", output_format = "pdf", output_dir = result_dir)
  plot_density(exprs_before, gse_id, stage = "Before", output_dir = result_dir)
  plot_density(exprs_after, gse_id, stage = "After", output_dir = result_dir)
  plot_pca(exprs_after, pheno_data, gse_id, output_dir = result_dir)
  
  # 5. Batch Correction
  res <- correct_batch(exprs_after, pheno_data, gse_id)
  exprs_batch <- res$exprs_corrected
  pheno_data <- res$pheno_data
  
  # 6. Filter Low Expression
  exprs_filter <- filter_low(exprs_batch, gse_id)
  
  # 7. Probe Gene Mapping
  mapping_list <- gene_mapping(exprs_filter, pheno_data, gse_id)
  mapping <- mapping_list$mapping
  exprs_mapped <- mapping_list$exprs_mapped
  
  # 8. Collapse_probes
  exprs_final <- collapse_probes(exprs_mapped, gse_id)
  
  # 9. Filter Protein Coding
  exprs_final <- filter_coding(exprs_final, gse_id)
  
  # 10. DEG Analysis
  deg_limma <- run_limma_DEG(exprs_final, pheno_data, gse_id, output_dir = result_dir, 
                             pval_cutoff = 0.05, logfc_cutoff = 1.0)
  
  # 11. Volcano Plot
  plot_volcano(deg_limma, gse_id, output_dir = result_dir, p_cutoff = 0.05, fc_cutoff = 1.0)
  
  # 12. Save result in list
  gse_results[[gse_id]] <- list(
    pheno_data = pheno_data,
    exprs_final = exprs_final,
    deg_limma = deg_limma,
    gse_id = gse_id
  )
  }, error = function(e) {
    cat("Error in", gse_id, ":", conditionMessage(e), "\n")
  })
}

# Log end
sink(type = "message")
close(error_con)        
sink() # Close sink  

# Export combined DEG lists (from each dataset, deg_limma)
combined_deg <- do.call(rbind, lapply(names(gse_results), function(gse_id) {
  df <- gse_results[[gse_id]]$deg_limma
  df$Gene    <- rownames(df)
  df$Dataset <- gse_id
  df[df$Regulation %in% c("Up", "Down"), ]
})) %>%
  dplyr::select(Regulation, Gene, Dataset)
rownames(combined_deg) <- NULL

meta_dir <- file.path(getwd(), "meta_results")
write.csv(combined_deg, file = file.path(meta_dir, "Combined_DEG.csv"), row.names = FALSE)

# 13. Combined Density plot with Multiple dataset, after normalizationn
plot_density(stage = "After", all_results = gse_results, gse_id = "All_Combined", output_dir = meta_dir)

# 14. Hierarchical Clustering Heatmap
hierarchical_clustering(gse_results, output_dir = meta_dir)

# 15. Meta Analysis for MetaVolcanoR
meta_result <- meta_rem_volcano(
  gse_results = gse_results,
  output_dir = meta_dir,
  pval_cutoff = 0.05,
  fc_cutoff = 0.5,
  top_n = 30,
  hypothesis_min = 7
)
