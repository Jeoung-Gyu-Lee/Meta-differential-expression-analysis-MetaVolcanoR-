#' Perform Random-Effects Meta-Analysis and MetaVolcano Plot
#'
#' This function performs meta-analysis across multiple datasets using MetaVolcanoR (REM model)
#' and generates a MetaVolcano plot.
#'
#' @param gse_results A named list of GSE analysis results (e.g., output of per-dataset DEG analysis)
#' @param output_dir Directory where results and plots will be saved (default = meta_dir)
#' @param pval_cutoff Significance threshold for meta-analysis p-value (default = 0.05)
#' @param fc_cutoff Log2 fold change threshold for Up/Down classification (default = 1.0)
#' @param top_n Number of top DEGs to label in the volcano plot (default = 30)
#' @param hypothesis_min Minimum number of datasets a gene must appear in to be included (default = 2)
#'
#' @return A data.frame of meta-analysis results with regulation labels
#' @export
meta_rem_volcano <- function(gse_results,
                             output_dir = meta_dir,
                             pval_cutoff = 0.05,
                             fc_cutoff = 1.0,
                             top_n = 30,
                             hypothesis_min = 2) {
      
  # 1. Reuse the already loaded DEG result objects from gse_results
  limma_list <- lapply(gse_results, function(x) x$deg_limma)
  
  # Add a Gene column to each data.frame
  gene_list_df <- lapply(limma_list, function(df) {
    df$Gene <- rownames(df)
    df
  })
  
  # Extract genes per dataset
  gene_all <- lapply(gene_list_df, function(df) df$Gene)
  cat("gene counts of each dataset:", sapply(gene_all, length), "\n")
  
  # 2. Identify genes present in at least `hypothesis_min` datasets
  gene_count <- table(unlist(gene_all))
  gene_hypo <- names(gene_count)[gene_count >= hypothesis_min]
  
  cat(
    "Genes present in ≥ hypothesis_min datasets", hypothesis_min, ":", sum(gene_count >= hypothesis_min),
    "/ Overall gene number:", length(gene_count), "\n"
  )
  
  # 3. Filter by selected genes and clean entries
  gene_list_df <- lapply(gene_list_df, function(df) {
    df <- df[df$Gene %in%gene_hypo, ]
    df <- df %>% 
      dplyr::filter(!is.na(t), !is.na(logFC), t != 0, is.finite(logFC), is.finite(t))
    return(df)
  })
  
  # 4. Prepare list for REM
  diffexp_list <- lapply(gene_list_df, function(df) {
    df <- df %>%
      dplyr::mutate(
        SE = abs(logFC / t),
        CI_low = logFC - 1.96 * SE,
        CI_high = logFC + 1.96 * SE,
        Variance = SE^2
      ) %>%
      dplyr::filter(!is.na(SE), !is.na(CI_low), !is.na(CI_high), !is.na(Variance)) %>%
      dplyr::select(Gene, logFC, adj.P.Val, CI_low, CI_high, Variance) %>%
      dplyr::rename(Symbol = Gene, Log2FC = logFC, pvalcol = adj.P.Val)
  })
  names(diffexp_list) <- names(gene_list_df)
  
  # 5. Run MetaVolcano REM
  meta_degs_rem <- rem_mv(
    diffexp       = diffexp_list,
    pcriteria     = "pvalcol",
    foldchangecol = "Log2FC",
    genenamecol   = "Symbol",
    llcol         = "CI_low",
    rlcol         = "CI_high",
    vcol          = "Variance",
    cvar          = FALSE,
    metathr       = pval_cutoff,
    jobname       = "",
    outputfolder  = output_dir,
    draw          = "PDF",
    ncores        = 1
  )
  cat("Total genes with meta-analysis:", nrow(meta_degs_rem@metaresult), "\n")
  
  # 6. Extract and classify results
  meta_results <- as.data.frame(meta_degs_rem@metaresult) %>%
    mutate(
      regulation = dplyr::case_when(
        randomP < pval_cutoff & randomSummary >= fc_cutoff  ~ "Up",
        randomP < pval_cutoff & randomSummary <= -fc_cutoff ~ "Down",
        TRUE                                                   ~ "NS"
      )
    )
  
  print(table(meta_results$regulation))
  write.csv(gene_count, file.path(output_dir, "_Gene_All.csv"), row.names = FALSE) # For downstream analysis
  write.csv(meta_results, file.path(output_dir, "_MetaResults.csv"), row.names = FALSE)
  
  # 7. Volcano Plot using ggplot
  meta_df <- meta_results[!is.na(meta_results$randomP), ]
  
  # Grouping
  meta_df$group <- "NS"
  meta_df$group[meta_df$randomP < pval_cutoff & abs(meta_df$randomSummary) < fc_cutoff] <- "WeakDE"
  meta_df$group[meta_df$randomP < pval_cutoff & meta_df$randomSummary >= fc_cutoff] <- "Up"
  meta_df$group[meta_df$randomP < pval_cutoff & meta_df$randomSummary <= -fc_cutoff] <- "Down"
  
  # Label top N genes
  top_genes <- meta_df[meta_df$group %in% c("Up", "Down"), ]
  top_genes <- top_genes[order(top_genes$randomP), ][1:min(top_n, nrow(top_genes)), ]
  
  pdf(file.path(output_dir, "MetaVolcano.pdf"), width = 13, height = 11)
  ggplot(meta_df, aes(x = randomSummary, y = -log10(randomP))) +
    geom_errorbarh(data = subset(meta_df, group %in% c("Up", "Down")),
                   aes(xmin = randomCi.lb, xmax = randomCi.ub, color = group),
                   height = 0.1, alpha = 0.4) +
    geom_point(aes(color = group)) +
    geom_text_repel(data = top_genes, aes(label = Symbol), size = 4) +
    scale_color_manual(values = c("Up" = "tomato",
                                  "Down" = "skyblue",
                                  "WeakDE" = "gray70",
                                  "NS" = "gray85"),
                       limits = c("Up", "Down", "WeakDE", "NS")) +
    labs(title = "MetaVolcano",
         x = "Log2 Fold Change",
         y = "-log10(p-value)",
         color = "Group") +
    theme_minimal(base_size = 16) +
    theme(
      text = element_text(face = "bold"),
      panel.border = element_blank(),
      axis.line.x.bottom = element_line(color = "black", linewidth = 0.5),
      axis.line.y.left = element_line(color = "black", linewidth = 0.5),
      axis.line.x.top = element_blank(),
      axis.line.y.right = element_blank()
    )
  dev.off()
  return(meta_results)
}