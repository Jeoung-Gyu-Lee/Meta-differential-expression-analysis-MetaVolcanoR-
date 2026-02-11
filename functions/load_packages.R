#' Load all required packages for DEG and meta-analysis pipeline
#'
#' This function loads all necessary CRAN and Bioconductor packages
#' for microarray normalization, visualization, annotation, and meta-analysis.
#'
#' @export
load_packages <- function() {
  pkgs <- c(
    "pd.hgu133plus2.hs.entrezg",
    "hgu133plus2hsentrezg.db",
    "primeviewhsentrezg.db",
    "pd.primeview.hs.entrezg",
    "hta20hsentrezg.db",
    "pd.hta20.hs.entrezg",
    "illuminaHumanv3.db",
    "limma",
    "R.utils",
    "GEOquery",
    "openxlsx",
    "dplyr",
    "stringr",
    "pheatmap",
    "EnhancedVolcano",
    "ggplot2",
    "tidyr",
    "ggfortify",
    "oligo",
    "tibble",
    "plotly",
    "AnnotationDbi",
    "WGCNA",
    "matrixStats",
    "grid",
    "sva",
    "biomaRt",
    "ggrepel",
    "RColorBrewer",
    "MetaVolcanoR"
  )
  
  sapply(pkgs, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  })
  
  cat("✅ All packages loaded successfully.\n")
}
