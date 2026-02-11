# ---------------------------------------------------------------------------- #
# 1. Install and Load Packages
# ---------------------------------------------------------------------------- #
# Installs necessary CRAN, Bioconductor, and local .tar.gz packages if not already installed.
install_packages <- function() {
 
# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Local Custom Packages (.tar.gz) as platform type
install.packages("hgu133plus2hsentrezg.db_25.0.0.tar.gz", repos = NULL, type = "source") # Brainarray Custom CDF (hgu133plus2)
install.packages("pd.hgu133plus2.hs.entrezg_25.0.0.tar.gz", repos = NULL, type = "source") # Platform design (hgu133plus2)
install.packages("primeviewhsentrezg.db_25.0.0.tar.gz", repos = NULL, type = "source")  # Brainarray Custom CDF (primeview)
install.packages("pd.primeview.hs.entrezg_25.0.0.tar.gz", repos = NULL, type = "source")  # Platform design (primeview)
install.packages("hta20hsentrezg.db_25.0.0.tar.gz", repos = NULL, type = "source")  # Brainarray Custom CDF (Affymetrix Human Transcriptome Array 2.0)
install.packages("pd.hta20.hs.entrezg_25.0.0.tar.gz", repos = NULL, type = "source")  # Platform design (Affymetrix Human Transcriptome Array 2.0)
BiocManager::install("illuminaHumanv3.db", dependencies = TRUE, force = TRUE) # Annotation database (illuminaHumanv3)

# CRAN Packages
install.packages("R.utils", dependencies = TRUE)
install.packages("stringr", dependencies = TRUE)
install.packages("openxlsx", dependencies = TRUE)
install.packages("dplyr", dependencies = TRUE)
install.packages("pheatmap", dependencies = TRUE)
install.packages("ggplot2", dependencies = TRUE)
install.packages("tidyr", dependencies = TRUE)
install.packages("ggfortify", dependencies = TRUE)
install.packages("tibble", dependencies = TRUE)
install.packages("plotly", dependencies = TRUE)
install.packages("WGCNA", dependencies = TRUE)
install.packages("matrixStats", dependencies = TRUE)
install.packages("reshape2", dependencies = TRUE)
install.packages("ggrepel", dependencies = TRUE)
install.packages("RColorBrewer", dependencies = TRUE)
install.packages("metap", dependencies = TRUE) # Dependency packages for "MetaVolcanoR"
install.packages("topconfects_1.25.0.tar.gz", repos = NULL, type = "source") # Dependency packages for "MetaVolcanoR"
install.packages("MetaVolcanoR_1.10.0.tar.gz", repos = NULL, type = "source")

# Bioconductor Packages
BiocManager::install("GEOquery", dependencies = TRUE)
BiocManager::install("org.Hs.eg.db", dependencies = TRUE)
BiocManager::install("AnnotationDbi", dependencies = TRUE, force = TRUE)
BiocManager::install("limma", dependencies = TRUE)
BiocManager::install("EnhancedVolcano", dependencies = TRUE)
BiocManager::install("oligo", dependencies = TRUE, force = TRUE)
BiocManager::install("sva", dependencies = TRUE)
BiocManager::install("biomaRt", dependencies = TRUE)
BiocManager::install("impute")  # for install 'WGCNA'
}