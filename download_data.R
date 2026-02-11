#' Download and prepare GEO expression data and metadata
#' Downloads raw files (CEL/txt), decompresses them, and loads phenotype metadata using GEOquery.
#'
#' @param gse_id GEO Series ID (e.g., "GSE136043")
#'
#' @return A list containing:
#'   \item{pheno_data}{Phenotype metadata (data.frame)}
#'   \item{raw_files}{Character vector of decompressed raw file paths}
#'   \item{gse_dir}{The directory where files were saved}
#'   \item{result_dir}{The directory where results will be saved}
#' @export
download_data <- function(gse_id) {
  
  # Set data and result directories
  base_dir   <- getwd()
  gse_dir    <- file.path(base_dir, "data", gse_id)
  result_dir <- file.path(base_dir, "results", gse_id)
  
  if (!dir.exists(gse_dir)) {
    dir.create(gse_dir, recursive = TRUE)
  }
  if (!dir.exists(result_dir)) {
    dir.create(result_dir, recursive = TRUE)
  }
  
  # Save current working directory and change to gse_dir
  old_wd <- base_dir
  setwd(gse_dir)
  
  # Download raw data in compressed format (.tar, .gz, etc.)
  com_files <- list.files(gse_dir, pattern = "\\.(tar|tar\\.gz|tgz|gz)$", full.names = TRUE, recursive = TRUE)
  
  if (length(com_files) == 0) {
    options(timeout = 1000)
    getGEOSuppFiles(gse_id, baseDir = gse_dir)
    cat(paste0("Downloaded raw data for ", gse_id, "\n"))
  } else {
    cat(paste0("Raw data for ", gse_id, " already exist. Skipping download.\n"))
  }
  
  # Decompress files
  comp_files <- list.files(path = gse_dir, pattern = "\\.(tar|tar\\.gz|tgz|gz)$",
                           full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
  
  if (length(comp_files) > 0) {
    for (file in comp_files) {
      if (grepl("\\.(tar|tar\\.gz|tgz)$", file, ignore.case = TRUE)) {
        untar(file, exdir = gse_id)
      } else if (grepl("\\.gz$", file, ignore.case = TRUE)) {
        gunzip(file, remove = TRUE)
      }
    }
    cat(paste0("Decompressed files for ", gse_id, "\n"))
  } else {
    cat("No .tar/.gz files found — skipping.\n")
  }
  
  # Decompress any remaining .CEL.gz or .txt.gz
  gz_files <- list.files(path = gse_dir, pattern = "\\.CEL\\.gz$|\\.txt\\.gz$",
                         full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
  
  if (length(gz_files) > 0) {
    for (file in gz_files) {
      gunzip(file, remove = TRUE)
    }
    cat(paste0("Decompressed ", length(gz_files), " .gz files (CEL/txt)\n"))
  } else {
    cat("All CEL/txt files already decompressed or not found.\n")
  }
  
  # List of raw files (.CEL or .txt)
  raw_files <- list.files(path = gse_dir, pattern = "\\.CEL$|\\.txt$",
                          full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
  
  # Load metadata using GEOquery
  gse <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = FALSE)
  
  if (inherits(gse, "ExpressionSet")) {
    pheno_data <- pData(gse)
  } else if (is.list(gse) && all(sapply(gse, inherits, "ExpressionSet"))) {
    pheno_data <- do.call(rbind, lapply(gse, pData))
  } else {
    stop(paste0("getGEO() did not return an ExpressionSet for ", gse_id))
  }
  
  rownames(pheno_data) <- pheno_data$geo_accession
  cat(paste0("Loaded phenotype data for ", gse_id,
             " → Rows: ", nrow(pheno_data), ", Columns: ", ncol(pheno_data), "\n"))
  
  # Restore original working directory
  setwd(old_wd)
  
  return(list(
    pheno_data = pheno_data,
    raw_files = raw_files,
    gse_dir = gse_dir,
    result_dir = result_dir
  ))
}
