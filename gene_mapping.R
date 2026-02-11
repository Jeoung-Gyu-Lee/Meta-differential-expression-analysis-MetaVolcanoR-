#' Map probe IDs to gene symbols and merge with expression data
#'
#' Performs platform-specific annotation for Affymetrix, Illumina, or Agilent data,
#' and returns merged expression matrix with SYMBOL (and ENTREZID where available).
#'
#' @param exprs_filter Filtered expression matrix (rows = probes, columns = samples)
#' @param pheno_data Phenotype metadata that includes `platform_id` column
#' @param gse_id GEO Series ID (for logging)
#'
#' @return A data.frame with gene annotation merged to expression values
#' @export
gene_mapping <- function(exprs_filter, pheno_data, gse_id) {
  
  # Get platform_id from pheno_data
  platform_id <- unique(pheno_data$platform_id)
  
  if (length(platform_id) != 1) {
    stop("pheno_data must contain exactly one unique platform_id.")
  }
 
  # 1. Affymetrix / Illumina
  affy_illum_mapping <- function(platform_id, exprs_filter, gse_id) {
    db_map <- list(
      GPL570   = "hgu133plus2hsentrezg.db",
      GPL96    = "hgu133ahsentrezg.db",
      GPL14951 = "illuminaHumanv4.db",
      GPL15207 = "primeviewhsentrezg.db",
      GPL17586 = "hta20hsentrezg.db",
      GPL6884  = "illuminaHumanv3.db",
      GPL6883  = "illuminaHumanv3.db"
    )
    
    db_name <- db_map[[platform_id]]
    if (is.null(db_name)) stop(paste("Unsupported platform:", platform_id))
    
    mapping <- AnnotationDbi::select(
      x        = get(db_name),
      keys     = rownames(exprs_filter),
      columns  = c("SYMBOL", "ENTREZID"),
      keytype  = "PROBEID"
    )
    
    cat("Mapping results (before NA removal):", nrow(mapping), "\n")
    cat("Probes with missing SYMBOL:", sum(is.na(mapping$SYMBOL)), "\n")
    
    mapping <- mapping[!is.na(mapping$SYMBOL), ]
    return(mapping)
  }
  
  
  # 2. Agilent
  agilent_mapping <- function(platform_id, gse_id) {
    gpl <- GEOquery::getGEO(platform_id, destdir = getwd())
    annot <- Table(gpl)
    
    probe_col <- if ("ProbeID" %in% colnames(annot)) "ProbeID" else if ("ID" %in% colnames(annot)) "ID" else NULL
    if (is.null(probe_col)) stop("No valid probe ID column found (ProbeID or ID).")
    
    gene_col <- grep("gene.*symbol|symbol.*gene", colnames(annot), ignore.case = TRUE, value = TRUE)
    if (length(gene_col) == 0) stop("No gene symbol column found.")
    
    control_col <- grep("control", colnames(annot), ignore.case = TRUE, value = TRUE)
    if (length(control_col) > 0) {
      colnames(annot)[which(colnames(annot) == control_col[1])] <- "ControlType"
      if (is.character(annot$ControlType)) {
        annot <- annot[annot$ControlType == FALSE, ]
      } else {
        annot <- annot[annot$ControlType == 0, ]
      }
    }
    
    mapping <- annot[, c(probe_col, gene_col[1])]
    colnames(mapping) <- c("PROBEID", "SYMBOL")
    mapping <- mapping[!is.na(mapping$SYMBOL) & mapping$SYMBOL != "", ]
    
    cat("Agilent mapped probes for", gse_id, ":", nrow(mapping), "\n")
    return(mapping)
  }
  
  
  # 3. Choose mapping strategy by platform
  if (platform_id %in% c("GPL570", "GPL96", "GPL14951", "GPL15207", "GPL17586", "GPL6884", "GPL6883")) {
    mapping <- affy_illum_mapping(platform_id, exprs_filter, gse_id)
  } else if (platform_id %in% c("GPL20115", "GPL13497")) {
    mapping <- agilent_mapping(platform_id, gse_id)
  } else {
    stop(paste("Platform", platform_id, "is not supported for annotation."))
  }
  
  
  # 4. Merge mapping with expression matrix
  exprs_filter_df <- as.data.frame(exprs_filter)
  exprs_filter_df$PROBEID <- rownames(exprs_filter_df)
  
  exprs_mapped <- merge(mapping, exprs_filter_df, by = "PROBEID")
  cat("Final merged gene-level matrix for", gse_id, ":", nrow(exprs_mapped), "rows\n")
  
  return(list(mapping = mapping, exprs_mapped = exprs_mapped))
}
