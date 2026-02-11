#' Load and normalize expression matrix based on platform type
#'
#' This wrapper dispatches to the appropriate platform-specific function:
#' - Affymetrix: oligo::rma
#' - Illumina: limma::neqc
#' - Agilent: limma::normalizeBetweenArrays
#'
#' @param gse_id GEO Series ID
#' @param pheno_data phenotype data.frame (must include platform_id column)
#' @param raw_files list of raw file paths (e.g. CEL, txt)
#'
#' @return A list containing:
#'   \item{exprs_before}{Log2-transformed raw matrix}
#'   \item{exprs_after}{Normalized expression matrix}
#' @export
normalize_platform <- function(gse_id, pheno_data, raw_files) {
  
  platform_id <- unique(pheno_data$platform_id)
  
  if (length(platform_id) != 1) {
    stop("Error: platform_id must be a single unique value in pheno_data$platform_id")
  }
  
  print_range <- function(mat, label) {
    cat(paste0("Expression matrix value range (", label, ") ", gse_id, 
               ": Min = ", round(min(mat, na.rm = TRUE), 2), 
               ", Max = ", round(max(mat, na.rm = TRUE), 2), "\n"))
  }
  
  # 1. Affymetrix
  if (platform_id %in% c("GPL570", "GPL96", "GPL15207", "GPL17586")) {
    affy_pkg <- list(
      GPL570   = "pd.hgu133plus2.hs.entrezg",
      GPL96    = "pd.hgu133a.hs.entrezg",
      GPL15207 = "pd.primeview.hs.entrezg",
      GPL17586 = "pd.hta20.hs.entrezg"
    )
    pkgname <- affy_pkg[[platform_id]]
    
    # Reorder pheno_data with raw file
    cel_ids    <- basename(raw_files) %>% stringr::str_extract("GSM[0-9]+")
    matched_idx <- match(cel_ids, rownames(pheno_data))
    if (any(is.na(matched_idx))) {
      stop("Some CEL files have no matching rows in pheno_data.")
    }
    pheno_data <- pheno_data[matched_idx, , drop = FALSE]
    
    # Read raw file using 'read.celfiles'
    cel_data <- suppressWarnings(read.celfiles(
      raw_files, pkgname = pkgname,
      phenoData = new("AnnotatedDataFrame", data = pheno_data)
    ))
    
    exprs_before <- log2(intensity(cel_data))
    exprs_after  <- exprs(rma(cel_data))
    
    print_range(exprs_before, "Before")
    print_range(exprs_after, "After")
    
    return(list(exprs_before = exprs_before, exprs_after = exprs_after))
  }
  
  # 2. Illumina
  if (platform_id %in% c("GPL14951", "GPL6884", "GPL6883")) {
    txt_file_path <- raw_files[grepl("non[-_]?normalized.*\\.txt$", raw_files, ignore.case = TRUE)][1]
    
    if (is.na(txt_file_path) || !file.exists(txt_file_path)) {
      stop(paste("No valid Illumina TXT file found for", gse_id))
    }
    
    find_header <- function(path) {
      lines <- readLines(path, warn = FALSE)
      for (i in seq_along(lines)) {
        if (grepl("^ID_REF", lines[[i]])) return(i)
      }
      stop("No valid header found.")
    }
    
    # Read raw file using 'read.delim'
    header_line <- find_header(txt_file_path)
    exprs_raw <- read.delim(txt_file_path, skip = header_line - 1, header = TRUE, quote = "", check.names = FALSE)
    
    col_name <- colnames(exprs_raw)
    detection_cols <- grep("Detection", col_name, value = TRUE)
    exprs_cols <- setdiff(col_name, c("ID_REF", detection_cols))
    id_col <- col_name[1]
    sample_names <- gsub("\\.AVG_Signal$", "", exprs_cols)
    
    # Match pheno_data column with exprssion column
    # Detect which column in pheno_data matches sample_names
    matching_cols <- Filter(function(col) {
      all(sample_names %in% as.character(pheno_data[[col]]))
    }, colnames(pheno_data))
    
    #If no exact match, try partial match
    if (length(matching_cols) != 1) {
      matching_cols <- Filter(function(col) {
        all(sapply(sample_names, function(sn) {
          any(grepl(sn, pheno_data[[col]], fixed = TRUE))
        }))
      }, colnames(pheno_data))
    }
    if (length(matching_cols) != 1) {
      stop("Unable to uniquely identify the matching column in pheno_data.")
    }
    
    # Use matching column to create cleaned sample_id
    if (all(grepl("^\\w+_\\w+$", pheno_data[[matching_cols]]))) {
      pheno_data$sample_id <- as.character(pheno_data[[matching_cols]])  # use directly
    } else {
      pheno_data$sample_id <- sub(".*(\\b\\w+_\\w+)\\s*\\(.*", "\\1", pheno_data[[matching_cols]])
    }
    
    # Match sample_names with sample_id (exact match only)
    matched_idx <- sapply(sample_names, function(sn) {
      matches <- which(pheno_data$sample_id == sn)
      if (length(matches) != 1) {
        stop(paste("Sample", sn, "matches", length(matches), "rows in pheno_data — should match exactly 1."))
      }
      matches
    })
    
    # Reorder pheno_data with sample names
    pheno_data   <- pheno_data[matched_idx, , drop = FALSE]
    exprs_before <- exprs_before[, matched_idx, drop = FALSE]
    exprs_after  <- exprs_after[, matched_idx, drop = FALSE]
    
    # Construct required components for EListRaw object
    E <- as.matrix(exprs_raw[, exprs_cols])
    rownames(pheno_data) <- pheno_data$geo_accession
    colnames(E) <- pheno_data$geo_accession
    
    D <- as.matrix(exprs_raw[, detection_cols])
    probe_id <- exprs_raw[, id_col, drop = FALSE]
    
    # Create raw_elist for 'limma::neqc'
    raw_elist <- new("EListRaw")
    raw_elist$E <- E
    raw_elist$other$Detection <- D
    raw_elist$probe_id <- probe_id
    raw_elist$targets <- pheno_data
    
    rownames(raw_elist$E) <- raw_elist$probe_id[[id_col]]
    
    exprs_before <- log2(raw_elist$E)
    norm_elist <- limma::neqc(raw_elist)
    exprs_after <- norm_elist$E
    
    # Reorder expression matrix with sample names
    exprs_before <- exprs_before[, matched_idx, drop = FALSE]
    exprs_after  <- exprs_after[,  matched_idx, drop = FALSE]
    
    print_range(exprs_before, "Before")
    print_range(exprs_after, "After")
    
    return(list(exprs_before = exprs_before, exprs_after = exprs_after))
  }
  
  # 3. Agilent
  if (platform_id %in% c("GPL20115", "GPL13497")) {
    # Read raw file using 'read.maimages'
    agilent_data <- read.maimages(
      files = raw_files,
      source = "agilent",
      green.only = TRUE,
      other.columns = c("gIsSaturated", "gIsFeatNonUnifOL", "gIsWellAboveBG")
    )
    
    exprs_before <- log2(agilent_data$E)
    exprs_after <- normalizeBetweenArrays(agilent_data, method = "quantile")$E
    
    keep <- agilent_data$genes$ControlType == 0 # Remove rows of contorl probe
    exprs_before <- exprs_before[keep, ]
    exprs_after  <- exprs_after[keep, ]
    agilent_data$genes <- agilent_data$genes[keep, ]
    
    # Match pheno_data column with exprssion column
    sample_names <- sub("_.+", "", basename(raw_files))
    colnames(exprs_before) <- colnames(exprs_after) <- sample_names
    rownames(exprs_before) <- rownames(exprs_after) <- agilent_data$genes$ProbeName
    
    # Match pheno_data column with exprssion column
    # Detect which column in pheno_data matches sample_names
    matching_cols <- Filter(function(col) {
      all(sample_names %in% as.character(pheno_data[[col]]))
    }, colnames(pheno_data))
    
    #If no exact match, try partial match
    if (length(matching_cols) != 1) {
      matching_cols <- Filter(function(col) {
        all(sapply(sample_names, function(sn) {
          any(grepl(sn, pheno_data[[col]], fixed = TRUE))
        }))
      }, colnames(pheno_data))
    }
    if (length(matching_cols) != 1) {
      stop("Unable to uniquely identify the matching column in pheno_data.")
    }
    
    # Use matching column to create cleaned sample_id
    if (all(grepl("^\\w+_\\w+$", pheno_data[[matching_cols]]))) {
      pheno_data$sample_id <- as.character(pheno_data[[matching_cols]]) 
    } else {
      pheno_data$sample_id <- sub(".*(\\b\\w+_\\w+)\\s*\\(.*", "\\1", pheno_data[[matching_cols]])
    }
    
    # Match sample_names with sample_id (exact match only)
    matched_idx <- sapply(sample_names, function(sn) {
      matches <- which(pheno_data$sample_id == sn)
      if (length(matches) != 1) {
        stop(paste("Sample", sn, "matches", length(matches), "rows in pheno_data — should match exactly 1."))
      }
      matches
    })
    
    # Reorder pheno_data and expression matrix with sample names
    pheno_data    <- pheno_data[matched_idx, , drop=FALSE]
    exprs_before  <- exprs_before[, matched_idx, drop=FALSE]
    exprs_after   <- exprs_after[,  matched_idx, drop=FALSE]
    
    print_range(exprs_before, "Before")
    print_range(exprs_after, "After")
    
    return(list(exprs_before = exprs_before, exprs_after = exprs_after))
  }
  
  stop(paste0("Unsupported platform ID '", platform_id, "' for GSE ID: ", gse_id))
}
