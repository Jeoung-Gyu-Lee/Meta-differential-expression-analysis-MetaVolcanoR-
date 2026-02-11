#' Filter expression matrix to retain only protein-coding genes
#'
#' Uses Ensembl BioMart via `biomaRt` to identify the biotype of input gene symbols,
#' and returns an expression matrix filtered to protein-coding genes only.
#'
#' @param exprs_matrix A gene expression matrix with rownames as HGNC gene symbols
#' @param gse_id Optional GEO dataset ID for logging (default = NULL)
#'
#' @return A filtered expression matrix (only protein-coding genes)
#' @export
filter_coding <- function(exprs_matrix, gse_id) {
  # Connect to Ensembl
  ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Query gene biotypes
  biotype <- biomaRt::getBM(
    attributes = c("hgnc_symbol", "gene_biotype"),
    filters = "hgnc_symbol",
    values = rownames(exprs_matrix),
    mart = ensembl
  )
  
  # Filter to protein-coding
  protein_genes <- biotype$hgnc_symbol[biotype$gene_biotype == "protein_coding"]
  exprs_filtered <- exprs_matrix[rownames(exprs_matrix) %in% protein_genes, ]
  
  cat(paste0(
    "Filtered by biotype",
    if (!is.null(gse_id)) paste0(" [", gse_id, "]") else "",
    ": Retained only protein_coding genes (n = ", nrow(exprs_filtered), ")\n"
  ))
  
  return(exprs_filtered)
}
