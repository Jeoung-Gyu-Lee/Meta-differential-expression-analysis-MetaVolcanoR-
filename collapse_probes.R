#' Resolve probe-to-gene mapping ambiguity and collapse redundant probes
#'
#' Handles 1:N and N:1 probe-to-gene mappings. Removes probes mapped to multiple genes,
#' and collapses multiple probes mapped to the same gene using max variance strategy.
#'
#' @param exprs_mapped A data.frame containing mapped expression data including PROBEID, SYMBOL, and expression columns
#' @param gse_id GEO dataset ID for logging and messages
#'
#' @return A matrix of collapsed expression values (genes x samples)
#' @export
collapse_probes <- function(exprs_mapped, gse_id) {
  # Check required columns
  if (!all(c("PROBEID", "SYMBOL") %in% colnames(exprs_mapped))) {
    stop("exprs_mapped must contain 'PROBEID' and 'SYMBOL' columns")
  }
  
  # Remove 1:N ambiguous mappings
  # Remove probes that map to more than two genes
  probe_dup_gene <- mapping %>%
    group_by(PROBEID) %>%
    filter(dplyr::n_distinct(SYMBOL) > 1) %>%
    pull(PROBEID) %>%
    unique()
  cat("Number of rows before removing 1:N mappings:", length(probe_dup_gene), "\n")
  
  # Both mapping and exprs_mapped should be filtered
  mapping <- mapping[!mapping$PROBEID %in% probe_dup_gene, ]
  exprs_mapped <- exprs_mapped[exprs_mapped$PROBEID %in% mapping$PROBEID, ]
  
  cat(paste0("Remaining unique probe mappings for ", gse_id, ": ", nrow(exprs_mapped), "\n"))
  
  # Collapse N:1 mappings using collapseRows
  expr_collapse <- exprs_mapped %>%
    dplyr::select(-any_of(c("PROBEID", "SYMBOL", "ENTREZID"))) %>%
    as.matrix()
  rownames(expr_collapse) <- exprs_mapped$PROBEID
  
  set.seed(123) # ensure reproducibility of probe collapsing
  
  collapsed_gene <- collapseRows(
    dat = expr_collapse,
    rowGroup = exprs_mapped$SYMBOL,
    rowID = exprs_mapped$PROBEID,
    method = "maxRowVariance",
    connectivityBasedCollapsing = FALSE
  )
  
  exprs_collapsed <- collapsed_gene$datETcollapsed
  cat(paste0("Remaining unique probes after collapseRows for ", gse_id, ": ", nrow(exprs_collapsed), "\n"))

  return(exprs_collapsed)
}