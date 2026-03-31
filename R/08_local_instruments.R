# R/08_local_instruments.R
# ============================================================
# IndepMR — Local instruments helpers
# Uses SNP table loaded from PATH_SNP_CSV (DATA$snp)
# Returns local instruments in a clean format (SNP, beta, se)
# ============================================================

local_instruments_available <- function(snp_tbl) {
  if (is.null(snp_tbl) || nrow(snp_tbl) == 0) return(FALSE)
  need <- c("GCST", "SNP", "beta.exposure", "se.exposure")
  all(need %in% names(snp_tbl))
}

get_local_instruments <- function(snp_tbl, gcst, p1 = NULL) {
  if (is.null(snp_tbl) || nrow(snp_tbl) == 0) return(NULL)
  if (!all(c("GCST","SNP","beta.exposure","se.exposure") %in% names(snp_tbl))) return(NULL)
  
  gcst <- stringr::str_trim(as.character(gcst))
  if (is.na(gcst) || gcst == "") return(NULL)
  
  df <- snp_tbl %>%
    dplyr::filter(GCST == gcst) %>%
    dplyr::mutate(
      SNP  = as.character(SNP),
      beta = suppressWarnings(as.numeric(.data[["beta.exposure"]])),
      se   = suppressWarnings(as.numeric(.data[["se.exposure"]]))
    )
  
  # Optional: filter by p-value threshold if available
  if (!is.null(p1) && "pval.exposure" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(pval = suppressWarnings(as.numeric(.data[["pval.exposure"]]))) %>%
      dplyr::filter(!is.na(pval), pval <= p1)
  }
  
  df %>%
    dplyr::filter(!is.na(SNP), SNP != "", !is.na(beta), !is.na(se)) %>%
    dplyr::distinct(SNP, .keep_all = TRUE)
}