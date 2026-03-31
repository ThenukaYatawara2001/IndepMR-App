# R/02_load_data.R
# ============================================================
# IndepMR — Load data once (merged + wide + SNP instruments + groups)
# Robust column standardization
# Safe numeric parsing and missing-column tolerance
# ============================================================

standardize_snp_tbl <- function(snp_tbl) {
  if (is.null(snp_tbl) || nrow(snp_tbl) == 0) return(NULL)
  
  nm <- names(snp_tbl)
  
  # ---- GCST column ----
  if (!("GCST" %in% nm)) {
    alt <- intersect(nm, c("STUDY_ACCESSION", "Study Accession", "STUDY ACCESSION", "gcst", "GCST"))
    if (length(alt) > 0) names(snp_tbl)[match(alt[1], nm)] <- "GCST"
  }
  
  nm <- names(snp_tbl)
  
  # ---- SNP column ----
  if (!("SNP" %in% nm)) {
    alt2 <- intersect(nm, c("SNPS", "rsid", "RSID", "variant", "snp", "SNP"))
    if (length(alt2) > 0) names(snp_tbl)[match(alt2[1], nm)] <- "SNP"
  }
  
  nm <- names(snp_tbl)
  
  # ---- beta / se / pval exposure ----
  if (!("beta.exposure" %in% nm)) {
    altb <- intersect(nm, c("beta", "BETA", "Beta", "b", "B"))
    if (length(altb) > 0) names(snp_tbl)[match(altb[1], nm)] <- "beta.exposure"
  }
  
  nm <- names(snp_tbl)
  if (!("se.exposure" %in% nm)) {
    alts <- intersect(nm, c("se", "SE", "Se", "stderr", "StdErr"))
    if (length(alts) > 0) names(snp_tbl)[match(alts[1], nm)] <- "se.exposure"
  }
  
  nm <- names(snp_tbl)
  if (!("pval.exposure" %in% nm)) {
    altp <- intersect(nm, c("pval", "P", "p", "p_value", "p.value", "PVAL"))
    if (length(altp) > 0) names(snp_tbl)[match(altp[1], nm)] <- "pval.exposure"
  }
  
  # must have at least GCST + SNP
  if (!all(c("GCST", "SNP") %in% names(snp_tbl))) return(NULL)
  
  snp_tbl %>%
    dplyr::mutate(
      GCST = stringr::str_trim(as.character(GCST)),
      SNP  = stringr::str_trim(as.character(SNP))
    ) %>%
    dplyr::filter(!is.na(GCST), GCST != "", !is.na(SNP), SNP != "")
}

load_data <- function() {
  
  # ---- Load merged + wide ----
  merged <- readxl::read_excel(PATH_MERGED_XLSX)
  wide   <- read.csv(PATH_WIDE_CSV, stringsAsFactors = FALSE, check.names = FALSE)
  
  # ---- Validate merged columns ----
  required_merged <- c(
    "PUBMEDID","FIRST AUTHOR","DATE","JOURNAL","STUDY",
    "DISEASE/TRAIT","MAPPED_TRAIT","STUDY ACCESSION"
  )
  missing_m <- setdiff(required_merged, names(merged))
  if (length(missing_m) > 0) {
    stop(
      paste0(
        "merged_GWAS_OpenGWAS.xlsx missing columns: ",
        paste(missing_m, collapse = ", ")
      )
    )
  }
  
  # Ensure EXTRACTED_TERMS exists
  if (!("EXTRACTED_TERMS" %in% names(merged))) merged$EXTRACTED_TERMS <- NA_character_
  
  # Standardize types
  merged <- merged %>%
    dplyr::mutate(
      PUBMEDID          = as.character(PUBMEDID),
      `FIRST AUTHOR`    = as.character(`FIRST AUTHOR`),
      DATE              = as.character(DATE),
      JOURNAL           = as.character(JOURNAL),
      STUDY             = as.character(STUDY),
      `DISEASE/TRAIT`   = as.character(`DISEASE/TRAIT`),
      `MAPPED_TRAIT`    = as.character(`MAPPED_TRAIT`),
      `STUDY ACCESSION` = stringr::str_trim(as.character(`STUDY ACCESSION`)),
      EXTRACTED_TERMS   = as.character(EXTRACTED_TERMS)
    )
  
  # ---- Load local SNP instruments (optional) ----
  snp_tbl <- NULL
  if (!is.null(PATH_SNP_CSV) && file.exists(PATH_SNP_CSV)) {
    raw <- read.csv(PATH_SNP_CSV, stringsAsFactors = FALSE, check.names = FALSE)
    snp_tbl <- standardize_snp_tbl(raw)
  }
  
  # ---- Load groups (optional) ----
  groups_tbl <- NULL
  if (!is.null(PATH_GROUPS_XLSX) && file.exists(PATH_GROUPS_XLSX)) {
    groups_tbl <- readxl::read_excel(PATH_GROUPS_XLSX)
    
    acc_col <- NULL
    for (c in c("STUDY ACCESSION","STUDY_ACCESSION","Study Accession","GCST")) {
      if (c %in% names(groups_tbl)) { acc_col <- c; break }
    }
    
    if (!is.null(acc_col) && all(c("TOP_CATEGORY_AUTO","SUB_GROUP_AUTO") %in% names(groups_tbl))) {
      groups_tbl <- groups_tbl %>%
        dplyr::mutate(
          GCST = stringr::str_trim(as.character(.data[[acc_col]])),
          TOP_CATEGORY_AUTO = as.character(TOP_CATEGORY_AUTO),
          SUB_GROUP_AUTO    = as.character(SUB_GROUP_AUTO)
        ) %>%
        dplyr::filter(!is.na(GCST), GCST != "") %>%
        dplyr::distinct(GCST, TOP_CATEGORY_AUTO, SUB_GROUP_AUTO)
    } else {
      groups_tbl <- NULL
    }
  }
  
  # ---- SNP counts (optional display stats) ----
  snp_counts <- NULL
  if (!is.null(snp_tbl) && nrow(snp_tbl) > 0 &&
      all(c("beta.exposure", "se.exposure") %in% names(snp_tbl))) {
    
    beta_num <- suppressWarnings(as.numeric(gsub(",", "", snp_tbl$beta.exposure)))
    se_num   <- suppressWarnings(as.numeric(gsub(",", "", snp_tbl$se.exposure)))
    
    tmp_tbl <- snp_tbl %>%
      dplyr::mutate(beta_num = beta_num, se_num = se_num)
    
    snp_counts <- tmp_tbl %>%
      dplyr::filter(!is.na(beta_num), !is.na(se_num)) %>%
      dplyr::group_by(GCST) %>%
      dplyr::summarise(Valid_SNPs = dplyr::n_distinct(SNP), .groups = "drop")
  }
  
  list(
    merged = merged,
    wide   = wide,
    snp    = snp_tbl,
    groups = groups_tbl,
    snp_counts = snp_counts
  )
}