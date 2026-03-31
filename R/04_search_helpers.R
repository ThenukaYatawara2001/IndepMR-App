# R/04_search_helpers.R
# ============================================================
# IndepMR — Search helpers for merged + wide tables
# Minor robustness improvements:
#   - accepts MAPPED_TRAIT or mapped_trait
#   - trims split accessions properly
# ============================================================

friendly_merged_cols <- function(df) {
  df %>%
    dplyr::rename(
      `PubMed ID`        = PUBMEDID,
      `First Author`     = `FIRST AUTHOR`,
      `Publication Date` = DATE,
      `Journal`          = JOURNAL,
      `Study Title`      = STUDY,
      `Disease / Trait`  = `DISEASE/TRAIT`,
      `Study Accession`  = `STUDY ACCESSION`,
      `Mapped Trait`     = MAPPED_TRAIT,
      `Population`       = EXTRACTED_TERMS
    )
}

fetch_by_accessions <- function(merged, accessions) {
  accessions <- unique(stringr::str_trim(as.character(accessions)))
  accessions <- accessions[accessions != "" & !is.na(accessions)]
  merged %>% dplyr::filter(`STUDY ACCESSION` %in% accessions)
}

#
# ====================== Helper overrides ======================
#
# The original definitions of friendly_merged_cols only rename
# columns, discarding the raw COHORT field.  The following
# override retains the COHORT column (if present) to enable cohort
# overlap calculations in downstream modules.  It also ensures
# consistent behaviour by replacing the prior definition when this
# file is sourced.

friendly_merged_cols <- function(df) {
  # Rename core columns to user-friendly names. Use a temporary variable `out`
  # so that we can append additional columns below.
  out <- df %>%
    dplyr::rename(
      `PubMed ID`        = PUBMEDID,
      `First Author`     = `FIRST AUTHOR`,
      `Publication Date` = DATE,
      `Journal`          = JOURNAL,
      `Study Title`      = STUDY,
      `Disease / Trait`  = `DISEASE/TRAIT`,
      `Study Accession`  = `STUDY ACCESSION`,
      `Mapped Trait`     = MAPPED_TRAIT,
      `Population`       = EXTRACTED_TERMS
    )
  # Preserve the raw COHORT column if present so that cohort overlap
  # calculations can be performed downstream. 
  if ("COHORT" %in% names(df)) {
    out$COHORT <- df$COHORT
  }
  out
}

#
# Cohort helpers
#
# The COHORT field is a pipe-separated string listing cohort names.
# These helpers split such strings into individual cohort labels and
# determine whether two studies share at least one cohort.

# Given a COHORT string (e.g. "BioME|GERA|MESA"), return a character
# vector of individual cohort names.  Empty or missing values yield
# an empty vector.  All names are converted to upper-case to allow
# case-insensitive comparisons.
parse_cohort <- function(cohort_str) {
  if (is.null(cohort_str) || is.na(cohort_str) || cohort_str == "") {
    return(character(0))
  }
  parts <- unlist(strsplit(as.character(cohort_str), "\\|", fixed = TRUE))
  parts <- trimws(parts)
  parts <- toupper(parts)
  parts[nchar(parts) > 0]
}

# Compare two cohort strings and report overlap status.
# - If either side has no cohorts, return "UNKNOWN".
# - If there is at least one cohort name in common, return "YES".
# - Otherwise return "NO".
check_overlap <- function(exp_cohort, out_cohort) {
  exp_vec <- parse_cohort(exp_cohort)
  out_vec <- parse_cohort(out_cohort)
  if (length(exp_vec) == 0 || length(out_vec) == 0) {
    return("UNKNOWN")
  }
  if (length(intersect(exp_vec, out_vec)) > 0) {
    return("YES")
  }
  "NO"
}

wide_keyword_to_accessions <- function(wide, keyword) {
  keyword <- stringr::str_trim(keyword)
  if (is.na(keyword) || keyword == "") return(character(0))
  
  kw <- stringr::str_to_lower(keyword)
  
  # allow MAPPED_TRAIT or mapped_trait
  trait_col <- NULL
  if ("MAPPED_TRAIT" %in% names(wide)) trait_col <- "MAPPED_TRAIT"
  if (is.null(trait_col) && "mapped_trait" %in% names(wide)) trait_col <- "mapped_trait"
  
  col_hits <- names(wide)[stringr::str_detect(stringr::str_to_lower(names(wide)), stringr::fixed(kw))]
  accessions <- character(0)
  
  if (!is.null(trait_col)) {
    trait_rows <- which(stringr::str_detect(stringr::str_to_lower(wide[[trait_col]]), stringr::fixed(kw)))
    if (length(trait_rows) > 0) {
      token_cols <- setdiff(names(wide), trait_col)
      vals <- unlist(wide[trait_rows, token_cols, drop = FALSE], use.names = FALSE)
      vals <- vals[!is.na(vals) & vals != ""]
      accessions <- c(accessions, unlist(stringr::str_split(vals, ",")))
    }
  }
  
  if (length(col_hits) > 0) {
    vals <- unlist(wide[, col_hits, drop = FALSE], use.names = FALSE)
    vals <- vals[!is.na(vals) & vals != ""]
    accessions <- c(accessions, unlist(stringr::str_split(vals, ",")))
  }
  
  accessions <- stringr::str_trim(accessions)
  accessions <- accessions[accessions != "" & !is.na(accessions)]
  unique(accessions)
}

publication_summary <- function(df_rows) {
  df_rows %>%
    dplyr::group_by(PUBMEDID) %>%
    dplyr::summarise(
      `First Author`       = dplyr::first(`FIRST AUTHOR`),
      `Publication Date`   = dplyr::first(DATE),
      `Journal`            = dplyr::first(JOURNAL),
      `Study Title`        = dplyr::first(STUDY),
      `# Study Accessions` = dplyr::n_distinct(`STUDY ACCESSION`),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(`# Study Accessions`))
}

publication_details <- function(df_rows, pubmed_id, snp_counts = NULL) {
  details <- df_rows %>%
    dplyr::filter(PUBMEDID == pubmed_id) %>%
    dplyr::distinct(
      PUBMEDID, `FIRST AUTHOR`, DATE, JOURNAL, STUDY,
      `DISEASE/TRAIT`, `MAPPED_TRAIT`, `STUDY ACCESSION`, EXTRACTED_TERMS, COHORT
    ) %>%
    dplyr::arrange(`STUDY ACCESSION`)
  
  if (!is.null(snp_counts) && nrow(snp_counts) > 0) {
    details <- details %>%
      dplyr::left_join(snp_counts, by = c(`STUDY ACCESSION` = "GCST")) %>%
      dplyr::rename(`Valid SNPs` = Valid_SNPs)
  }
  
  friendly_merged_cols(details)
}

# Override cohort helpers to handle duplicates more flexibly
# Split a pipe-delimited COHORT string into tokens, trimming whitespace and converting to upper case
parse_cohort <- function(cohort_str) {
  if (is.null(cohort_str) || is.na(cohort_str) || cohort_str == "") return(character(0))
  tokens <- unlist(strsplit(as.character(cohort_str), "\\|"))
  tokens <- toupper(trimws(tokens))
  tokens[nchar(tokens) > 0]
}

# Determine whether two COHORT strings share any token
check_overlap <- function(exp_cohort, out_cohort) {
  # Handle raw values directly: if either COHORT entry is literally "multiple"
  # (case-insensitive), treat that study as overlapping with any other.
  if (!is.null(exp_cohort) && !is.na(exp_cohort) &&
      tolower(trimws(as.character(exp_cohort))) == "multiple") {
    return("YES")
  }
  if (!is.null(out_cohort) && !is.na(out_cohort) &&
      tolower(trimws(as.character(out_cohort))) == "multiple") {
    return("YES")
  }
  
  exp_tokens <- parse_cohort(exp_cohort)
  out_tokens <- parse_cohort(out_cohort)
  # If either side is empty, you can't determine overlap
  if (length(exp_tokens) == 0 || length(out_tokens) == 0) return("UNKNOWN")
  # If any token overlaps between exposure and outcome, return YES
  if (any(exp_tokens %in% out_tokens)) return("YES")
  # Optionally: mark YES if the outcome's COHORT string itself has duplicates
  if (any(duplicated(out_tokens))) return("YES")
  "NO"
}