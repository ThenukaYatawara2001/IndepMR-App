# R/04_search_helpers.R
# ============================================================
# IndepMR — Search helper functions
#
# Provides:
#   friendly_merged_cols()       — rename columns for display
#   fetch_by_accessions()        — filter merged by GCST IDs
#   wide_keyword_to_accessions() — trait keyword → GCST IDs
#   publication_summary()        — one row per PubMed ID
#   publication_details()        — per-study detail table
#
# COHORT column is intentionally excluded from all display
# tables.  Ancestry independence is determined solely from
# EXTRACTED_TERMS via 03_independence_engine.R.
# ============================================================

# Rename core columns to user-friendly display names.
# COHORT is excluded — it is no longer used for any purpose.
friendly_merged_cols <- function(df) {
  # Only rename columns that actually exist
  rename_map <- c(
    PUBMEDID          = "PubMed ID",
    `FIRST AUTHOR`    = "First Author",
    DATE              = "Publication Date",
    JOURNAL           = "Journal",
    STUDY             = "Study Title",
    `DISEASE/TRAIT`   = "Disease / Trait",
    `STUDY ACCESSION` = "Study Accession",
    MAPPED_TRAIT      = "Mapped Trait",
    EXTRACTED_TERMS   = "Population"
  )

  available <- names(rename_map)[names(rename_map) %in% names(df)]
  for (old in available) {
    names(df)[names(df) == old] <- rename_map[old]
  }

  # Drop COHORT column entirely if present
  df[, setdiff(names(df), "COHORT"), drop = FALSE]
}

# Filter merged data frame to rows matching the given accession IDs
fetch_by_accessions <- function(merged, accessions) {
  accessions <- unique(stringr::str_trim(as.character(accessions)))
  accessions <- accessions[accessions != "" & !is.na(accessions)]
  merged %>% dplyr::filter(`STUDY ACCESSION` %in% accessions)
}

# Map a trait keyword to GWAS Catalog study accession IDs via the wide table
wide_keyword_to_accessions <- function(wide, keyword) {
  keyword <- stringr::str_trim(keyword)
  if (is.na(keyword) || keyword == "") return(character(0))

  kw <- stringr::str_to_lower(keyword)

  trait_col <- NULL
  if ("MAPPED_TRAIT"  %in% names(wide)) trait_col <- "MAPPED_TRAIT"
  if (is.null(trait_col) &&
      "mapped_trait" %in% names(wide))  trait_col <- "mapped_trait"

  col_hits  <- names(wide)[stringr::str_detect(
    stringr::str_to_lower(names(wide)), stringr::fixed(kw))]

  accessions <- character(0)

  if (!is.null(trait_col)) {
    trait_rows <- which(stringr::str_detect(
      stringr::str_to_lower(wide[[trait_col]]), stringr::fixed(kw)))
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
  unique(accessions[accessions != "" & !is.na(accessions)])
}

# One-row-per-PubMed-ID summary for the search results table
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

# Per-study detail table for a selected PubMed ID
# COHORT is deliberately excluded from the output.
publication_details <- function(df_rows, pubmed_id, snp_counts = NULL) {
  keep_cols <- intersect(
    c("PUBMEDID", "FIRST AUTHOR", "DATE", "JOURNAL", "STUDY",
      "DISEASE/TRAIT", "MAPPED_TRAIT", "STUDY ACCESSION", "EXTRACTED_TERMS"),
    names(df_rows)
  )

  details <- df_rows %>%
    dplyr::filter(PUBMEDID == pubmed_id) %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(keep_cols))) %>%
    dplyr::arrange(`STUDY ACCESSION`)

  if (!is.null(snp_counts) && nrow(snp_counts) > 0) {
    details <- details %>%
      dplyr::left_join(snp_counts, by = c("STUDY ACCESSION" = "GCST")) %>%
      dplyr::rename(`Valid SNPs` = Valid_SNPs)
  }

  friendly_merged_cols(details)
}
