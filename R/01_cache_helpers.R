# R/01_cache_helpers.R
# ============================================================
# IndepMR — Cache helpers (Excel + filenames + directories)
# Robust temp-write + fallback copy
# Optional "tag" so Tab3/Tab5/Tab7 caches never overwrite each other
# ============================================================

ensure_dir <- function(d) {
  if (is.null(d) || is.na(d) || d == "") return(invisible(FALSE))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(TRUE)
}

safe_filename <- function(x) {
  x <- as.character(x %||% "")
  x <- trimws(x)
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (nchar(x) == 0) x <- "NA"
  x
}

pair_snp_xlsx_path <- function(exp_gcst, out_gcst, tag = NULL) {
  if (!exists("PATH_CACHE_DIR", inherits = TRUE)) stop("PATH_CACHE_DIR is not defined (global.R).")
  ensure_dir(PATH_CACHE_DIR)
  
  tag_part <- ""
  if (!is.null(tag) && !is.na(tag) && trimws(tag) != "") tag_part <- paste0("_", safe_filename(tag))
  
  file.path(
    PATH_CACHE_DIR,
    paste0("SNP_effects", tag_part, "_", safe_filename(exp_gcst), "_to_", safe_filename(out_gcst), ".xlsx")
  )
}

outcome_lookup_xlsx_path <- function(out_gcst, tag = NULL) {
  if (!exists("PATH_CACHE_DIR", inherits = TRUE)) stop("PATH_CACHE_DIR is not defined (global.R).")
  ensure_dir(PATH_CACHE_DIR)
  
  tag_part <- ""
  if (!is.null(tag) && !is.na(tag) && trimws(tag) != "") tag_part <- paste0("_", safe_filename(tag))
  
  file.path(PATH_CACHE_DIR, paste0("Outcome_lookup", tag_part, "_", safe_filename(out_gcst), ".xlsx"))
}

# ------------------------------------------------------------
# Safer write: temp → rename (fallback copy)
# - removes temp file if anything fails
# ------------------------------------------------------------
save_df_to_xlsx <- function(path, sheet = "Sheet1", df) {
  ensure_dir(dirname(path))
  
  tmp <- paste0(path, ".tmp.xlsx")
  if (file.exists(tmp)) suppressWarnings(file.remove(tmp))
  on.exit({ if (file.exists(tmp)) suppressWarnings(file.remove(tmp)) }, add = TRUE)
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheet)
  
  if (is.null(df) || nrow(df) == 0) {
    openxlsx::writeData(wb, sheet, data.frame(Message = "No data", stringsAsFactors = FALSE))
  } else {
    openxlsx::writeData(wb, sheet, df)
  }
  
  openxlsx::saveWorkbook(wb, tmp, overwrite = TRUE)
  
  ok <- suppressWarnings(file.rename(tmp, path))
  if (!isTRUE(ok)) {
    ok2 <- suppressWarnings(file.copy(tmp, path, overwrite = TRUE))
    suppressWarnings(file.remove(tmp))
    ok <- isTRUE(ok2)
  }
  
  if (!isTRUE(ok)) stop("❌ Could not write Excel cache file: ", path)
  invisible(path)
}

read_xlsx_safe <- function(path, sheet = 1) {
  if (is.null(path) || is.na(path) || path == "" || !file.exists(path)) return(NULL)
  
  out <- tryCatch(openxlsx::read.xlsx(path, sheet = sheet), error = function(e) NULL)
  if (is.null(out)) return(NULL)
  
  as.data.frame(out, stringsAsFactors = FALSE)
}