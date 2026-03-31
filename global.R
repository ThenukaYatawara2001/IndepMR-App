# global.R
# ============================================================
# IndepMR — GLOBAL: packages, paths, data loading
# ✅ Updated for FOUR OpenGWAS tokens:
#    TAB3, TAB5, TAB7_MY, TAB7_X
# ============================================================

suppressPackageStartupMessages({
  library(shiny)
  library(DT)
  library(dplyr)
  library(stringr)
  library(readxl)
  library(openxlsx)
  library(DiagrammeR)
  library(sf)
  library(rnaturalearth)
  library(ggplot2)
  library(plotly)
  library(here)
  
  # MR
  library(TwoSampleMR)
  library(ieugwasr)
  
  # PDFs (no LaTeX)
  library(grid)
  library(gridExtra)
  
  # MVMR (optional)
  suppressWarnings(suppressMessages({
    if (requireNamespace("MVMR", quietly = TRUE)) library(MVMR)
  }))
})

# ============================================================
# Helper: Safe data.frame wrapper (prevents Shiny crashes)
# Used in Tab 3, Tab 5, Tab 6 (MVMR), Tab 7 (Mediation)
# ============================================================
safe_df <- function(x, label = "Computation failed") {
  if (inherits(x, "error")) {
    return(data.frame(
      Message = paste0(label, ": ", conditionMessage(x)),
      stringsAsFactors = FALSE
    ))
  }
  if (is.null(x)) {
    return(data.frame(
      Message = paste0(label, ": No output"),
      stringsAsFactors = FALSE
    ))
  }
  if (!is.data.frame(x)) {
    return(data.frame(
      Message = paste0(label, ": Invalid output"),
      stringsAsFactors = FALSE
    ))
  }
  x
}

# ============================================================
# 0) FILE PATHS (EDIT THESE)
# ============================================================
PATH_MERGED_XLSX <- here("data","merged_GWAS_OpenGWAS.xlsx")
PATH_WIDE_CSV    <- here("data","wide_merge.csv")
PATH_SNP_CSV     <- here("data","snp_level_ALL_instruments.csv")
PATH_GROUPS_XLSX <- here("data","paths_AND_SUBGROUP.xlsx")
PATH_CACHE_DIR   <- here("data","SNP_CACHE_EXCEL")

# ============================================================
# OpenGWAS auth (FOUR TOKEN CHECK)
# ============================================================
# You must set these in ~/.Renviron and restart R:
# OPENGWAS_JWT_TAB3=...
# OPENGWAS_JWT_TAB5=...
# OPENGWAS_JWT_TAB7_MY=...
# OPENGWAS_JWT_TAB7_X=...

required_tokens <- c(
  "OPENGWAS_JWT_TAB3",
  "OPENGWAS_JWT_TAB5",
  "OPENGWAS_JWT_TAB7_MY",
  "OPENGWAS_JWT_TAB7_X"
)

missing_tokens <- required_tokens[Sys.getenv(required_tokens) == ""]
if (length(missing_tokens) > 0) {
  stop(
    "❌ Missing OpenGWAS JWT(s): ", paste(missing_tokens, collapse = ", "), "\n",
    "Add them to ~/.Renviron and restart R.\n\n",
    "Example:\n",
    "OPENGWAS_JWT_TAB3=xxxxxxxx\n",
    "OPENGWAS_JWT_TAB5=xxxxxxxx\n",
    "OPENGWAS_JWT_TAB7_MY=xxxxxxxx\n",
    "OPENGWAS_JWT_TAB7_X=xxxxxxxx\n"
  )
}

# ============================================================
# Source all R scripts
# ============================================================
# (Order matters: utils first)
source(here("R", "00_utils.R"))
source(here("R", "01_cache_helpers.R"))
source(here("R", "02_load_data.R"))
source(here("R", "03_independence_engine.R"))
source(here("R", "04_search_helpers.R"))
source(here("R", "05_search_module.R"))
source(here("R", "06_mr_helpers.R"))
source(here("R", "07_pdf_helpers.R"))
source(here("R", "08_local_instruments.R"))
source(here("R", "server_original.R")) 


# ============================================================
# Load DATA once
# ============================================================
DATA <- load_data()
ensure_dir(PATH_CACHE_DIR)