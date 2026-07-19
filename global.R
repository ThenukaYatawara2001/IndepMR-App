# global.R
# ============================================================
# IndepMR — Global configuration
# Sourced once by app.R before the UI is built.
# Loads all R modules, defines data paths, and loads DATA.
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
  library(TwoSampleMR)
  library(ieugwasr)
  library(grid)
  library(gridExtra)
})

# ============================================================
# File paths  (relative to project root = getwd())
# ============================================================

APP_DIR          <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
PATH_MERGED_XLSX <- file.path(APP_DIR, "data", "Dataset01.xlsx")
PATH_SNP_RDS     <- file.path(APP_DIR, "data", "Dataset02.rds")
PATH_WIDE_CSV    <- file.path(APP_DIR, "data", "Dataset03.csv")
PATH_CACHE_DIR   <- file.path(APP_DIR, "data", "SNP_CACHE_EXCEL")

# ============================================================
# OpenGWAS token check (informational — user sets token in UI)
# ============================================================

required_tokens <- c(
  "OPENGWAS_JWT_TAB3",
  "OPENGWAS_JWT_TAB5",
  "OPENGWAS_JWT_TAB7_MY",
  "OPENGWAS_JWT_TAB7_X"
)
missing_tokens <- required_tokens[Sys.getenv(required_tokens) == ""]
if (length(missing_tokens) > 0) {
  message("OpenGWAS JWT not found at startup. User can enter token inside the app.")
}

# ============================================================
# Source all R modules
# ============================================================

source(file.path(APP_DIR, "R", "00_utils.R"),             local = FALSE)
source(file.path(APP_DIR, "R", "01_cache_helpers.R"),     local = FALSE)
source(file.path(APP_DIR, "R", "02_load_data.R"),         local = FALSE)
source(file.path(APP_DIR, "R", "03_independence_engine.R"),local = FALSE)
source(file.path(APP_DIR, "R", "04_search_helpers.R"),    local = FALSE)
source(file.path(APP_DIR, "R", "05_search_module.R"),     local = FALSE)
source(file.path(APP_DIR, "R", "06_mr_helpers.R"),        local = FALSE)
source(file.path(APP_DIR, "R", "08_local_instruments.R"), local = FALSE)

if (!exists("searchServer", mode = "function")) {
  stop("searchServer() was not loaded. Check R/05_search_module.R.")
}

# ============================================================
# Load DATA — wrapped so a data error shows inside the app,
# not as a blank "No UI" crash before the UI is registered.
# ============================================================

ensure_dir(PATH_CACHE_DIR)

DATA <- tryCatch(
  load_data(),
  error = function(e) {
    message("DATA LOAD ERROR: ", conditionMessage(e))
    # Return a minimal stub so the UI still renders
    list(
      merged     = data.frame(),
      wide       = data.frame(),
      snp        = NULL,
      snp_counts = NULL,
      load_error = conditionMessage(e)
    )
  }
)
