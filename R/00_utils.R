# R/00_utils.R
# ============================================================
# IndepMR — Utils: OpenGWAS JWT handling + safe API execution
# upports 4 logical tokens:
#    "tab3", "tab5", "tab7_my", "tab7_x"
# Also accepts a raw JWT string directly (case-sensitive!)
# ============================================================

`%||%` <- function(a, b) {
  if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b
}

# ------------------------------------------------------------
# OpenGWAS token mapping (aliases -> ~/.Renviron variable names)
# ------------------------------------------------------------
token_env_name <- function(token_alias) {
  token_alias <- tolower(trimws(as.character(token_alias %||% "")))
  switch(
    token_alias,
    "tab3"    = "OPENGWAS_JWT_TAB3",
    "tab5"    = "OPENGWAS_JWT_TAB5",
    "tab7_my" = "OPENGWAS_JWT_TAB7_MY",
    "tab7_x"  = "OPENGWAS_JWT_TAB7_X",
    NA_character_
  )
}

# ------------------------------------------------------------
# Resolve logical token name → real JWT
# - If token matches known alias: read from Sys.getenv()
# - Else: assume user passed raw JWT directly (DO NOT lowercase!)
# ------------------------------------------------------------
get_opengwas_jwt <- function(token = "tab3") {
  token_raw <- trimws(as.character(token %||% ""))
  if (token_raw == "") {
    stop("❌ token is empty. Use 'tab3', 'tab5', 'tab7_my', 'tab7_x' or a raw JWT.")
  }
  
  envn <- token_env_name(token_raw)
  
  if (!is.na(envn)) {
    jwt <- Sys.getenv(envn)
    jwt <- as.character(jwt %||% "")
    if (jwt == "") {
      stop(
        "❌ OpenGWAS JWT missing for token alias = ", token_raw, "\n\n",
        "Set these in ~/.Renviron and restart R:\n",
        "  OPENGWAS_JWT_TAB3=...\n",
        "  OPENGWAS_JWT_TAB5=...\n",
        "  OPENGWAS_JWT_TAB7_MY=...\n",
        "  OPENGWAS_JWT_TAB7_X=...\n"
      )
    }
    return(jwt)
  }
  
  # Otherwise treat as raw JWT (case-sensitive)
  token_raw
}

# ------------------------------------------------------------
# Temporarily set OPENGWAS_JWT and evaluate expr safely
# ------------------------------------------------------------
with_opengwas_jwt <- function(jwt, expr) {
  old <- Sys.getenv("OPENGWAS_JWT")
  on.exit(Sys.setenv(OPENGWAS_JWT = old), add = TRUE)
  
  Sys.setenv(OPENGWAS_JWT = jwt)
  force(expr)
}

# ------------------------------------------------------------
# Safe evaluation wrapper (returns either value or error object)
# IMPORTANT: this captures errors without killing Shiny app
# ------------------------------------------------------------
gwas_call_eval <- function(expr_sub, envir, retry = 1) {
  
  if (requireNamespace("ieugwasr", quietly = TRUE)) {
    try(ieugwasr::check_reset(), silent = TRUE)
  }
  
  out <- tryCatch(
    eval(expr_sub, envir = envir),
    error = function(e) e
  )
  
  if (!inherits(out, "error")) return(out)
  
  resp <- attr(out, "response")
  
  if (!is.null(resp) && requireNamespace("ieugwasr", quietly = TRUE)) {
    try(ieugwasr::set_reset(resp), silent = TRUE)
    
    if (retry > 0) {
      try(ieugwasr::check_reset(), silent = TRUE)
      return(gwas_call_eval(expr_sub, envir, retry = retry - 1))
    }
  }
  
  out
}

gwas_call_with_token <- function(token, expr, retry = 1) {
  jwt <- get_opengwas_jwt(token)
  expr_sub <- substitute(expr)
  with_opengwas_jwt(jwt,
                    gwas_call_eval(expr_sub, envir = parent.frame(), retry = retry))
}

# ------------------------------------------------------------
# Utility: result → data.frame (or an error message df)
# ------------------------------------------------------------
as_df_or_msg <- function(x) {
  if (inherits(x, "error")) {
    return(data.frame(Message = conditionMessage(x), stringsAsFactors = FALSE))
  }
  if (is.null(x)) {
    return(data.frame(Message = "No result", stringsAsFactors = FALSE))
  }
  as.data.frame(x)
}

# Optional helper for Shiny: get a readable error string
err_msg <- function(x) {
  if (inherits(x, "error")) return(conditionMessage(x))
  NULL
}

# ============================================================
# Cohort Utilities   
# ============================================================

split_cohorts <- function(x) {
  if (is.null(x) || is.na(x) || x == "") return(character(0))
  
  parts <- unlist(strsplit(as.character(x), "\\|"))
  parts <- trimws(parts)
  parts <- parts[parts != ""]
  
  unique(parts)
}

detect_cohort_overlap <- function(cohort_exp, cohort_out) {
  
  exp_vec <- split_cohorts(cohort_exp)
  out_vec <- split_cohorts(cohort_out)
  
  if (length(exp_vec) == 0 || length(out_vec) == 0) {
    return("UNKNOWN")
  }
  
  overlap <- intersect(exp_vec, out_vec)
  
  if (length(overlap) > 0) {
    return("YES")
  } else {
    return("NO")
  }
}