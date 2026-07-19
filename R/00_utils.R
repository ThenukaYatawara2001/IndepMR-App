# R/00_utils.R
# ============================================================
# IndepMR — Utility functions
# ============================================================

# Null-coalescing operator
`%||%` <- function(a, b) {
  if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b
}

# ── OpenGWAS JWT token management ───────────────────────────

token_env_name <- function(token_alias) {
  token_alias <- tolower(trimws(as.character(token_alias %||% "")))
  switch(token_alias,
    "tab3"    = "OPENGWAS_JWT_TAB3",
    "tab5"    = "OPENGWAS_JWT_TAB5",
    "tab7_my" = "OPENGWAS_JWT_TAB7_MY",
    "tab7_x"  = "OPENGWAS_JWT_TAB7_X",
    NA_character_
  )
}

get_opengwas_jwt <- function(token = "tab3") {
  token_raw <- trimws(as.character(token %||% ""))
  if (token_raw == "")
    stop("token is empty. Use 'tab3','tab5','tab7_my','tab7_x' or a raw JWT.")
  envn <- token_env_name(token_raw)
  if (!is.na(envn)) {
    jwt <- Sys.getenv(envn)
    if (jwt == "")
      stop("OpenGWAS JWT missing for token alias = ", token_raw)
    return(jwt)
  }
  token_raw
}

with_opengwas_jwt <- function(jwt, expr) {
  old <- Sys.getenv("OPENGWAS_JWT")
  on.exit(Sys.setenv(OPENGWAS_JWT = old), add = TRUE)
  Sys.setenv(OPENGWAS_JWT = jwt)
  force(expr)
}

# ── Safe OpenGWAS API call ───────────────────────────────────

gwas_call_eval <- function(expr_sub, envir, retry = 1) {
  if (requireNamespace("ieugwasr", quietly = TRUE))
    try(ieugwasr::check_reset(), silent = TRUE)
  out <- tryCatch(eval(expr_sub, envir = envir), error = function(e) e)
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
  jwt      <- get_opengwas_jwt(token)
  expr_sub <- substitute(expr)
  with_opengwas_jwt(jwt,
    gwas_call_eval(expr_sub, envir = parent.frame(), retry = retry))
}

# ── Data-frame helpers ────────────────────────────────────────

as_df_or_msg <- function(x) {
  if (inherits(x, "error"))
    return(data.frame(Message = conditionMessage(x), stringsAsFactors = FALSE))
  if (is.null(x))
    return(data.frame(Message = "No result", stringsAsFactors = FALSE))
  as.data.frame(x)
}

err_msg <- function(x) {
  if (inherits(x, "error")) conditionMessage(x) else NULL
}
