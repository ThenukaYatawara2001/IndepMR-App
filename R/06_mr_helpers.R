# ============================================================
# IndepMR — MR Helpers (Corrected & Safe Version)
# ============================================================

# ============================================================
# MR method labels
# ============================================================

mr_method_choices <- function() {
  c(
    "Wald ratio"                     = "mr_wald_ratio",
    "Inverse variance weighted (IVW)" = "mr_ivw",
    "MR Egger"                       = "mr_egger_regression",
    "Weighted median"                = "mr_weighted_median",
    "Simple median"                  = "mr_simple_median",
    "Weighted mode"                  = "mr_weighted_mode",
    "Simple mode"                    = "mr_simple_mode",
    "IVW (multiplicative RE)"        = "mr_ivw_mre",
    "IVW (fixed effects)"            = "mr_ivw_fe"
  )
}

# ============================================================
# Harmonise mapping
# ============================================================

harmonise_action_value <- function(x) {
  if (identical(x, "pos")) return(1)
  if (identical(x, "align")) return(2)
  if (identical(x, "exclude")) return(3)
  2
}

# ============================================================
# P-value helpers
# ============================================================

sig_label <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  if (is.na(p)) return("Unknown")
  if (p < 0.001) return("*** (p<0.001)")
  if (p < 0.01)  return("** (p<0.01)")
  if (p < 0.05)  return("* (p<0.05)")
  "NS (p>=0.05)"
}

find_p_cols <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(character(0))
  nm <- names(df)
  nm_lower <- tolower(nm)
  
  idx <- which(grepl("^pr\\(|p\\s*value|p_value|pval", nm_lower))
  if (length(idx) > 0) return(nm[idx[1]])
  character(0)
}

round_df_display <- function(df, digits = 4) {
  if (is.null(df) || nrow(df) == 0) return(df)
  out <- df
  num_cols <- names(out)[vapply(out, is.numeric, logical(1))]
  for (cc in num_cols) out[[cc]] <- round(out[[cc]], digits)
  out
}

add_significance_column <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  pcol <- find_p_cols(df)
  if (length(pcol) == 0) {
    df$Significance <- "Unknown"
    return(df)
  }
  df$Significance <- vapply(df[[pcol]], sig_label, character(1))
  df
}

# ============================================================
# Allele alignment (SAFE VERSION)
# ============================================================

is_palindromic <- function(a1, a2) {
  paste0(toupper(a1), toupper(a2)) %in% c("AT","TA","CG","GC")
}

align_to_ref <- function(beta, ea, oa, ref_ea, ref_oa) {
  
  ea <- toupper(ea); oa <- toupper(oa)
  ref_ea <- toupper(ref_ea); ref_oa <- toupper(ref_oa)
  
  same <- ea == ref_ea & oa == ref_oa
  flip <- ea == ref_oa & oa == ref_ea
  keep <- same | flip
  
  beta2 <- beta
  beta2[flip] <- -beta2[flip]
  
  list(beta = beta2, flipped = flip, keep = keep)
}

# ============================================================
# MVMR IVW table cleaner
# ============================================================

process_ivw_table <- function(df) {
  
  df <- as.data.frame(df)
  
  if (nrow(df) == 0 || "Message" %in% names(df)) {
    return(data.frame(Message = "No IVW MVMR results."))
  }
  
  if (!is.null(rownames(df)) && any(rownames(df) != "")) {
    df$Exposure <- rownames(df)
    rownames(df) <- NULL
  } else {
    df$Exposure <- paste0("Exposure", seq_len(nrow(df)))
  }
  
  rename_map <- c(
    "Estimate"   = "b",
    "Std. Error" = "se",
    "Pr(>|t|)"   = "pval"
  )
  
  for (nm in names(rename_map)) {
    if (nm %in% names(df)) {
      names(df)[names(df) == nm] <- rename_map[[nm]]
    }
  }
  
  round_df_display(df, 4)
}

# ============================================================
# Token-aware in-memory caches
# ============================================================

OUTCOMES_CACHE <- new.env(parent = emptyenv())
ID_CACHE       <- new.env(parent = emptyenv())

token_key <- function(token = "tab3") {
  token <- tolower(trimws(as.character(token %||% "")))
  if (token == "") token <- "tab3"
  paste0("TKN__", token)
}

# ============================================================
# Get available outcomes (cached per token)
# ============================================================

get_available_outcomes_cached <- function(token = "tab3") {
  
  tkey <- token_key(token)
  obj  <- paste0("df__", tkey)
  
  if (exists(obj, envir = OUTCOMES_CACHE, inherits = FALSE)) {
    return(get(obj, envir = OUTCOMES_CACHE))
  }
  
  df <- gwas_call_with_token(token, TwoSampleMR::available_outcomes())
  
  if (inherits(df, "error")) return(df)
  
  assign(obj, df, envir = OUTCOMES_CACHE)
  df
}

# ============================================================
# Resolve OpenGWAS ID (SAFE VERSION)
# ============================================================

resolve_opengwas_id <- function(gcst, token = "tab3") {
  
  gcst <- stringr::str_trim(as.character(gcst))
  if (is.na(gcst) || gcst == "") return(NA_character_)
  
  tkey <- token_key(token)
  key  <- paste0("ID__", tkey, "__", gcst)
  
  if (exists(key, envir = ID_CACHE, inherits = FALSE)) {
    return(get(key, envir = ID_CACHE))
  }
  
  all_studies <- get_available_outcomes_cached(token)
  
  if (inherits(all_studies, "error") || is.null(all_studies) || nrow(all_studies) == 0) {
    assign(key, NA_character_, envir = ID_CACHE)
    return(NA_character_)
  }
  
  if ("id" %in% names(all_studies)) {
    
    # exact match
    if (gcst %in% all_studies$id) {
      assign(key, gcst, envir = ID_CACHE)
      return(gcst)
    }
    
    # match GCST inside id
    hit <- all_studies$id[grepl(paste0("\\b", gcst, "\\b"), all_studies$id)]
    if (length(hit) > 0) {
      assign(key, hit[1], envir = ID_CACHE)
      return(hit[1])
    }
  }
  
  assign(key, NA_character_, envir = ID_CACHE)
  NA_character_
}

# ============================================================
# Hash helper (SAFE cache key)
# ============================================================

hash_vec <- function(x) {
  x <- sort(unique(as.character(x)))
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(x))
  }
  paste0(length(x), "_", substr(paste0(x, collapse="|"), 1, 200))
}

# ============================================================
# Chunk SNPs
# ============================================================

chunk_snps <- function(snps, size = 100) {
  snps <- unique(stringr::str_trim(as.character(snps)))
  snps <- snps[snps != "" & !is.na(snps)]
  split(snps, ceiling(seq_along(snps) / size))
}

# ============================================================
# Safe chunked extraction with correct cache hashing
# ============================================================

extract_outcome_chunked_cached <- function(
    snps,
    outcome_id,
    token,
    cache_dir,
    proxies = TRUE,
    rsq = 0.8,
    palindromes = TRUE,
    maf_threshold = 0.3
) {
  
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  
  chunks <- chunk_snps(snps, size = 100)
  out_list <- list()
  
  for (i in seq_along(chunks)) {
    
    subdir <- file.path(cache_dir, outcome_id)
    dir.create(subdir, showWarnings = FALSE, recursive = TRUE)
    
    chunk_hash <- hash_vec(c(
      chunks[[i]],
      outcome_id,
      proxies, rsq, palindromes, maf_threshold
    ))
    
    cache_file <- file.path(subdir, paste0("chunk_", i, "_", chunk_hash, ".rds"))
    
    if (file.exists(cache_file)) {
      out_list[[i]] <- readRDS(cache_file)
      next
    }
    
    res <- tryCatch(
      gwas_call_with_token(
        token,
        TwoSampleMR::extract_outcome_data(
          snps = chunks[[i]],
          outcomes = outcome_id,
          proxies = proxies,
          rsq = rsq,
          palindromes = palindromes,
          maf_threshold = maf_threshold
        )
      ),
      error = function(e) e
    )
    
    if (inherits(res, "error")) {
      msg <- conditionMessage(res)
      
      if (grepl("rate limit|429|too many requests", tolower(msg))) {
        reset_time_str <- stringr::str_extract(
          msg,
          "\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}\\.\\d+"
        )
        return(structure(
          list(error = TRUE, message = msg, reset_time = reset_time_str),
          class = "opengwas_limit"
        ))
      }
      
      return(res)
    }
    
    saveRDS(res, cache_file)
    out_list[[i]] <- res
  }
  
  dplyr::bind_rows(out_list)
}