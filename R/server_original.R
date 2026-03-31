# IndepMR — Server Logic
#
# This file contains the server logic for the IndepMR Shiny application.
extract_countries_from_population <- function(pop_string) {
  
  if (is.null(pop_string) || is.na(pop_string) || pop_string == "")
    return(NULL)
  
  pops <- unlist(strsplit(tolower(pop_string), ","))
  pops <- trimws(pops)
  
  countries <- unlist(population_to_country[pops])
  countries <- unique(countries[!is.na(countries)])
  
  return(countries)
}

get_trait_label <- function(gcst_id, data_tbl) {
  
  if (is.null(gcst_id) || gcst_id == "") return(NA)
  
  trait_name <- data_tbl %>%
    dplyr::filter(`STUDY ACCESSION` == gcst_id) %>%
    dplyr::pull(MAPPED_TRAIT) %>%
    unique()
  
  trait_name <- trait_name[!is.na(trait_name)][1]
  
  if (is.na(trait_name) || length(trait_name) == 0) {
    trait_name <- "Unknown trait"
  }
  
  paste0(trait_name, "  || id: ebi-a-", gcst_id)
}
server <- function(input, output, session) {
  
  opengwas_reset_rv <- reactiveVal(NULL)
  
  countdown_tick <- reactiveTimer(1000)
  
  handle_opengwas_limit <- function(obj, rv) {
    if (inherits(obj, "opengwas_limit")) {
      
      if (!is.null(obj$reset_time)) {
        rv(as.POSIXct(obj$reset_time, tz = "Asia/Colombo"))
      }
      
      shiny::showNotification(
        "⛔ OpenGWAS rate limit reached. Countdown started.",
        type = "warning",
        duration = NULL
      )
      return(TRUE)
    }
    FALSE
  }
  
  # -----------------------------
  # OpenGWAS countdown output
  # -----------------------------
  output$opengwas_countdown <- renderText({
    countdown_tick()
    
    rt <- opengwas_reset_rv()
    
    # Case 1: No rate limit active
    if (is.null(rt)) return("")
    
    now <- Sys.time()
    secs_left <- as.numeric(difftime(rt, now, units = "secs"))
    
    # Case 2: Rate limit but reset time unknown
    if (!is.finite(secs_left)) {
      return(
        paste(
          "⛔ OpenGWAS rate limit reached.",
          "\nReset time was not provided by OpenGWAS.",
          "\nPlease wait a few minutes before retrying."
        )
      )
    }
    
    # Case 3: Reset passed
    if (secs_left <= 0) {
      opengwas_reset_rv(NULL)
      return("OpenGWAS allowance reset. You can continue.")
    }
    
    # Case 4: Valid countdown
    mins <- floor(secs_left / 60)
    secs <- floor(secs_left %% 60)
    
    paste0(
      "⏳ OpenGWAS allowance resets in ",
      mins, " min ", secs, " sec\n",
      "Reset time (Sri Lanka): ",
      format(rt, "%Y-%m-%d %H:%M:%S")
    )
  })
  
  render_pair_box <- function(lhs, rhs, title = NULL) {
    renderUI({
      req(lhs(), rhs())
      
      div(
        style = "
        padding: 16px 20px;
        border-radius: 22px;
        background: rgba(8, 18, 40, 0.75);
        border: 1px solid rgba(255,255,255,0.15);
        box-shadow: inset 0 0 0 1px rgba(255,255,255,0.05);
      ",
        
        div(
          style = "display:flex; align-items:center; justify-content:space-between;",
          
          div(
            style = "font-size:24px; font-weight:900; letter-spacing:0.5px;",
            paste0(lhs(), " \u2192 ", rhs())
          ),
          
          div(
            style = "font-size:20px; opacity:0.7;",
            HTML("&#9662;")   # ▼ arrow
          )
        )
      )
    })
  }
  # Helper: Busy modal
  show_busy_modal <- function(title="Performing...", msg="Please wait. This can take a while.") {
    showModal(modalDialog(
      title = title,
      div(class = "spinWrap",
          div(class = "spinner"),
          div(tags$b(msg))
      ),
      footer = NULL,
      easyClose = FALSE
    ))
  }
  # ==========================================================
  # TAB 7 PATH CACHE (avoid repeated OpenGWAS calls)
  # ==========================================================
  TAB7_PATH_CACHE <- new.env(parent = emptyenv())
  
  .path_key <- function(exp_gcst, out_gcst, token_alias) {
    paste0(token_alias, " :: ", str_trim(exp_gcst), "  →  ", str_trim(out_gcst))
  }
  .cache_get <- function(key) {
    if (exists(key, envir = TAB7_PATH_CACHE, inherits = FALSE)) get(key, envir = TAB7_PATH_CACHE) else NULL
  }
  .cache_set <- function(key, obj) assign(key, obj, envir = TAB7_PATH_CACHE)
  
  # ---- SAFE IVW runner (NO input$... used inside!) ----
  run_ivw_pair_tab7 <- function(
    exp_gcst, out_gcst, token_alias,
    p1 = 5e-8,
    do_clump = TRUE, r2 = 0.01, kb = 10000,
    proxies = TRUE, proxy_rsq = 0.8,
    allow_palindromes = TRUE, maf_threshold = 0.30,
    harm_action = "align",
    do_steiger = TRUE
  ) {
    
    key <- .path_key(exp_gcst, out_gcst, token_alias)
    
    e_id <- resolve_opengwas_id(exp_gcst, token_alias)
    o_id <- resolve_opengwas_id(out_gcst, token_alias)
    if (is.na(e_id) || is.na(o_id)) {
      out <- list(error = paste0("Could not resolve OpenGWAS IDs: ",
                                 exp_gcst, "→", out_gcst, " | e_id=", e_id, " o_id=", o_id))
      return(out)
    }
    
    ex_dat <- gwas_call_with_token(
      token_alias,
      extract_instruments(outcomes = e_id, p1 = p1, clump = do_clump, r2 = r2, kb = kb)
    )
    if (inherits(ex_dat, "error") || is.null(ex_dat) || nrow(as.data.frame(ex_dat)) == 0) {
      out <- list(error = paste0("extract_instruments failed: ",
                                 if (inherits(ex_dat,"error")) conditionMessage(ex_dat) else "No instruments returned."))
      .cache_set(key, out)
      return(out)
    }
    
    out_dat <- extract_outcome_chunked_cached(
      snps         = ex_dat$SNP,
      outcome_id  = o_id,
      token       = token_alias,
      cache_dir = file.path(
        tempdir(),
        paste0("mr_mediation_", as.integer(Sys.time()), "_", sample(1e6,1))
      ),
      proxies     = proxies,
      rsq         = proxy_rsq,
      palindromes = allow_palindromes,
      maf_threshold = maf_threshold
    )
    
    if (inherits(out_dat, "error") || is.null(out_dat) || nrow(as.data.frame(out_dat)) == 0) {
      out <- list(error = paste0("extract_outcome_data failed: ",
                                 if (inherits(out_dat,"error")) conditionMessage(out_dat) else "No outcome rows returned."))
      .cache_set(key, out)
      return(out)
    }
    
    hdat <- tryCatch(
      harmonise_data(exposure_dat = ex_dat, outcome_dat = out_dat, action = harmonise_action_value(harm_action)),
      error = function(e) e
    )
    if (inherits(hdat, "error") || is.null(hdat) || nrow(as.data.frame(hdat)) == 0) {
      out <- list(error = paste0("harmonise_data failed: ",
                                 if (inherits(hdat,"error")) conditionMessage(hdat) else "No harmonised rows."))
      .cache_set(key, out)
      return(out)
    }
    
    mr_ivw <- tryCatch(mr(hdat, method_list = c("mr_ivw")), error = function(e) e)
    if (inherits(mr_ivw, "error")) {
      out <- list(error = paste0("mr(ivw) failed: ", conditionMessage(mr_ivw)))
      .cache_set(key, out)
      return(out)
    }
    
    steiger_tbl <- NULL
    if (isTRUE(do_steiger)) {
      st <- tryCatch(directionality_test(hdat), error = function(e) e)
      if (!inherits(st, "error")) steiger_tbl <- as.data.frame(st)
    }
    
    out <- list(
      mr_ivw = as.data.frame(mr_ivw),
      steiger = steiger_tbl,
      harmonised = hdat
    )
    .cache_set(key, out)
    out
  }
  
  # --------------------------------------------------------------------
  # MODULES (searchServer for Exposure and Outcome)
  # Pass snp_counts = DATA$snp_counts as required by new design
  # --------------------------------------------------------------------
  
  # Exposure selection module
  exposure_mod <- searchServer(
    id = "exposure",
    merged = DATA$merged,
    wide = DATA$wide,
    snp_counts = DATA$snp_counts,
    status_fun_reactive = NULL,
    filter_only_indep_reactive = NULL,
    enable_selection = TRUE
  )
  
  # Reactive for outcome independence vs selected exposures
  outcome_status_fun <- reactive({
    exp1 <- unique(str_trim(as.character(exposure_mod$selected_gcst())))
    exp1 <- exp1[exp1 != "" & !is.na(exp1)]
    
    function(out_gcst_vec) {
      vapply(out_gcst_vec, function(oid) {
        if (length(exp1) == 0 || is.na(oid) || oid == "") return("Unknown")
        ok <- is_independent_from_all(oid, exp1[1], DATA$merged)
        if (is.na(ok)) "Unknown" else if (ok) "Independent" else "Not independent"
      }, character(1))
    }
  })
  
  outcome_only_indep <- reactive({ isTRUE(input$show_only_independent) })
  
  # Reactive tracking the second checkbox for outcome overlap filtering
  outcome_only_non_overlap <- reactive({ isTRUE(input$show_non_overlap_outcome) })
  
  # Function to compute cohort overlap between selected exposures and candidate outcomes.
  # Returns a function that accepts a vector of GCST IDs and returns "YES", "NO"
  # or "UNKNOWN" for each.
  outcome_overlap_fun <- reactive({
    # Identify selected exposure GCST IDs
    exp_gcst <- unique(stringr::str_trim(as.character(exposure_mod$selected_gcst())))
    exp_gcst <- exp_gcst[exp_gcst != "" & !is.na(exp_gcst)]
    if (length(exp_gcst) == 0) return(NULL)
    # Build a union of cohort names across all selected exposures using the COHORT column.
    # For each exposure GCST, collect all non-empty COHORT strings and parse each one;
    # then flatten and deduplicate.  This avoids relying on only the first non-empty
    # entry when there are multiple records per GCST.
    exp_cohorts <- unique(unlist(lapply(exp_gcst, function(gid) {
      rows <- which(DATA$merged$`STUDY ACCESSION` == gid)
      if (length(rows) == 0) return(character(0))
      co_strs <- DATA$merged$COHORT[rows]
      co_strs <- co_strs[!is.na(co_strs) & co_strs != ""]
      if (length(co_strs) == 0) return(character(0))
      unique(unlist(lapply(co_strs, parse_cohort)))
    })))
    if (length(exp_cohorts) == 0) {
      return(function(out_gcst_vec) rep("UNKNOWN", length(out_gcst_vec)))
    }
    # Closure computing overlap for each candidate outcome GCST
    function(out_gcst_vec) {
      vapply(out_gcst_vec, function(oid) {
        rows <- which(DATA$merged$`STUDY ACCESSION` == oid)
        if (length(rows) == 0) return("UNKNOWN")
        co_strs <- DATA$merged$COHORT[rows]
        co_strs <- co_strs[!is.na(co_strs) & co_strs != ""]
        if (length(co_strs) == 0) return("UNKNOWN")
        out_coh <- unique(unlist(lapply(co_strs, parse_cohort)))
        # Treat presence of 'MULTIPLE' in either candidate or reference as overlapping
        if ("MULTIPLE" %in% out_coh || "MULTIPLE" %in% exp_cohorts) {
          return("YES")
        }
        if (length(intersect(exp_cohorts, out_coh)) > 0) "YES" else "NO"
      }, character(1))
    }
  })
  
  # Outcome selection module (passes snp_counts).  We supply the
  # overlap_fun_reactive and filter_only_non_overlap_reactive arguments
  # to support same‑ancestry and non‑overlap filtering in the search module.
  outcome_mod <- searchServer(
    id = "outcome",
    merged = DATA$merged,
    wide = DATA$wide,
    snp_counts = DATA$snp_counts,
    status_fun_reactive = outcome_status_fun,
    filter_only_indep_reactive = outcome_only_indep,
    overlap_fun_reactive = outcome_overlap_fun,
    filter_only_non_overlap_reactive = outcome_only_non_overlap,
    enable_selection = TRUE
  )
  
  observeEvent(outcome_mod$selected_gcst(), {
    
    ids <- outcome_mod$selected_gcst()
    if (length(ids) == 0) return()
    
    outcome_id <- ids[1]
    
    row <- DATA$merged %>%
      dplyr::filter(`STUDY ACCESSION` == outcome_id)
    
    if (nrow(row) == 0) return()
    
    outcome_type <- row$`Data Type`[1]
    
    if (!is.na(outcome_type) && outcome_type == "Binary") {
      
      showNotification(
        paste0(
          "⚠ Outcome ", outcome_id, " (",
          row$MAPPED_TRAIT[1],
          ") is binary. Effects are in log-odds scale. ",
          "Use exp(beta) to interpret as Odds Ratio."
        ),
        type = "warning",
        duration = 8
      )
      
    }
    
  })
  
  mr_rv <- reactiveValues(
    status = "Select exposure(s) and outcome(s), then click Run MR.",
    pairs = character(0),
    results_by_pair = list(),
    snp_effects_by_pair = list()
  )
  observeEvent(
    list(exposure_mod$selected_gcst(), outcome_mod$selected_gcst()),
    {
      
      mr_rv$results_by_pair <- list()
      mr_rv$snp_effects_by_pair <- list()
      mr_rv$pairs <- character(0)
      mr_rv$status <- "Selection changed. Click Run MR."
    },
    ignoreInit = TRUE
  )
  # ==========================================================
  # TAB 7 — Reference IDs (Exposure1 + Outcome)
  # ==========================================================
  tab7_ref_ids <- reactive({
    
    # Updated reference set for mediator independence: only use the first selected exposure GCST.
    exp1 <- unique(str_trim(as.character(exposure_mod$selected_gcst())))[1]
    ref <- exp1
    ref <- ref[!is.na(ref) & ref != ""]
    if (length(ref) == 0) return(character(0))
    ref
  })
  # --------------------------------------------------------------------
  # TAB 3: MR ENGINE
  # Manage caching of SNP effect tables and results for each pair
  # --------------------------------------------------------------------
  OUTCOME_SNP_CACHE <- new.env(parent = emptyenv())
  
  output$mr_status <- renderText({ mr_rv$status })
  
  # Helper to build per-pair SNP effect table from harmonised data
  build_snp_effects_table <- function(hdat) {
    h <- as.data.frame(hdat)
    need <- c("SNP", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome")
    if (!all(need %in% names(h))) return(NULL)
    h %>%
      transmute(
        SNP         = as.character(SNP),
        Exposure_beta = as.numeric(beta.exposure),
        Exposure_se   = as.numeric(se.exposure),
        Outcome_beta  = as.numeric(beta.outcome),
        Outcome_se    = as.numeric(se.outcome)
      ) %>%
      filter(!is.na(SNP), SNP != "") %>%
      distinct(SNP, .keep_all = TRUE)
  }
  
  # Caching functions for SNP effects and outcome lookups
  persist_outcome_lookup <- function(out_gcst, df_lookup) {
    if (is.null(df_lookup) || nrow(df_lookup) == 0) return(invisible())
    df2 <- df_lookup %>%
      transmute(
        SNP         = as.character(SNP),
        Outcome_beta = as.numeric(Outcome_beta),
        Outcome_se   = as.numeric(Outcome_se)
      ) %>%
      filter(!is.na(SNP), SNP != "", !is.na(Outcome_beta), !is.na(Outcome_se)) %>%
      distinct(SNP, .keep_all = TRUE)
    save_df_to_xlsx(outcome_lookup_xlsx_path(out_gcst), sheet = "Outcome_lookup", df = df2)
  }
  
  persist_pair_snp_effects <- function(exp_gcst, out_gcst, snp_eff) {
    if (is.null(snp_eff) || nrow(snp_eff) == 0) return(invisible())
    save_df_to_xlsx(pair_snp_xlsx_path(exp_gcst, out_gcst), sheet = "SNP_effects", df = snp_eff)
  }
  
  update_outcome_cache <- function(out_gcst, snp_eff_pair) {
    if (is.null(snp_eff_pair) || nrow(snp_eff_pair) == 0) return(invisible())
    key <- paste0("OUT__", out_gcst)
    new_rows <- snp_eff_pair %>%
      transmute(SNP, Outcome_beta, Outcome_se) %>%
      filter(!is.na(SNP), SNP != "", !is.na(Outcome_beta), !is.na(Outcome_se)) %>%
      distinct(SNP, .keep_all = TRUE)
    
    if (!exists(key, envir = OUTCOME_SNP_CACHE, inherits = FALSE)) {
      assign(key, new_rows, envir = OUTCOME_SNP_CACHE)
      persist_outcome_lookup(out_gcst, new_rows)
      return(invisible())
    }
    old <- get(key, envir = OUTCOME_SNP_CACHE, inherits = FALSE)
    comb <- bind_rows(old, new_rows) %>% distinct(SNP, .keep_all = TRUE)
    assign(key, comb, envir = OUTCOME_SNP_CACHE)
    persist_outcome_lookup(out_gcst, comb)
  }
  
  # Main MR event: run two-sample MR across selected pairs
  observeEvent(input$run_mr, {
    run_cache_dir <- file.path(
      tempdir(),
      paste0("mr_tab3_", as.integer(Sys.time()), "_", sample(1e6,1))
    )
    
    dir.create(run_cache_dir, recursive = TRUE, showWarnings = FALSE)
    show_busy_modal("Performing MR Exposure 1 \u2192 Outcome...", "Running MR Exposure 1 \u2192 Outcome now…")
    on.exit({ removeModal() }, add = TRUE)
    
    # Gather selected exposure and outcome GCST IDs
    exp_gcst <- unique(str_trim(as.character(exposure_mod$selected_gcst())))
    exp_gcst <- exp_gcst[exp_gcst != "" & !is.na(exp_gcst)]
    out_gcst <- unique(str_trim(as.character(outcome_mod$selected_gcst())))
    out_gcst <- out_gcst[out_gcst != "" & !is.na(out_gcst)]
    
    if (length(exp_gcst) == 0 || length(out_gcst) == 0) {
      mr_rv$status <- "Please select at least 1 Exposure and 1 Outcome (Tabs 1–2)."
      mr_rv$pairs <- character(0)
      mr_rv$results_by_pair <- list()
      mr_rv$snp_effects_by_pair <- list()
      return()
    }
    
    ao <- get_available_outcomes_cached("tab3")
    if (inherits(ao, "error")) {
      mr_rv$status <- paste0("OpenGWAS error (available_outcomes): ", conditionMessage(ao))
      return()
    }
    
    do_clump  <- isTRUE(input$do_clump)
    p1        <- input$p_thresh
    r2        <- input$clump_r2
    kb        <- input$clump_kb
    use_proxies <- isTRUE(input$use_proxies)
    proxy_rsq  <- input$proxy_rsq
    allow_pal  <- isTRUE(input$allow_palindromes)
    maf_thr    <- input$maf_thresh
    harm_action <- harmonise_action_value(input$harm_action)
    method_list <- input$mr_methods
    if (is.null(method_list) || length(method_list) == 0) method_list <- c("mr_ivw")
    
    pair_keys <- as.vector(outer(exp_gcst, out_gcst, FUN = function(e, o) paste0(e, "  →  ", o)))
    mr_rv$pairs <- pair_keys
    mr_rv$results_by_pair <- list()
    mr_rv$snp_effects_by_pair <- list()
    mr_rv$status <- paste0("Performing MR for ", length(pair_keys), " pair(s)...")
    
    withProgress(message = "Performing MR (OpenGWAS calls)", value = 0, {
      total <- length(pair_keys)
      idx <- 0
      
      for (e_gcst in exp_gcst) {
        for (o_gcst in out_gcst) {
          idx <- idx + 1
          incProgress(1 / total, detail = paste0("Pair ", idx, "/", total, " : ", e_gcst, " → ", o_gcst))
          key <- paste0(e_gcst, "  →  ", o_gcst)
          
          e_id <- resolve_opengwas_id(e_gcst, "tab3")
          o_id <- resolve_opengwas_id(o_gcst, "tab3")
          
          if (is.na(e_id) || is.na(o_id)) {
            mr_rv$results_by_pair[[key]] <- list(
              error = paste0("Could not resolve OpenGWAS ID. Exposure: ", e_gcst, " -> ", e_id,
                             " | Outcome: ", o_gcst, " -> ", o_id)
            )
            next
          }
          
          ex_dat <- gwas_call_with_token(
            "tab3",
            extract_instruments(
              outcomes = e_id,
              p1 = p1,
              clump = do_clump,
              r2 = r2,
              kb = kb
            )
          )
          if (inherits(ex_dat, "error") || is.null(ex_dat) || nrow(as.data.frame(ex_dat)) == 0) {
            mr_rv$results_by_pair[[key]] <- list(
              error = paste0("extract_instruments failed for ", e_id, ": ",
                             if (inherits(ex_dat, "error")) conditionMessage(ex_dat) else "No instruments returned.")
            )
            next
          }
          
          out_dat <- extract_outcome_chunked_cached(
            snps        = ex_dat$SNP,
            outcome_id = o_id,
            token      = "tab3",
            cache_dir = run_cache_dir,
            proxies    = use_proxies,
            rsq        = proxy_rsq,
            palindromes= allow_pal,
            maf_threshold = maf_thr
          )
          
          if (handle_opengwas_limit(out_dat, opengwas_reset_rv)) {
            
            mr_rv$status <- "⛔ OpenGWAS rate limit reached. Please wait for reset."
            
            mr_rv$results_by_pair[[key]] <- list(
              error = "OpenGWAS rate limit reached. Retry after allowance reset."
            )
            
            return()
          }
          
          # ❌ REAL errors only (not rate limit)
          if (inherits(out_dat, "error") ||
              is.null(out_dat) ||
              (is.data.frame(out_dat) && nrow(out_dat) == 0)) {
            
            out <- list(
              error = paste0(
                "extract_outcome_data failed: ",
                if (inherits(out_dat, "error"))
                  conditionMessage(out_dat)
                else
                  "No outcome rows returned."
              )
            )

            mr_rv$results_by_pair[[key]] <- out
            next
          }
          
          hdat <- tryCatch(harmonise_data(exposure_dat = ex_dat, outcome_dat = out_dat, action = harm_action),
                           error = function(e) e)
          if (inherits(hdat, "error") || is.null(hdat) || nrow(as.data.frame(hdat)) == 0) {
            mr_rv$results_by_pair[[key]] <- list(
              error = paste0("harmonise_data failed: ",
                             if (inherits(hdat, "error")) conditionMessage(hdat) else "No harmonised rows.")
            )
            next
          }
          
          mr_res <- tryCatch(mr(hdat, method_list = method_list), error = function(e) e)
          mr_or  <- tryCatch(generate_odds_ratios(mr_res), error = function(e) e)
          
          het    <- tryCatch(mr_heterogeneity(hdat, method_list = c("mr_egger_regression", "mr_ivw")),
                             error = function(e) e)
          pleio  <- tryCatch(mr_pleiotropy_test(hdat), error = function(e) e)
          dirres <- tryCatch(directionality_test(hdat), error = function(e) e)
          
          single <- tryCatch(mr_singlesnp(hdat), error = function(e) e)
          loo    <- tryCatch(mr_leaveoneout(hdat), error = function(e) e)
          
          p_scatter <- tryCatch(mr_scatter_plot(mr_res, hdat), error = function(e) e)
          p_forest  <- tryCatch(mr_forest_plot(single), error = function(e) e)
          p_loo     <- tryCatch(mr_leaveoneout_plot(loo), error = function(e) e)
          p_funnel  <- tryCatch(mr_funnel_plot(single), error = function(e) e)
          
          snp_eff <- build_snp_effects_table(hdat)
          
          mr_rv$results_by_pair[[key]] <- list(
            exposure_gcst = e_gcst, outcome_gcst = o_gcst,
            exposure_id = e_id, outcome_id = o_id,
            n_instruments = nrow(as.data.frame(ex_dat)),
            n_harmonised  = nrow(as.data.frame(hdat)),
            snps_inst = unique(as.character(ex_dat$SNP)),
            harmonised = hdat,
            mr_or = mr_or, heterogeneity = het, pleiotropy = pleio,
            direction = dirres,
            plot_scatter = p_scatter, plot_forest = p_forest,
            plot_loo = p_loo, plot_funnel = p_funnel
          )
          
          mr_rv$snp_effects_by_pair[[key]] <- snp_eff
          persist_pair_snp_effects(e_gcst, o_gcst, snp_eff)
          update_outcome_cache(o_gcst, snp_eff)
        }
      }
    })
    
    if (length(mr_rv$pairs) > 0) {
      updateSelectInput(
        session,
        "pair_pick",
        choices = mr_rv$pairs,
        selected = mr_rv$pairs[[1]]
      )
    }
    
    mr_rv$status <- "Finished  (Done. Pick a pair to view results.)"
    session$sendCustomMessage("playSound", list())
  })
  
  # UI select for pair (Tab 3)
  output$pair_picker_ui <- renderUI({
    
    if (is.null(mr_rv$results_by_pair) || length(mr_rv$results_by_pair) == 0) {
      return(
        tags$p(
          style = "opacity:0.7; font-style:italic;",
          "No MR results yet."
        )
      )
    }
    
    selectInput(
      "pair_pick",
      "Choose Exposure → Outcome pair",
      choices = names(mr_rv$results_by_pair)
    )
  })
  
  # --------------------------------
  # TAB 3: download handlers
  # --------------------------------
  
  output$dl_tab3_mr_xlsx <- downloadHandler(
    filename = function() "TAB3_MR_Results.xlsx",
    content = function(file) {
      r <- current_pair_res()
      req(!is.null(r), is.null(r$error))
      openxlsx::write.xlsx(as_df_or_msg(r$mr_or), file)
    }
  )
  
  output$dl_tab3_sens_xlsx <- downloadHandler(
    filename = function() "TAB3_Sensitivity_Tables.xlsx",
    content = function(file) {
      r <- current_pair_res()
      req(!is.null(r), is.null(r$error))
      
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "Heterogeneity")
      openxlsx::writeData(wb, "Heterogeneity", as_df_or_msg(r$heterogeneity))
      
      openxlsx::addWorksheet(wb, "Pleiotropy")
      openxlsx::writeData(wb, "Pleiotropy", as_df_or_msg(r$pleiotropy))
      
      openxlsx::addWorksheet(wb, "Direction")
      openxlsx::writeData(wb, "Direction", as_df_or_msg(r$direction))
      
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  output$dl_tab3_snp_xlsx <- downloadHandler(
    filename = function() "TAB3_SNP_Effects.xlsx",
    content = function(file) {
      df <- current_pair_snp_table()
      if (is.null(df) || nrow(df) == 0) {
        df <- data.frame(Message = "No SNP effects table.")
      }
      openxlsx::write.xlsx(df, file)
    }
  )
  
  # Helpers to fetch current results by pair
  current_pair_res <- reactive({
    req(input$pair_pick)
    mr_rv$results_by_pair[[input$pair_pick]]
  })
  # ===============================
  # TAB 3: has-results flag
  # ===============================
  output$tab3_has_results <- reactive({
    # No MR run yet
    if (is.null(mr_rv$results_by_pair)) return(FALSE)
    
    # No results stored yet
    if (length(mr_rv$results_by_pair) == 0) return(FALSE)
    
    # If pair not selected yet
    if (is.null(input$pair_pick)) return(FALSE)
    
    # If selected pair has no results
    if (!input$pair_pick %in% names(mr_rv$results_by_pair)) return(FALSE)
    
    r <- mr_rv$results_by_pair[[input$pair_pick]]
    
    # If error happened
    if (is.null(r) || !is.null(r$error)) return(FALSE)
    
    TRUE
  })
  outputOptions(output, "tab3_has_results", suspendWhenHidden = FALSE)
  current_pair_snp_table <- reactive({
    req(input$pair_pick)
    mr_rv$snp_effects_by_pair[[input$pair_pick]]
  })
  
  # ---- MR plots
  output$plot_scatter <- renderPlot({
    r <- current_pair_res()
    if (is.null(r) || !is.null(r$error)) return(invisible())
    if (inherits(r$plot_scatter, "error")) return(invisible())
    if (is.list(r$plot_scatter) && length(r$plot_scatter) >= 1)
      print(r$plot_scatter[[1]])
    else
      print(r$plot_scatter)
  })
  output$plot_forest <- renderPlot({
    r <- current_pair_res()
    if (is.null(r) || !is.null(r$error)) return(invisible())
    if (inherits(r$plot_forest, "error")) return(invisible())
    if (is.list(r$plot_forest) && length(r$plot_forest) >= 1)
      print(r$plot_forest[[1]])
    else
      print(r$plot_forest)
  })
  output$plot_loo <- renderPlot({
    r <- current_pair_res()
    if (is.null(r) || !is.null(r$error)) return(invisible())
    if (inherits(r$plot_loo, "error")) return(invisible())
    if (is.list(r$plot_loo) && length(r$plot_loo) >= 1)
      print(r$plot_loo[[1]])
    else
      print(r$plot_loo)
  })
  output$plot_funnel <- renderPlot({
    r <- current_pair_res()
    if (is.null(r) || !is.null(r$error)) return(invisible())
    if (inherits(r$plot_funnel, "error")) return(invisible())
    if (is.list(r$plot_funnel) && length(r$plot_funnel) >= 1)
      print(r$plot_funnel[[1]])
    else
      print(r$plot_funnel)
  })
  
  output$plot_tab7_diagram <- renderGrViz({
    
    req(tab7_rv$res_mvmr, tab7_rv$res_2smr)
    
    if ("Message" %in% names(tab7_rv$res_mvmr)) return(NULL)
    
    # -----------------------------
    # Extract GCST IDs
    # -----------------------------
    exp_gcst <- tab7_rv$res_mvmr$Exposure[1]
    med_gcst <- tab7_rv$res_mvmr$Mediator[1]
    out_gcst <- tab7_rv$res_mvmr$Outcome[1]
    
    # -----------------------------
    # Convert GCST → Trait Name
    # -----------------------------
    get_trait <- function(gcst) {
      x <- DATA$merged %>%
        dplyr::filter(`STUDY ACCESSION` == gcst) %>%
        dplyr::pull(MAPPED_TRAIT)
      if (length(x) == 0 || is.na(x[1])) return(gcst)
      x[1]
    }
    
    exp_trait <- get_trait(exp_gcst)
    med_trait <- get_trait(med_gcst)
    out_trait <- get_trait(out_gcst)
    
    exp_label <- paste0(exp_trait)
    med_label <- paste0(med_trait)
    out_label <- paste0(out_trait)
    # -----------------------------
    # Extract Effects
    # -----------------------------
    beta_total  <- as.numeric(tab7_rv$res_mvmr$Total_X_to_Y)
    beta_direct <- as.numeric(tab7_rv$res_mvmr$Direct_X_to_Y_adjM)
    
    beta_xm  <- as.numeric(tab7_rv$res_2smr$X_to_M)
    beta_my  <- as.numeric(tab7_rv$res_2smr$M_to_Y_adjX)
    
    p_direct <- as.numeric(tab7_rv$res_mvmr$P_DIFF)
    p_prod   <- as.numeric(tab7_rv$res_2smr$P_PROD)
    
    # Significance colors
    col_direct <- ifelse(p_direct < 0.05, "darkgreen", "firebrick")
    col_ind    <- ifelse(p_prod < 0.05, "darkgreen", "firebrick")
    
    DiagrammeR::grViz(paste0("
    digraph mediation {
      
      graph [layout = dot, rankdir = LR, bgcolor = 'white']
      
      node [shape = box, style = filled, fillcolor = '#f4f6f8',
            fontname = Helvetica, fontsize = 12]
      
      X [label = '", exp_label, "']
      M [label = '", med_label, "']
      Y [label = '", out_label, "']
      
      # X → M
      X -> M [
        label = 'X→M\\nβ = ", round(beta_xm,3), "',
        color = 'black',
        fontsize = 11
      ]
      
      # M → Y
      M -> Y [
        label = 'M→Y|X\\nβ = ", round(beta_my,3), "',
        color = 'black',
        fontsize = 11
      ]
      
      # Direct X → Y
      X -> Y [
        label = 'Direct\\nβ = ", round(beta_direct,3), "',
        color = '", col_direct, "',
        penwidth = 2,
        fontsize = 11
      ]
      
      # Indirect summary text
      Indirect [shape = note,
        label = 'Indirect (Product)\\nβ = ", 
                             round(beta_xm * beta_my,3),
                             "\\np = ", signif(p_prod,3), "']
      
      {rank = same; X; M; Y}
    }
  "))
  })
  
    # ---- MR results tables
  output$tbl_mr_results <- renderDT({
    
    r <- current_pair_res()
    
    if (is.null(r)) {
      return(datatable(
        data.frame(Message = "No results."),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    if (!is.null(r$error)) {
      return(datatable(
        data.frame(Message = r$error),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    df <- as_df_or_msg(r$mr_or)
    
    if (nrow(df) == 0 || "Message" %in% names(df)) {
      return(datatable(
        data.frame(Message = "No MR results."),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    # ---------------------------------------
    # Add trait labels
    # ---------------------------------------
    df$id.exposure <- get_trait_label(r$exposure_gcst, DATA$merged)
    df$id.outcome  <- get_trait_label(r$outcome_gcst,  DATA$merged)
    
    df$exposure <- NULL
    df$outcome  <- NULL
    
    # ---------------------------------------
    # Convert numeric
    # ---------------------------------------
    df$b  <- as.numeric(df$b)
    df$se <- as.numeric(df$se)
    
    # ---------------------------------------
    # Detect if outcome is binary
    # ---------------------------------------
    is_binary_outcome <- FALSE
    
    out_ids <- unique(stringr::str_trim(as.character(outcome_mod$selected_gcst())))
    out_ids <- out_ids[out_ids != "" & !is.na(out_ids)]
    
    if (length(out_ids) > 0) {
      row <- DATA$merged %>%
        dplyr::filter(`STUDY ACCESSION` == out_ids[1])
      
      if (nrow(row) > 0) {
        is_binary_outcome <- identical(as.character(row$`Data Type`[1]), "Binary")
      }
    }
    
    z <- 1.96
    
    # ---------------------------------------
    # Build Beta CI (always shown)
    # ---------------------------------------
    lci_beta <- df$b - z * df$se
    uci_beta <- df$b + z * df$se
    
    df$`Beta (95% CI)` <- sprintf(
      "%.4f (%.4f, %.4f)",
      df$b, lci_beta, uci_beta
    )
    
    # ---------------------------------------
    # Build OR CI (only if binary outcome)
    # ---------------------------------------
    if (is_binary_outcome) {
      
      or     <- exp(df$b)
      lci_or <- exp(lci_beta)
      uci_or <- exp(uci_beta)
      
      df$`Odds Ratio (95% CI)` <- sprintf(
        "%.4f (%.4f, %.4f)",
        or, lci_or, uci_or
      )
    }
    
    # ---------------------------------------
    # Add significance column
    # ---------------------------------------
    df <- add_significance_column(df)
    
    # ---------------------------------------
    # Format p-value (4 dp)
    # ---------------------------------------
    if ("pval" %in% names(df)) {
      df$pval <- sprintf("%.4f", as.numeric(df$pval))
    }
    
    # ---------------------------------------
    # Keep final display columns
    # (include Standard Error)
    # ---------------------------------------
    keep_cols <- c(
      "id.exposure",
      "id.outcome",
      "method",
      "nsnp",
      "Beta (95% CI)",
      "se"
    )
    
    if (is_binary_outcome) {
      keep_cols <- c(keep_cols, "Odds Ratio (95% CI)")
    }
    
    keep_cols <- c(keep_cols, "pval", "Significance")
    keep_cols <- intersect(keep_cols, names(df))
    
    df <- df[, keep_cols, drop = FALSE]
    
    # ---------------------------------------
    # Clean column names
    # ---------------------------------------
    rename_map <- c(
      id.exposure = "Exposure (Trait | GCST ID)",
      id.outcome  = "Outcome (Trait | GCST ID)",
      method      = "MR Method",
      nsnp        = "Number of SNPs",
      se          = "Standard Error",
      pval        = "P-value",
      Significance = "Significance"
    )
    
    old_names <- names(df)
    new_names <- ifelse(old_names %in% names(rename_map),
                        rename_map[old_names],
                        old_names)
    names(df) <- new_names
    
    # Round standard error to 4 dp (optional but nice)
    if ("Standard Error" %in% names(df)) {
      df$`Standard Error` <- round(as.numeric(df$`Standard Error`), 4)
    }
    
    datatable(
      df,
      rownames = FALSE,
      options = list(
        pageLength = 10,
        autoWidth = TRUE,
        scrollX = TRUE,
        dom = "tip"
      )
    )
  })
  
  output$tbl_het <- renderDT({
    
    r <- current_pair_res()
    if (is.null(r) || !is.null(r$error)) {
      return(datatable(
        data.frame(Message="No results."),
        options=list(dom="t"),
        rownames=FALSE
      ))
    }
    
    df <- as_df_or_msg(r$heterogeneity)
    
    # Add trait labels
    df$id.exposure <- get_trait_label(r$exposure_gcst, DATA$merged)
    df$id.outcome  <- get_trait_label(r$outcome_gcst,  DATA$merged)
    
    df$exposure <- NULL
    df$outcome  <- NULL
    
    df <- add_significance_column(df)
    
    # Keep only important columns
    keep_cols <- intersect(c(
      "id.exposure","id.outcome",
      "method","Q","Q_df","Q_pval","Significance"
    ), names(df))
    
    df <- df[, keep_cols, drop=FALSE]
    
    # Rename columns (clean like MR table)
    rename_map <- c(
      id.exposure = "Exposure (Trait | GCST ID)",
      id.outcome  = "Outcome (Trait | GCST ID)",
      method      = "MR Method",
      Q           = "Cochran’s Q Statistic",
      Q_df        = "Degrees of Freedom",
      Q_pval      = "Q-test P-value",
      Significance = "Heterogeneity Significance"
    )
    
    colnames(df) <- rename_map[colnames(df)]
    
    df <- round_df_display(df,4)
    
    datatable(df, rownames=FALSE,
              options=list(pageLength=10, scrollX=TRUE, dom="tip"))
  })
  output$tbl_pleio <- renderDT({
    
    r <- current_pair_res()
    if (is.null(r) || !is.null(r$error)) {
      return(datatable(
        data.frame(Message="No results."),
        options=list(dom="t"),
        rownames=FALSE
      ))
    }
    
    df <- as_df_or_msg(r$pleiotropy)
    
    df$id.exposure <- get_trait_label(r$exposure_gcst, DATA$merged)
    df$id.outcome  <- get_trait_label(r$outcome_gcst,  DATA$merged)
    
    df$exposure <- NULL
    df$outcome  <- NULL
    
    df <- add_significance_column(df)
    
    keep_cols <- intersect(c(
      "id.exposure","id.outcome",
      "egger_intercept","se","pval","Significance"
    ), names(df))
    
    df <- df[, keep_cols, drop=FALSE]
    
    rename_map <- c(
      id.exposure = "Exposure (Trait | GCST ID)",
      id.outcome  = "Outcome (Trait | GCST ID)",
      egger_intercept = "MR-Egger Intercept",
      se = "Standard Error",
      pval = "P-value",
      Significance = "Pleiotropy Significance"
    )
    
    colnames(df) <- rename_map[colnames(df)]
    
    df <- round_df_display(df,4)
    
    datatable(df, rownames=FALSE,
              options=list(pageLength=10, scrollX=TRUE, dom="tip"))
  })
  output$tbl_dir <- renderDT({
    
    r <- current_pair_res()
    if (is.null(r) || !is.null(r$error)) {
      return(datatable(
        data.frame(Message="No results."),
        options=list(dom="t"),
        rownames=FALSE
      ))
    }
    
    df <- as_df_or_msg(r$direction)
    
    df$id.exposure <- get_trait_label(r$exposure_gcst, DATA$merged)
    df$id.outcome  <- get_trait_label(r$outcome_gcst,  DATA$merged)
    
    df$exposure <- NULL
    df$outcome  <- NULL
    
    df <- add_significance_column(df)
    
    keep_cols <- intersect(c(
      "id.exposure","id.outcome",
      "snp_r2.exposure","snp_r2.outcome",
      "correct_causal_direction",
      "steiger_pval","Significance"
    ), names(df))
    
    df <- df[, keep_cols, drop=FALSE]
    
    rename_map <- c(
      id.exposure = "Exposure (Trait | GCST ID)",
      id.outcome  = "Outcome (Trait | GCST ID)",
      snp_r2.exposure = "R² (Exposure)",
      snp_r2.outcome  = "R² (Outcome)",
      correct_causal_direction = "Direction Supported?",
      steiger_pval = "Steiger P-value",
      Significance = "Direction Significance"
    )
    
    colnames(df) <- rename_map[colnames(df)]
    
    df <- round_df_display(df,4)
    
    datatable(df, rownames=FALSE,
              options=list(pageLength=10, scrollX=TRUE, dom="tip"))
  })
  
  # Summarise Egger intercept for display
  output$pleio_summary <- renderText({
    r <- current_pair_res()
    if (is.null(r) || !is.null(r$error)) return("")
    df <- as_df_or_msg(r$pleiotropy)
    if (is.null(df) || nrow(df) == 0) return("")
    intercept <- NA; se <- NA; pval <- NA
    if ("egger_intercept" %in% names(df)) {
      intercept <- df$egger_intercept[1]
      if ("se" %in% names(df)) se <- df$se[1]
      if ("pval" %in% names(df)) pval <- df$pval[1]
    } else if ("intercept" %in% names(df)) {
      intercept <- df$intercept[1]
      if ("se" %in% names(df)) se <- df$se[1]
      if ("pval" %in% names(df)) pval <- df$pval[1]
    }
    intercept_txt <- ifelse(is.na(intercept), "Unknown", format(as.numeric(intercept), digits = 4))
    se_txt <- ifelse(is.na(se), "Unknown", format(as.numeric(se), digits = 4))
    pval_txt <- ifelse(is.na(pval), "Unknown", format(as.numeric(pval), digits = 4))
    paste0("Egger regression intercept: ", intercept_txt, "\n",
           "Standard error: ", se_txt, "\n",
           "Directionality p-value: ", pval_txt)
  })
  
  # ==========================================================
  # TAB 5A: Full MR engine for Exposure 2 -> Outcome (Tab3-style)
  # Uses token alias "tab5"
  # ==========================================================
  
  tab5_mr_rv <- reactiveValues(
    status = "Pick Exposure 2 and Outcome, then click Run Tab 5 MR.",
    res = NULL,
    snp_eff = NULL,
    singlesnp = NULL,
    loo = NULL,
    plots = list()
  )
  output$tab5_mr_status <- renderText(tab5_mr_rv$status)
  # ---- UI flag: does Tab 5 have results? (used to hide empty plots)
  output$tab5_has_results <- reactive({
    
    !is.null(tab5_mr_rv$res) &&
      is.null(tab5_mr_rv$res$error)
    
  })
  outputOptions(output, "tab5_has_results", suspendWhenHidden = FALSE)
  # --------------------------------------------------------------------
  # Tab 3 downloads (MR results, sensitivity, plots, SNP table) 
  # definitions omitted here due to length but kept unchanged from single-file version
  # (These include dl_mr_xlsx, dl_sens_xlsx, dl_plot_scatter, dl_plot_forest,
  #  dl_plot_loo, dl_plot_funnel, dl_snp_table_xlsx)
  # --------------------------------------------------------------------
  
  output$dl_tab5_snp_table_xlsx <- downloadHandler(
    filename = function() "TAB5_SNP_Effects.xlsx",
    content = function(file) {
      df <- tab5_mr_rv$snp_eff
      if (is.null(df) || nrow(df)==0) df <- data.frame(Message="No SNP effects table.")
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "SNP_effects")
      openxlsx::writeData(wb, "SNP_effects", df)
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  output$dl_tab5_mr_xlsx <- downloadHandler(
    filename = function() "TAB5_MR_Results.xlsx",
    content = function(file) {
      r <- tab5_mr_rv$res
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "MR_results")
      openxlsx::writeData(wb, "MR_results", as_df_or_msg(r$mr_or))
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  output$dl_tab5_sens_xlsx <- downloadHandler(
    filename = function() "TAB5_Sensitivity_Tables.xlsx",
    content = function(file) {
      r <- tab5_mr_rv$res
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "Heterogeneity")
      openxlsx::writeData(wb, "Heterogeneity", as_df_or_msg(r$heterogeneity))
      openxlsx::addWorksheet(wb, "Pleiotropy")
      openxlsx::writeData(wb, "Pleiotropy", as_df_or_msg(r$pleiotropy))
      openxlsx::addWorksheet(wb, "Directionality")
      openxlsx::writeData(wb, "Directionality", as_df_or_msg(r$direction))
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  # ---- Plot PDFs
  .output_plot_pdf <- function(plot_obj, file) {
    grDevices::pdf(file, width = 11, height = 8.5)
    on.exit(grDevices::dev.off(), add=TRUE)
    if (is.list(plot_obj) && length(plot_obj) >= 1) print(plot_obj[[1]]) else print(plot_obj)
  }
  # -------------------------------
  # TAB 3 plot downloads
  # -------------------------------
  
  output$dl_plot_scatter <- downloadHandler(
    filename = function() "TAB3_Method_Comparison.pdf",
    content = function(file) {
      r <- current_pair_res()
      req(!is.null(r), is.null(r$error))
      .output_plot_pdf(r$plot_scatter, file)
    }
  )
  
  output$dl_plot_forest <- downloadHandler(
    filename = function() "TAB3_Single_SNP_Analysis.pdf",
    content = function(file) {
      r <- current_pair_res()
      req(!is.null(r), is.null(r$error))
      .output_plot_pdf(r$plot_forest, file)
    }
  )
  
  output$dl_plot_loo <- downloadHandler(
    filename = function() "TAB3_Leave_One_Out.pdf",
    content = function(file) {
      r <- current_pair_res()
      req(!is.null(r), is.null(r$error))
      .output_plot_pdf(r$plot_loo, file)
    }
  )
  
  output$dl_plot_funnel <- downloadHandler(
    filename = function() "TAB3_Funnel_Plot.pdf",
    content = function(file) {
      r <- current_pair_res()
      req(!is.null(r), is.null(r$error))
      .output_plot_pdf(r$plot_funnel, file)
    }
  )
  
  output$dl_tab5_plot_scatter <- downloadHandler(
    filename = function() "TAB5_Method_Comparison.pdf",
    content = function(file) {
      p <- tab5_mr_rv$plots$scatter
      if (is.null(p) || inherits(p,"error")) { grDevices::pdf(file); grDevices::plot.new(); text(0.5,0.5,"No plot"); grDevices::dev.off(); return() }
      .output_plot_pdf(p, file)
    }
  )
  
  output$dl_tab5_plot_forest <- downloadHandler(
    filename = function() "TAB5_Single_SNP_Forest.pdf",
    content = function(file) {
      p <- tab5_mr_rv$plots$forest
      if (is.null(p) || inherits(p,"error")) { grDevices::pdf(file); grDevices::plot.new(); text(0.5,0.5,"No plot"); grDevices::dev.off(); return() }
      .output_plot_pdf(p, file)
    }
  )
  
  output$dl_tab5_plot_loo <- downloadHandler(
    filename = function() "TAB5_Leave_One_Out.pdf",
    content = function(file) {
      p <- tab5_mr_rv$plots$loo
      if (is.null(p) || inherits(p,"error")) { grDevices::pdf(file); grDevices::plot.new(); text(0.5,0.5,"No plot"); grDevices::dev.off(); return() }
      .output_plot_pdf(p, file)
    }
  )
  
  output$dl_tab5_plot_funnel <- downloadHandler(
    filename = function() "TAB5_Funnel.pdf",
    content = function(file) {
      p <- tab5_mr_rv$plots$funnel
      if (is.null(p) || inherits(p,"error")) { grDevices::pdf(file); grDevices::plot.new(); text(0.5,0.5,"No plot"); grDevices::dev.off(); return() }
      .output_plot_pdf(p, file)
    }
  )
  # --------------------------------------------------------------------
  # TAB 4: REPEATED SNP
  # Manage repeated SNP table and counts by categories
  # --------------------------------------------------------------------
  rep_filter_rv <- reactiveValues(cat = NULL, sub = NULL)
  observeEvent(input$clear_rep_filters, {
    rep_filter_rv$cat <- NULL
    rep_filter_rv$sub <- NULL
  })
  
  repeated_tbl_raw <- reactive({
    r <- current_pair_res()
    snp_eff <- current_pair_snp_table()
    
    if (is.null(r) || !is.null(r$error)) {
      return(data.frame(Message = "Run Tab 3 MR and pick a valid pair first (no errors)."))
    }
    if (is.null(snp_eff) || nrow(snp_eff) == 0) {
      return(data.frame(Message = "No SNP effects table for current pair. Run MR and ensure harmonised SNPs exist."))
    }
    if (is.null(DATA$snp) || nrow(DATA$snp) == 0 || !all(c("GCST", "SNP") %in% names(DATA$snp))) {
      return(data.frame(Message = "snp_level_ALL_instruments.csv missing or must contain GCST and SNP columns."))
    }
    
    # ---- VALID SNP COUNTS (same logic as Tab 1 / Tab 2) ----
    valid_snp_counts <- DATA$snp %>%
      mutate(
        GCST = str_trim(as.character(GCST)),
        SNP  = str_trim(as.character(SNP))
      ) %>%
      filter(GCST != "", SNP != "", !is.na(GCST), !is.na(SNP)) %>%
      group_by(GCST) %>%
      summarise(
        Valid_SNPs = n_distinct(SNP),
        .groups = "drop"
      )
    snps_use <- unique(str_trim(as.character(snp_eff$SNP)))
    snps_use <- snps_use[snps_use != "" & !is.na(snps_use)]
    if (length(snps_use) == 0) return(data.frame(Message = "Current SNP table has no SNP IDs."))
    
    # Determine IDs to exclude from candidate set (exposure and outcome) and
    # reference IDs for independence testing.  Exclude both exposure and outcome GCSTs
    # from the candidate search, but only use exposure GCSTs as the reference for
    # independence checks (updated rule).
    exclude_ids <- unique(c(r$exposure_gcst, r$outcome_gcst))
    exclude_ids <- exclude_ids[exclude_ids != "" & !is.na(exclude_ids)]
    
    cand <- DATA$snp %>%
      filter(SNP %in% snps_use) %>%
      filter(!(GCST %in% exclude_ids)) %>%
      group_by(GCST) %>%
      summarise(
        overlap_snps = n_distinct(SNP),
        overlap_snp_examples = paste(head(sort(unique(SNP)), 10), collapse = ", "),
        .groups = "drop"
      ) %>%
      arrange(desc(overlap_snps), GCST)
    
    cand <- cand %>%
      left_join(valid_snp_counts, by = "GCST")
    
    
    if (nrow(cand) == 0) return(data.frame(Message = "No repeated-SNP studies found for the current MR SNP set."))
    
    # Independence is evaluated using only exposure GCSTs as reference
    ref_ids <- unique(r$exposure_gcst)
    ref_ids <- ref_ids[ref_ids != "" & !is.na(ref_ids)]
    cand$Independence <- vapply(cand$GCST, function(cid) {
      ok <- is_independent_from_all(cid, ref_ids, DATA$merged)
      if (is.na(ok)) "Unknown" else if (ok) "Independent" else "Not independent"
    }, character(1))
    
    # Compute cohort overlap status between candidate GCSTs and the reference set (exposure + outcome).
    # We build a set of cohort names from the reference IDs using parse_cohort() and then
    # check whether each candidate's cohort shares any of these names.
    # Build reference cohorts from both exposure and outcome GCSTs (exclude_ids)
    ref_cohorts <- unique(unlist(lapply(exclude_ids, function(gid) {
      rows <- which(DATA$merged$`STUDY ACCESSION` == gid)
      if (length(rows) == 0) return(character(0))
      co_strs <- DATA$merged$COHORT[rows]
      co_strs <- co_strs[!is.na(co_strs) & co_strs != ""]
      if (length(co_strs) == 0) return(character(0))
      unique(unlist(lapply(co_strs, parse_cohort)))
    })))
    if (length(ref_cohorts) == 0) {
      cand$`Cohort Overlap` <- rep("UNKNOWN", nrow(cand))
    } else {
      cand$`Cohort Overlap` <- vapply(cand$GCST, function(cid) {
        rows <- which(DATA$merged$`STUDY ACCESSION` == cid)
        if (length(rows) == 0) return("UNKNOWN")
        co_strs <- DATA$merged$COHORT[rows]
        co_strs <- co_strs[!is.na(co_strs) & co_strs != ""]
        if (length(co_strs) == 0) return("UNKNOWN")
        candidate_coh <- unique(unlist(lapply(co_strs, parse_cohort)))
        # If either reference or candidate contains 'MULTIPLE', treat as overlapping.
        if ("MULTIPLE" %in% candidate_coh || "MULTIPLE" %in% ref_cohorts) {
          return("YES")
        }
        if (length(intersect(ref_cohorts, candidate_coh)) > 0) "YES" else "NO"
      }, character(1))
    }
    
    if (!is.null(DATA$groups) && nrow(DATA$groups) > 0) {
      cand <- cand %>% left_join(DATA$groups, by = c("GCST" = "GCST"))
    } else {
      cand$TOP_CATEGORY_AUTO <- NA_character_
      cand$SUB_GROUP_AUTO    <- NA_character_
    }
    
    # Apply independence filtering: when the first checkbox is selected, retain
    # only "Not independent" (same‑ancestry) candidates.  When the second
    # checkbox is also selected, retain only those with no cohort overlap ("NO").
    if (isTRUE(input$rep_only_indep)) {
      cand <- cand %>% dplyr::filter(Independence == "Not independent")
      if (!is.null(input$rep_non_overlap) && isTRUE(input$rep_non_overlap)) {
        cand <- cand %>% dplyr::filter(`Cohort Overlap` == "NO")
      }
    }
    if (!is.null(rep_filter_rv$cat)) {
      cand <- cand %>%
        filter(ifelse(is.na(TOP_CATEGORY_AUTO) | TOP_CATEGORY_AUTO == "", "Unknown",
                      TOP_CATEGORY_AUTO) == rep_filter_rv$cat)
    }
    if (!is.null(rep_filter_rv$sub)) {
      cand <- cand %>%
        filter(ifelse(is.na(SUB_GROUP_AUTO) | SUB_GROUP_AUTO == "", "Unknown",
                      SUB_GROUP_AUTO) == rep_filter_rv$sub)
    }
    
    # ------------------------------
    # Add mapped trait name here
    # ------------------------------
    
    if (!is.null(DATA$merged) &&
        "STUDY ACCESSION" %in% names(DATA$merged)) {
      
      trait_lookup <- DATA$merged %>%
        dplyr::select(
          GCST = `STUDY ACCESSION`,
          MAPPED_TRAIT
        ) %>%
        distinct()
      
      cand <- cand %>%
        left_join(trait_lookup, by = "GCST")
    }
    
    # Final table including cohort overlap
    cand %>%
      transmute(
        GCST,
        MAPPED_TRAIT,
        overlap_snps,
        overlap_snp_examples,
        Independence,
        `Cohort Overlap`,
        TOP_CATEGORY_AUTO,
        SUB_GROUP_AUTO,
        Valid_SNPs
      )
    
  })
  
  # Render repeated SNP table
  output$tbl_repeated <- renderDT({
    df <- repeated_tbl_raw()
    if ("Message" %in% names(df)) {
      return(datatable(df, rownames = FALSE, options = list(dom = "t")))
    }
    
    # Correct join using GCST
    df <- df %>%
      dplyr::left_join(
        DATA$merged %>%
          dplyr::select(
            GCST = `STUDY ACCESSION`,
            `Data Type`
          ) %>%
          distinct(),
        by = "GCST"
      )
    
    # Replace NA Data Type safely
    if (!"Data Type" %in% names(df)) {
      df$`Data Type` <- NA_character_
    }
    
    colnames(df) <- c(
      "Study ID (GCST)",
      "Trait Name",
      "Overlapping SNP Count",
      "Example Overlapping SNPs",
      "Independence Status",
      "Cohort Overlap",
      "Top Category",
      "Sub Category",
      "Total Valid SNPs in Study",
      "Data Type"
    )
    
    # ---- DISPLAY CHANGE ONLY ----
    df$`Independence Status` <- dplyr::recode(
      df$`Independence Status`,
      "Independent" = "Different",
      "Not independent" = "Same"
    )
    
    colnames(df)[colnames(df) == "Independence Status"] <- "Ancestry Status"
    # -----------------------------
    df <- df %>%
      dplyr::relocate(`Data Type`, .after = `Trait Name`)
    datatable(df,
              rownames = FALSE,
              selection = "multiple",
              options = list(pageLength = 12, scrollX = TRUE, dom = "tip")) %>%
      formatStyle(
        "Ancestry Status",
        fontWeight = "bold",
        color = styleEqual(
          c("Different", "Same", "Unknown"),
          c("#34d399", "#fb7185", "#fbbf24")
        )
      )
  })
  
  # Category & subcategory counts
  selected_repeated_gcst <- reactive({
    df <- repeated_tbl_raw()
    sel <- input$tbl_repeated_rows_selected
    if ("Message" %in% names(df)) return(character(0))
    if (is.null(sel) || length(sel) == 0) return(character(0))
    ids <- df$GCST[sel]
    ids <- unique(str_trim(as.character(ids)))
    ids[ids != "" & !is.na(ids)]
  })
  
  observeEvent(selected_repeated_gcst(), {
    
    ids <- selected_repeated_gcst()
    if (length(ids) == 0) return()
    
    check_id <- ids[1]
    
    row <- DATA$merged %>%
      dplyr::filter(`STUDY ACCESSION` == check_id)
    
    if (nrow(row) == 0) return()
    
    study_type <- row$`Data Type`[1]
    
    if (!is.na(study_type) && study_type == "Binary") {
      
      showNotification(
        paste0(
          "⚠ Exposure 2 / Mediator ", check_id, " (",
          row$MAPPED_TRAIT[1],
          ") in Related SNP table is binary. ",
          "X → M estimates are in log-odds scale. ",
          "Use exp(beta) for Odds Ratio interpretation."
        ),
        type = "warning",
        duration = 8
      )
      
    }
    
  })
  
  repeated_tbl_for_counts <- reactive({
    df <- repeated_tbl_raw()
    if ("Message" %in% names(df)) return(df)
    sel_ids <- selected_repeated_gcst()
    if (length(sel_ids) == 0) return(df)
    df %>% filter(GCST %in% sel_ids)
  })
  
  output$tbl_cat_counts <- renderDT({
    df <- repeated_tbl_for_counts()
    if ("Message" %in% names(df)) {
      return(datatable(data.frame(Message = "No summary (need repeated-SNP results)."),
                       rownames = FALSE, options = list(dom = "t")))
    }
    x <- df %>%
      mutate(TOP_CATEGORY_AUTO = ifelse(is.na(TOP_CATEGORY_AUTO) | TOP_CATEGORY_AUTO == "", "Unknown", TOP_CATEGORY_AUTO)) %>%
      group_by(TOP_CATEGORY_AUTO) %>%
      summarise(STUDY_COUNT = n_distinct(GCST), .groups = "drop") %>%
      arrange(desc(STUDY_COUNT), TOP_CATEGORY_AUTO)
    
    colnames(x) <- c(
      "Top Category",
      "Number of Studies"
    )
    
    datatable(x,
              rownames = FALSE,
              selection = "single",
              options = list(pageLength = 12, scrollX = TRUE, dom = "tip"))
    
  })
  
  output$tbl_subcat_counts <- renderDT({
    df <- repeated_tbl_for_counts()
    if ("Message" %in% names(df)) {
      return(datatable(data.frame(Message = "No summary (need repeated-SNP results)."),
                       rownames = FALSE, options = list(dom = "t")))
    }
    x <- df %>%
      mutate(SUB_GROUP_AUTO = ifelse(is.na(SUB_GROUP_AUTO) | SUB_GROUP_AUTO == "", "Unknown", SUB_GROUP_AUTO)) %>%
      group_by(SUB_GROUP_AUTO) %>%
      summarise(STUDY_COUNT = n_distinct(GCST), .groups = "drop") %>%
      arrange(desc(STUDY_COUNT), SUB_GROUP_AUTO)
    
    colnames(x) <- c(
      "Sub Category",
      "Number of Studies"
    )
    
    datatable(x,
              rownames = FALSE,
              selection = "single",
              options = list(pageLength = 12, scrollX = TRUE, dom = "tip"))
  })
  
  observeEvent(input$tbl_cat_counts_rows_selected, {
    df_counts <- isolate(repeated_tbl_for_counts())
    if ("Message" %in% names(df_counts)) return()
    x <- df_counts %>%
      mutate(TOP_CATEGORY_AUTO = ifelse(is.na(TOP_CATEGORY_AUTO) | TOP_CATEGORY_AUTO == "", "Unknown", TOP_CATEGORY_AUTO)) %>%
      group_by(TOP_CATEGORY_AUTO) %>%
      summarise(STUDY_COUNT = n_distinct(GCST), .groups = "drop") %>%
      arrange(desc(STUDY_COUNT), TOP_CATEGORY_AUTO)
    sel <- input$tbl_cat_counts_rows_selected
    if (is.null(sel) || length(sel) == 0) return()
    rep_filter_rv$cat <- x$TOP_CATEGORY_AUTO[sel[1]]
  })
  
  observeEvent(input$tbl_subcat_counts_rows_selected, {
    df_counts <- isolate(repeated_tbl_for_counts())
    if ("Message" %in% names(df_counts)) return()
    x <- df_counts %>%
      mutate(SUB_GROUP_AUTO = ifelse(is.na(SUB_GROUP_AUTO) | SUB_GROUP_AUTO == "", "Unknown", SUB_GROUP_AUTO)) %>%
      group_by(SUB_GROUP_AUTO) %>%
      summarise(STUDY_COUNT = n_distinct(GCST), .groups = "drop") %>%
      arrange(desc(STUDY_COUNT), SUB_GROUP_AUTO)
    sel <- input$tbl_subcat_counts_rows_selected
    if (is.null(sel) || length(sel) == 0) return()
    rep_filter_rv$sub <- x$SUB_GROUP_AUTO[sel[1]]
  })
  
  # --------------------------------------------------------------------
  # TAB 5: Exposure 2 selection (build SNP table)
  # --------------------------------------------------------------------
  
  output$exp1_auto_box <- renderText({
    exp1 <- unique(str_trim(as.character(exposure_mod$selected_gcst())))
    exp1 <- exp1[exp1 != "" & !is.na(exp1)]
    if (length(exp1) == 0) return("No Exposure 1 selected in Tab 1 yet.")
    paste(exp1, collapse = ", ")
  })
  
  exp2_candidates <- reactive({
    df <- repeated_tbl_raw()
    if ("Message" %in% names(df)) return(character(0))
    ids <- unique(df$GCST)
    ids <- ids[!is.na(ids) & ids != ""]
    sort(ids)
  })
  
  rv_exp2 <- reactiveValues(selected = character(0))
  
  output$exp2_pick_ui <- renderUI({
    
    ids <- exp2_candidates()
    if (length(ids) == 0) {
      return(div(class="note", "No candidates yet. Use Tab 3 + Tab 4 first."))
    }
    
    selectizeInput(
      "exp2_pick",
      label = NULL,
      choices = ids,
      selected = rv_exp2$selected,  
      multiple = FALSE,
      options = list(placeholder = "Exposure 2 from Tab 4")
    )
  })
  
  observeEvent(input$tbl_repeated_rows_selected, {
    
    df <- isolate(repeated_tbl_raw())
    sel <- input$tbl_repeated_rows_selected
    
    if (is.null(sel) || length(sel) == 0) return()
    if ("Message" %in% names(df)) return()
    
    pick <- str_trim(as.character(df$GCST[sel[1]]))
    if (is.na(pick) || pick == "") return()
    
    # STORE ONLY (like Tab 1 & Tab 2)
    rv_exp2$selected <- pick
    
  }, ignoreInit = TRUE)
  
  observeEvent(input$exp2_pick, {
    x <- str_trim(as.character(input$exp2_pick))
    if (is.na(x) || x == "") rv_exp2$selected <- character(0) else rv_exp2$selected <- x
  })
  
  output$exp2_pick_box <- renderText({
    if (length(rv_exp2$selected) == 0) return("No Exposure 2 selected yet.")
    rv_exp2$selected
  })
  
  exposure2_selected <- reactive({
    x <- rv_exp2$selected
    x <- unique(str_trim(as.character(x)))
    x <- x[x != "" & !is.na(x)]
    x
  })
  
  # ==========================================
  # Mediator selection (NEW unified structure)
  # ==========================================
  mediator_gcst <- reactive({
    req(input$mediator_source)
    
    if (input$mediator_source == "from_tabA") {
      x <- exposure2_selected()[1]
      return(x)
    }
    
    # new_search
    req(tab7_search)
    x <- tab7_search$selected_gcst()[1]
    x
  })
  
  # --------------------------------------------------------------------
  # Cohort presence flags
  #
  # Define separate reactive flags to indicate when cohort information is
  # available on the selected exposure(s) alone (has_cohort_exp) and when
  # cohort information exists for both the selected exposure(s) and
  # outcome(s) (has_cohort_union).  The former controls the visibility of
  # the non‑overlap checkbox in the outcome tab, while the latter controls
  # it for the repeated‑SNP and mediator tabs.  We also provide a
  # backwards-compatible alias has_cohort_pair (same as has_cohort_exp).
  # -----------------------------------------
  # Cohort presence logic (PURE reactive)
  # -----------------------------------------
  
  has_cohort_exp <- reactive({
    exp_gcst <- unique(stringr::str_trim(as.character(exposure_mod$selected_gcst())))
    exp_gcst <- exp_gcst[exp_gcst != "" & !is.na(exp_gcst)]
    if (length(exp_gcst) == 0) return(FALSE)
    
    any(sapply(exp_gcst, function(gid) {
      rows <- which(DATA$merged$`STUDY ACCESSION` == gid)
      if (length(rows) == 0) return(FALSE)
      co <- DATA$merged$COHORT[rows]
      any(!is.na(co) & trimws(as.character(co)) != "")
    }))
  })
  
  has_cohort_union <- reactive({
    
    # Exposure check
    exp_gcst <- unique(stringr::str_trim(as.character(exposure_mod$selected_gcst())))
    exp_gcst <- exp_gcst[exp_gcst != "" & !is.na(exp_gcst)]
    if (length(exp_gcst) == 0) return(FALSE)
    
    exp_has <- any(sapply(exp_gcst, function(gid) {
      rows <- which(DATA$merged$`STUDY ACCESSION` == gid)
      if (length(rows) == 0) return(FALSE)
      co <- DATA$merged$COHORT[rows]
      any(!is.na(co) & trimws(as.character(co)) != "")
    }))
    
    if (!exp_has) return(FALSE)
    
    # Outcome check
    out_gcst <- unique(stringr::str_trim(as.character(outcome_mod$selected_gcst())))
    out_gcst <- out_gcst[out_gcst != "" & !is.na(out_gcst)]
    if (length(out_gcst) == 0) return(FALSE)
    
    any(sapply(out_gcst, function(gid) {
      rows <- which(DATA$merged$`STUDY ACCESSION` == gid)
      if (length(rows) == 0) return(FALSE)
      co <- DATA$merged$COHORT[rows]
      any(!is.na(co) & trimws(as.character(co)) != "")
    }))
  })
  
  has_cohort_pair <- reactive({
    has_cohort_exp()
  })
  
  output$has_cohort_exp <- renderText({
    has_cohort_exp()
  })
  
  output$has_cohort_union <- renderText({
    has_cohort_union()
  })
  
  output$has_cohort_pair <- renderText({
    has_cohort_pair()
  })
  # ensure reactivity persists
  outputOptions(output, "has_cohort_exp", suspendWhenHidden = FALSE)
  outputOptions(output, "has_cohort_union", suspendWhenHidden = FALSE)
  outputOptions(output, "has_cohort_pair", suspendWhenHidden = FALSE)
  
  # --------------------------------------------------------------------
  # Warning notifications for missing cohort and 'Other' populations
  #
  # Helper: show warnings based on a set of GCST IDs and the role
  warn_for_gcst <- function(ids, role) {
    ids <- unique(stringr::str_trim(as.character(ids)))
    ids <- ids[ids != "" & !is.na(ids)]
    if (length(ids) == 0) return(invisible(NULL))
    rows <- DATA$merged$`STUDY ACCESSION` %in% ids
    # Check populations labelled 'Other'
    pops <- unique(DATA$merged$EXTRACTED_TERMS[rows])
    pops <- pops[!is.na(pops) & pops != ""]
    if (any(tolower(pops) == "other")) {
      shiny::showNotification(
        paste0("⚠️ ", role, " population is labelled 'Other'. Please consult the GWAS catalog for more details."),
        type = "warning",
        duration = 10
      )
    }
    # Check missing cohort
    missing_co <- any(sapply(ids, function(gid) {
      ridx <- which(DATA$merged$`STUDY ACCESSION` == gid)
      if (length(ridx) == 0) return(TRUE)
      all(is.na(DATA$merged$COHORT[ridx]) | trimws(as.character(DATA$merged$COHORT[ridx])) == "")
    }))
    if (missing_co) {
      msg <- switch(
        role,
        "Exposure"   = "⚠️ COHORT was not there for the selected exposure(s). Sample overlap can occur and the non‑overlap filter is disabled.",
        "Outcome"    = "⚠️ COHORT was not there for the selected outcome(s). Sample overlap may occur and the non‑overlap filter is disabled.",
        "Exposure 2" = "⚠️ COHORT was not there for the selected Exposure 2. Sample overlap may occur and the non‑overlap filter is disabled.",
        "Mediator"   = "⚠️ COHORT was not there for the selected mediator. Sample overlap may occur and the non‑overlap filter is disabled.",
        paste0("⚠️ COHORT was not there for the selected ", role, ". Sample overlap may occur and the non‑overlap filter is disabled.")
      )
      shiny::showNotification(
        msg,
        type = "warning",
        duration = 10
      )
    }
    invisible(NULL)
  }
  
  # Observe selected exposure, outcome, exposure2 and mediator to trigger warnings
  observeEvent(exposure_mod$selected_gcst(), {
    warn_for_gcst(exposure_mod$selected_gcst(), role = "Exposure")
  }, ignoreInit = TRUE)
  observeEvent(outcome_mod$selected_gcst(), {
    warn_for_gcst(outcome_mod$selected_gcst(), role = "Outcome")
  }, ignoreInit = TRUE)
  observeEvent(exposure2_selected(), {
    warn_for_gcst(exposure2_selected(), role = "Exposure 2")
  }, ignoreInit = TRUE)
  observeEvent(mediator_gcst(), {
    warn_for_gcst(mediator_gcst(), role = "Mediator")
  }, ignoreInit = TRUE)
  
  # Exposure 2 used in SNP table + MVMR should ALWAYS be the mediator
  exp2_for_mvmr <- reactive({
    x <- mediator_gcst()[1]
    x <- str_trim(as.character(x))
    if (is.na(x) || x == "") return(NA_character_)
    x
  })
  
  observeEvent(mediator_gcst(), {
    
    tab5_rv$df <- data.frame(
      Message = "Mediator changed. Rebuild SNP table.",
      stringsAsFactors = FALSE
    )
    
    mvmr_rv$res_ivw <- NULL
    mvmr_rv$res_strength <- NULL
    mvmr_rv$res_pleio <- NULL
    
    mvmr_rv$status <- "Mediator changed. Rebuild SNP table and rerun MVMR."
    
  }, ignoreInit = TRUE)
  # ============================================================
  # Tab 5 — reset cache + results when Exposure 2 or Outcome changes
  # ============================================================
  observeEvent(
    list(exposure2_selected(), outcome_mod$selected_gcst()),
    {
      
      tab5_mr_rv$res <- NULL
      tab5_mr_rv$snp_eff <- NULL
      tab5_mr_rv$status <- "Selection changed. Click Run MR (Exposure1→Mediator)"
    },
    ignoreInit = TRUE
  )
  # ============================================================
  # Tab 5: Exposure 1 → Mediator pair display
  # (UI box showing the MR study being run)
  # ============================================================
  
  tab5_pair_x <- reactive({
    req(exposure_mod$selected_gcst())
    unique(str_trim(as.character(exposure_mod$selected_gcst())))[1]
  })
  
  tab5_pair_m <- reactive({
    req(mediator_gcst())
    unique(str_trim(as.character(mediator_gcst())))[1]
  })
  
  output$tab5_pair_ui <- renderUI({
    
    req(tab5_pair_x(), tab5_pair_m())
    
    pair_choices <- paste0(
      tab5_pair_x(), "  →  ", tab5_pair_m()
    )
    
    selectInput(
      "tab5_pair",
      label = NULL,
      choices = pair_choices,
      selected = pair_choices,
      width = "100%"
    )
    
  })  
  
  # ==========================================================
  # TAB 5A: Full MR engine for Exposure 2 -> Outcome (Tab3-style)
  # ==========================================================
  
  observeEvent(input$run_mr_tab5, {
    run_cache_dir <- file.path(
      tempdir(),
      paste0("mr_tab5_", as.integer(Sys.time()), "_", sample(1e6,1))
    )
    
    dir.create(run_cache_dir, recursive = TRUE, showWarnings = FALSE)
    show_busy_modal("Performing Exposure 1 → Exposure 2 / Mediator MR...", "Running Exposure 1 → Exposure 2 / Mediator MR now…")
    on.exit({ removeModal() }, add = TRUE)
    
    exp2_gcst <- unique(str_trim(as.character(exposure_mod$selected_gcst())))[1]  # Exposure 1
    out_gcst  <- unique(str_trim(as.character(mediator_gcst())))[1]              # Mediator
    
    if (is.na(exp2_gcst) || exp2_gcst == "" || is.na(out_gcst) || out_gcst == "") {
      tab5_mr_rv$status <- "❌ Need Exposure 1(Tab 1) and  Mediator (Tab 4)."
      tab5_mr_rv$res <- NULL
      tab5_mr_rv$hdat <- NULL
      tab5_mr_rv$ivw_summary <- NULL
      tab5_mr_rv$meta <- NULL
      return()
    }
    
    e2_id <- resolve_opengwas_id(exp2_gcst, "tab5")
    o_id  <- resolve_opengwas_id(out_gcst,  "tab5")
    
    if (is.na(e2_id) || is.na(o_id)) {
      tab5_mr_rv$status <- paste0("❌ Could not resolve OpenGWAS IDs. Exp2=", e2_id, " Outcome=", o_id)
      tab5_mr_rv$res <- NULL
      tab5_mr_rv$hdat <- NULL
      tab5_mr_rv$ivw_summary <- NULL
      tab5_mr_rv$meta <- NULL
      return()
    }
    
    do_clump    <- isTRUE(input$tab5_do_clump)
    p1          <- input$tab5_p_thresh
    r2          <- input$tab5_clump_r2
    kb          <- input$tab5_clump_kb
    use_proxies <- isTRUE(input$tab5_use_proxies)
    proxy_rsq   <- input$tab5_proxy_rsq
    allow_pal   <- isTRUE(input$tab5_allow_palindromes)
    maf_thr     <- input$tab5_maf_thresh
    harm_action <- harmonise_action_value(input$tab5_harm_action)
    method_list <- input$tab5_mr_methods
    if (is.null(method_list) || length(method_list) == 0) method_list <- c("mr_ivw")
    
    tab5_mr_rv$status <- paste0("Running Tab 5 MR (X→M): ", exp2_gcst, " → ", out_gcst, " ...")
    
    ex2_dat <- gwas_call_with_token(
      "tab5",
      extract_instruments(outcomes = e2_id, p1 = p1, clump = do_clump, r2 = r2, kb = kb)
    )
    if (inherits(ex2_dat, "error") || is.null(ex2_dat) || nrow(as.data.frame(ex2_dat)) == 0) {
      tab5_mr_rv$status <- paste0("❌ extract_instruments failed: ",
                                  if (inherits(ex2_dat, "error")) conditionMessage(ex2_dat) else "No instruments.")
      tab5_mr_rv$res <- list(error = tab5_mr_rv$status)
      tab5_mr_rv$hdat <- NULL
      tab5_mr_rv$ivw_summary <- NULL
      tab5_mr_rv$meta <- NULL
      return()
    }
    tab5_mr_rv$snps_inst <- unique(as.character(ex2_dat$SNP))  
    
    out2_dat <- extract_outcome_chunked_cached(
      snps        = ex2_dat$SNP,
      outcome_id = o_id,
      token      = "tab5",
      cache_dir = run_cache_dir,
      proxies    = use_proxies,
      rsq        = proxy_rsq,
      palindromes= allow_pal,
      maf_threshold = maf_thr
    )
    
    # ⛔ RATE LIMIT (MUST BE FIRST)
    if (handle_opengwas_limit(out2_dat, opengwas_reset_rv)) {
      tab5_mr_rv$status <- "⛔ OpenGWAS limit reached. Please wait for reset."
      tab5_mr_rv$res <- NULL
      tab5_mr_rv$hdat <- NULL
      tab5_mr_rv$ivw_summary <- NULL
      tab5_mr_rv$meta <- NULL
      return()
    }
    
    # ❌ REAL ERRORS ONLY
    if (inherits(out2_dat, "error") ||
        is.null(out2_dat) ||
        (is.data.frame(out2_dat) && nrow(out2_dat) == 0)) {
      
      tab5_mr_rv$status <- paste0(
        "❌ extract_outcome_data failed: ",
        if (inherits(out2_dat, "error"))
          conditionMessage(out2_dat)
        else
          "No outcome rows."
      )
      
      tab5_mr_rv$res <- list(error = tab5_mr_rv$status)
      tab5_mr_rv$hdat <- NULL
      tab5_mr_rv$ivw_summary <- NULL
      tab5_mr_rv$meta <- NULL
      return()
    }
    
    hdat <- tryCatch(
      harmonise_data(exposure_dat = ex2_dat, outcome_dat = out2_dat, action = harm_action),
      error = function(e) e
    )
    if (inherits(hdat, "error") || is.null(hdat) || nrow(as.data.frame(hdat)) == 0) {
      tab5_mr_rv$status <- paste0("❌ harmonise_data failed: ",
                                  if (inherits(hdat, "error")) conditionMessage(hdat) else "No harmonised rows.")
      tab5_mr_rv$res <- list(error = tab5_mr_rv$status)
      tab5_mr_rv$hdat <- NULL
      tab5_mr_rv$ivw_summary <- NULL
      tab5_mr_rv$meta <- NULL
      return()
    }
    
    mr_res <- tryCatch(mr(hdat, method_list = method_list), error = function(e) e)
    mr_or  <- tryCatch(generate_odds_ratios(mr_res), error = function(e) e)
    het    <- tryCatch(mr_heterogeneity(hdat, method_list = c("mr_egger_regression", "mr_ivw")), error = function(e) e)
    pleio  <- tryCatch(mr_pleiotropy_test(hdat), error = function(e) e)
    dirres <- tryCatch(directionality_test(hdat), error = function(e) e)
    
    single <- tryCatch(mr_singlesnp(hdat), error = function(e) e)
    loo    <- tryCatch(mr_leaveoneout(hdat), error = function(e) e)
    
    p_scatter <- tryCatch(mr_scatter_plot(mr_res, hdat), error = function(e) e)
    p_forest  <- tryCatch(mr_forest_plot(single), error = function(e) e)
    p_loo     <- tryCatch(mr_leaveoneout_plot(loo), error = function(e) e)
    p_funnel  <- tryCatch(mr_funnel_plot(single), error = function(e) e)
    
    snp_eff <- build_snp_effects_table(hdat)
    
    # --- store full res + harmonised for reuse ---
    tab5_mr_rv$res <- list(
      exposure_gcst = exp2_gcst, outcome_gcst = out_gcst,
      exposure_id = e2_id, outcome_id = o_id,
      harmonised = hdat,
      mr_or = mr_or, heterogeneity = het, pleiotropy = pleio, direction = dirres
    )
    tab5_mr_rv$hdat <- hdat
    tab5_mr_rv$singlesnp <- single
    tab5_mr_rv$loo <- loo
    tab5_mr_rv$plots <- list(scatter=p_scatter, forest=p_forest, loo=p_loo, funnel=p_funnel)
    tab5_mr_rv$snp_eff <- snp_eff
    
    # IVW summary (for mediation reuse)
    ivw_row <- get_ivw_row(as_df_or_msg(mr_or))
    tab5_mr_rv$ivw_summary <- extract_bse(as_df_or_msg(mr_or))
    
    tab5_mr_rv$meta <- list(
      exp2_gcst = exp2_gcst,
      out_gcst  = out_gcst
    )
    
    tab5_mr_rv$status <- "Finished (MR done)"
    session$sendCustomMessage("playSound", list())
  })
  
  # ---- Tab5 plots
  output$tab5_plot_scatter <- renderPlot({
    r <- tab5_mr_rv$plots$scatter
    if (is.null(r) || inherits(r, "error")) return(invisible())
    if (is.list(r) && length(r) >= 1) print(r[[1]]) else print(r)
  })
  output$tab5_plot_forest <- renderPlot({
    r <- tab5_mr_rv$plots$forest
    if (is.null(r) || inherits(r, "error")) return(invisible())
    if (is.list(r) && length(r) >= 1) print(r[[1]]) else print(r)
  })
  output$tab5_plot_loo <- renderPlot({
    r <- tab5_mr_rv$plots$loo
    if (is.null(r) || inherits(r, "error")) return(invisible())
    if (is.list(r) && length(r) >= 1) print(r[[1]]) else print(r)
  })
  output$tab5_plot_funnel <- renderPlot({
    r <- tab5_mr_rv$plots$funnel
    if (is.null(r) || inherits(r, "error")) return(invisible())
    if (is.list(r) && length(r) >= 1) print(r[[1]]) else print(r)
  })
  
  # ---- Tab5 tables
  # ---- TAB 5: MR results table (Exposure 1 -> Mediator)
  output$tbl_tab5_mr_results <- renderDT({
    
    r <- tab5_mr_rv$res
    
    if (is.null(r)) {
      return(datatable(
        data.frame(Message = "No results yet."),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    if (!is.null(r$error)) {
      return(datatable(
        data.frame(Message = r$error),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    df <- as_df_or_msg(r$mr_or)
    
    if (nrow(df) == 0 || "Message" %in% names(df)) {
      return(datatable(
        data.frame(Message = "No MR results."),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    # ---------------------------------------
    # Add trait labels
    # ---------------------------------------
    df$id.exposure <- get_trait_label(r$exposure_gcst, DATA$merged)
    df$id.outcome  <- get_trait_label(r$outcome_gcst,  DATA$merged)
    
    # Remove default cols if present
    df$exposure <- NULL
    df$outcome  <- NULL
    
    # ---------------------------------------
    # Convert numeric
    # ---------------------------------------
    df$b  <- as.numeric(df$b)
    df$se <- as.numeric(df$se)
    
    # ---------------------------------------
    # Detect if Mediator is binary (Tab5 outcome = mediator)
    # ---------------------------------------
    is_binary_mediator <- FALSE
    med_id <- stringr::str_trim(as.character(r$outcome_gcst))
    
    if (!is.na(med_id) && med_id != "") {
      row <- DATA$merged %>%
        dplyr::filter(`STUDY ACCESSION` == med_id)
      
      if (nrow(row) > 0) {
        is_binary_mediator <- identical(as.character(row$`Data Type`[1]), "Binary")
      }
    }
    
    z <- 1.96
    
    # ---------------------------------------
    # Build Beta CI (always shown)
    # ---------------------------------------
    lci_beta <- df$b - z * df$se
    uci_beta <- df$b + z * df$se
    
    df$`Beta (95% CI)` <- sprintf(
      "%.4f (%.4f, %.4f)",
      df$b, lci_beta, uci_beta
    )
    
    # ---------------------------------------
    # Build OR CI (only if mediator is binary)
    # ---------------------------------------
    if (is_binary_mediator) {
      or     <- exp(df$b)
      lci_or <- exp(lci_beta)
      uci_or <- exp(uci_beta)
      
      df$`Odds Ratio (95% CI)` <- sprintf(
        "%.4f (%.4f, %.4f)",
        or, lci_or, uci_or
      )
    }
    
    # ---------------------------------------
    # Add significance column
    # ---------------------------------------
    df <- add_significance_column(df)
    
    # ---------------------------------------
    # Keep final display columns
    # ---------------------------------------
    keep_cols <- c(
      "id.exposure",
      "id.outcome",
      "method",
      "nsnp",
      "Beta (95% CI)",
      "se"
    )
    
    if (is_binary_mediator) {
      keep_cols <- c(keep_cols, "Odds Ratio (95% CI)")
    }
    
    keep_cols <- c(keep_cols, "pval", "Significance")
    keep_cols <- intersect(keep_cols, names(df))
    
    df <- df[, keep_cols, drop = FALSE]
    
    # ---------------------------------------
    # Rename columns nicely
    # ---------------------------------------
    rename_map <- c(
      id.exposure = "Exposure (Trait | GCST ID)",
      id.outcome  = "Mediator (Trait | GCST ID)",
      method      = "MR Method",
      nsnp        = "Number of SNPs",
      se          = "Standard Error",
      pval        = "P-value",
      Significance = "Significance"
    )
    
    old_names <- names(df)
    new_names <- ifelse(old_names %in% names(rename_map), rename_map[old_names], old_names)
    names(df) <- new_names
    
    # ---------------------------------------
    # Round p-value + numeric columns to 4 dp (display)
    # ---------------------------------------
    if ("P-value" %in% names(df)) {
      df$`P-value` <- round(as.numeric(df$`P-value`), 4)
    }
    if ("Standard Error" %in% names(df)) {
      df$`Standard Error` <- round(as.numeric(df$`Standard Error`), 4)
    }
    if ("Number of SNPs" %in% names(df)) {
      df$`Number of SNPs` <- as.integer(df$`Number of SNPs`)
    }
    
    datatable(
      df,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE, dom = "tip")
    )
  })
  
  output$tbl_tab5_het <- renderDT({
    
    r <- tab5_mr_rv$res
    
    if (is.null(r)) {
      return(datatable(data.frame(Message="No results."),
                       options=list(dom="t"), rownames=FALSE))
    }
    
    df <- as_df_or_msg(r$heterogeneity)
    
    # ADD TRAIT LABELS
    df$id.exposure <- get_trait_label(r$exposure_gcst, DATA$merged)
    df$id.outcome  <- get_trait_label(r$outcome_gcst,  DATA$merged)
    
    # REMOVE
    df$exposure <- NULL
    df$outcome  <- NULL
    
    df <- add_significance_column(df)
    df <- round_df_display(df,4)
    rename_map <- c(
      id.exposure = "Exposure (Trait | GCST ID)",
      id.outcome  = "Outcome (Trait | GCST ID)",
      method      = "MR Method",
      Q           = "Cochran’s Q Statistic",
      Q_df        = "Degrees of Freedom",
      Q_pval      = "Q-test P-value",
      Significance = "Heterogeneity Significance"
    )
    
    colnames(df) <- rename_map[colnames(df)]
    datatable(df, rownames=FALSE,
              options=list(pageLength=10, scrollX=TRUE, dom="tip"))
  })
  
  output$tbl_tab5_pleio <- renderDT({
    
    r <- tab5_mr_rv$res
    
    if (is.null(r)) {
      return(datatable(data.frame(Message="No results."),
                       options=list(dom="t"), rownames=FALSE))
    }
    
    df <- as_df_or_msg(r$pleiotropy)
    
    # ADD TRAIT LABELS
    df$id.exposure <- get_trait_label(r$exposure_gcst, DATA$merged)
    df$id.outcome  <- get_trait_label(r$outcome_gcst,  DATA$merged)
    
    # EMOVE
    df$exposure <- NULL
    df$outcome  <- NULL
    
    df <- add_significance_column(df)
    df <- round_df_display(df,4)
    rename_map <- c(
      id.exposure = "Exposure (Trait | GCST ID)",
      id.outcome  = "Outcome (Trait | GCST ID)",
      egger_intercept = "MR-Egger Intercept",
      se = "Standard Error",
      pval = "P-value",
      Significance = "Pleiotropy Significance"
    )
    
    colnames(df) <- rename_map[colnames(df)]
    datatable(df, rownames=FALSE,
              options=list(pageLength=10, scrollX=TRUE, dom="tip"))
  })
  
  output$tbl_tab5_dir <- renderDT({
    
    r <- tab5_mr_rv$res
    
    if (is.null(r)) {
      return(datatable(
        data.frame(Message="No results."),
        options=list(dom="t"),
        rownames=FALSE
      ))
    }
    
    df <- as_df_or_msg(r$direction)
    
    if (is.null(df) || nrow(df)==0 || "Message" %in% names(df)) {
      return(datatable(
        data.frame(Message="No directionality results."),
        options=list(dom="t"),
        rownames=FALSE
      ))
    }
    
    # Add trait labels
    df$id.exposure <- get_trait_label(r$exposure_gcst, DATA$merged)
    df$id.outcome  <- get_trait_label(r$outcome_gcst,  DATA$merged)
    
    # Build clean display table
    df_display <- df %>%
      dplyr::mutate(
        `Exposure (Trait | GCST ID)` = id.exposure,
        `Outcome (Trait | GCST ID)`  = id.outcome,
        `R² (Exposure)` = round(snp_r2.exposure, 4),
        `R² (Outcome)`  = round(snp_r2.outcome, 4),
        `Correct Direction?` = ifelse(correct_causal_direction, "Yes", "No"),
        `Steiger P-value` = signif(steiger_pval, 3)
      ) %>%
      dplyr::select(
        `Exposure (Trait | GCST ID)`,
        `Outcome (Trait | GCST ID)`,
        `R² (Exposure)`,
        `R² (Outcome)`,
        `Correct Direction?`,
        `Steiger P-value`
      )
    
    datatable(
      df_display,
      rownames=FALSE,
      options=list(pageLength=10, scrollX=TRUE, dom="tip")
    )
  })
  output$tab5_pleio_summary <- renderText({
    r <- tab5_mr_rv$res
    if (is.null(r)) return("")
    df <- as_df_or_msg(r$pleiotropy)
    if (is.null(df) || nrow(df) == 0 || "Message" %in% names(df)) return("")
    intercept <- df$egger_intercept[1] %||% df$intercept[1] %||% NA
    se <- df$se[1] %||% NA
    p <- df$pval[1] %||% NA
    paste0("Egger intercept: ", format(as.numeric(intercept), digits=4), "\n",
           "SE: ", format(as.numeric(se), digits=4), "\n",
           "P: ", format(as.numeric(p), digits=4))
  })
  
  tab5_rv <- reactiveValues(
    df = data.frame(Message = "Click 'Build SNP Table' to create it.", stringsAsFactors = FALSE)
  )
  # ==========================================================
  # TAB 5B — Download SNP Table for MVMR
  # ==========================================================
  output$dl_tab5_snp_mvmr_xlsx <- downloadHandler(
    
    filename = function() {
      paste0(
        "TAB5_MVMR_SNP_Table_",
        Sys.Date(),
        ".xlsx"
      )
    },
    
    content = function(file) {
      
      df <- tab5_rv$df
      
      if (is.null(df) || nrow(df) == 0) {
        df <- data.frame(Message = "No SNP table built yet.")
      }
      
      openxlsx::write.xlsx(df, file)
    }
  )
  # ==========================================================
  # TAB 5B: Build SNP table for MVMR (REUSE Tab5A when possible)
  # ==========================================================
  observeEvent(input$build_tab5, {
    
    show_busy_modal("Building SNP table...", "Building UNION-based MVMR table (Exp1 ∪ Exp2) …")
    on.exit({ removeModal() }, add = TRUE)
    
    exp1_gcst <- unique(str_trim(as.character(exposure_mod$selected_gcst())))[1]
    out_gcst  <- unique(str_trim(as.character(outcome_mod$selected_gcst())))[1]
    exp2_gcst <- unique(str_trim(as.character(exp2_for_mvmr())))[1]
    
    if (is.na(exp1_gcst) || exp1_gcst == "" ||
        is.na(out_gcst)  || out_gcst  == "" ||
        is.na(exp2_gcst) || exp2_gcst == "") {
      tab5_rv$df <- data.frame(Message = "Must select Exposure 1 (Tab1), Outcome (Tab2) and Exposure 2 (Tab5).")
      return()
    }
    
    # ---- Must have Tab3 result for Exposure1 -> Outcome (to reuse res1 metadata if available)
    key_pair <- paste0(exp1_gcst, "  →  ", out_gcst)
    res1 <- mr_rv$results_by_pair[[key_pair]]
    if (is.null(res1) || !is.null(res1$error)) {
      tab5_rv$df <- data.frame(Message = "Run Tab 3 MR for this Exposure 1 → Outcome pair first (and select the same pair).")
      return()
    }
    
    h1 <- as.data.frame(res1$harmonised)
    if (is.null(h1) || nrow(h1) == 0) {
      tab5_rv$df <- data.frame(Message = "Tab 3 harmonised data is empty for this pair.")
      return()
    }
    
    # ============================
    # UNION-based SNP table build
    # ============================
    
    # ---- 1) instrument SNP sets (prefer cached instruments)
    snps_exp1 <- res1$snps_inst %||% unique(as.character(h1$SNP))
    snps_exp1 <- unique(str_trim(snps_exp1))
    snps_exp1 <- snps_exp1[!is.na(snps_exp1) & snps_exp1 != ""]
    
    # Exposure2 instruments: prefer cached instruments (if you store them), otherwise extract instruments now
    snps_exp2 <- NULL
    can_reuse_exp2_inst <- FALSE
    
    if (!is.null(tab5_mr_rv$snps_inst) && !is.null(tab5_mr_rv$meta)) {
      if (identical(str_trim(tab5_mr_rv$meta$exp2_gcst), str_trim(exp2_gcst))) {
        snps_exp2 <- tab5_mr_rv$snps_inst
        can_reuse_exp2_inst <- TRUE
      }
    }
    
    if (!can_reuse_exp2_inst) {
      e2_id <- resolve_opengwas_id(exp2_gcst, "tab5")
      if (is.na(e2_id) || e2_id == "") {
        tab5_rv$df <- data.frame(Message = paste0("Could not resolve OpenGWAS ID for Exposure 2: ", exp2_gcst))
        return()
      }
      
      ex2_inst <- gwas_call_with_token(
        "tab5",
        extract_instruments(
          outcomes = e2_id,
          p1 = input$p_thresh %||% 5e-8,
          clump = isTRUE(input$do_clump),
          r2 = input$clump_r2,
          kb = input$clump_kb
        )
      )
      
      if (inherits(ex2_inst, "error") || is.null(ex2_inst) || nrow(as.data.frame(ex2_inst)) == 0) {
        tab5_rv$df <- data.frame(Message = paste0(
          "Exposure 2 extract_instruments failed: ",
          if (inherits(ex2_inst, "error")) conditionMessage(ex2_inst) else "No instruments returned."
        ))
        return()
      }
      
      snps_exp2 <- unique(str_trim(as.character(ex2_inst$SNP)))
      snps_exp2 <- snps_exp2[!is.na(snps_exp2) & snps_exp2 != ""]
      
      tab5_mr_rv$snps_inst <- snps_exp2
    }
    
    snps_all <- unique(c(snps_exp1, snps_exp2))
    if (length(snps_all) == 0) {
      tab5_rv$df <- data.frame(Message = "No SNPs found in instrument sets for Exposure 1 or Exposure 2.")
      return()
    }
    
    # ---- 2) helper: extract associations for ANY trait (exp1/exp2/outcome)
    fetch_assoc <- function(snps, gcst, token_alias, cache_dir,
                            proxies=TRUE, rsq=0.8, palindromes=TRUE, maf_threshold=0.30) {
      
      gid <- resolve_opengwas_id(gcst, token_alias)
      if (is.na(gid) || gid == "") {
        return(structure(simpleError(paste0("Could not resolve OpenGWAS ID: ", gcst)), class = "error"))
      }
      
      # Build args; ONLY pass cache_dir if your function supports it (avoids unused argument error)
      args <- list(
        snps = snps,
        outcome_id = gid,
        token = token_alias,
        proxies = proxies,
        rsq = rsq,
        palindromes = palindromes,
        maf_threshold = maf_threshold
      )
      fml <- names(formals(extract_outcome_chunked_cached))
      if ("cache_dir" %in% fml) args$cache_dir <- cache_dir
      
      do.call(extract_outcome_chunked_cached, args)
    }
    
    # ---- 3) fetch Outcome + both exposures for ALL union SNPs
    run_cache_dir2 <- file.path(tempdir(), paste0("mvmr_union_", as.integer(Sys.time()), "_", sample(1e6,1)))
    dir.create(run_cache_dir2, recursive = TRUE, showWarnings = FALSE)
    
    tok <- "tab5"
    use_proxies2 <- isTRUE(input$use_proxies)
    proxy_rsq2   <- input$proxy_rsq
    allow_pal2   <- isTRUE(input$allow_palindromes)
    maf_thr2     <- input$maf_thresh
    
    out_assoc <- fetch_assoc(snps_all, out_gcst, tok, run_cache_dir2,
                             proxies=use_proxies2, rsq=proxy_rsq2, palindromes=allow_pal2, maf_threshold=maf_thr2)
    if (handle_opengwas_limit(out_assoc, opengwas_reset_rv)) {
      tab5_rv$df <- data.frame(Message = "⛔ OpenGWAS rate limit reached while fetching Outcome associations.")
      return()
    }
    if (inherits(out_assoc, "error") || is.null(out_assoc) || nrow(as.data.frame(out_assoc)) == 0) {
      tab5_rv$df <- data.frame(Message = "Failed to fetch Outcome associations for union SNP set.")
      return()
    }
    
    x1_assoc <- fetch_assoc(snps_all, exp1_gcst, tok, run_cache_dir2,
                            proxies=use_proxies2, rsq=proxy_rsq2, palindromes=allow_pal2, maf_threshold=maf_thr2)
    if (handle_opengwas_limit(x1_assoc, opengwas_reset_rv)) {
      tab5_rv$df <- data.frame(Message = "⛔ OpenGWAS rate limit reached while fetching Exposure 1 associations.")
      return()
    }
    if (inherits(x1_assoc, "error") || is.null(x1_assoc) || nrow(as.data.frame(x1_assoc)) == 0) {
      tab5_rv$df <- data.frame(Message = "Failed to fetch Exposure 1 associations for union SNP set.")
      return()
    }
    
    x2_assoc <- fetch_assoc(snps_all, exp2_gcst, tok, run_cache_dir2,
                            proxies=use_proxies2, rsq=proxy_rsq2, palindromes=allow_pal2, maf_threshold=maf_thr2)
    if (handle_opengwas_limit(x2_assoc, opengwas_reset_rv)) {
      tab5_rv$df <- data.frame(Message = "⛔ OpenGWAS rate limit reached while fetching Exposure 2 associations.")
      return()
    }
    if (inherits(x2_assoc, "error") || is.null(x2_assoc) || nrow(as.data.frame(x2_assoc)) == 0) {
      tab5_rv$df <- data.frame(Message = "Failed to fetch Exposure 2 associations for union SNP set.")
      return()
    }
    
    # ---- 4) standardise columns
    Y <- as.data.frame(out_assoc) %>%
      transmute(
        SNP = str_trim(as.character(SNP)),
        EA_ref = effect_allele.outcome,
        OA_ref = other_allele.outcome,
        Outcome_beta = as.numeric(beta.outcome),
        Outcome_se   = as.numeric(se.outcome)
      ) %>%
      filter(!is.na(SNP), SNP != "") %>% distinct(SNP, .keep_all = TRUE)
    
    X1 <- as.data.frame(x1_assoc) %>%
      transmute(
        SNP = str_trim(as.character(SNP)),
        EA1 = effect_allele.outcome,
        OA1 = other_allele.outcome,
        B1  = as.numeric(beta.outcome),
        SE1 = as.numeric(se.outcome)
      ) %>%
      filter(!is.na(SNP), SNP != "") %>% distinct(SNP, .keep_all = TRUE)
    
    X2 <- as.data.frame(x2_assoc) %>%
      transmute(
        SNP = str_trim(as.character(SNP)),
        EA2 = effect_allele.outcome,
        OA2 = other_allele.outcome,
        B2  = as.numeric(beta.outcome),
        SE2 = as.numeric(se.outcome)
      ) %>%
      filter(!is.na(SNP), SNP != "") %>% distinct(SNP, .keep_all = TRUE)
    
    # ---- 5) join by SNP (UNION table)
    dfu <- tibble::tibble(SNP = snps_all) %>%
      left_join(Y,  by = "SNP") %>%
      left_join(X1, by = "SNP") %>%
      left_join(X2, by = "SNP")
    
    # Drop SNPs missing any needed association
    dfu <- dfu %>%
      filter(
        is.finite(Outcome_beta), is.finite(Outcome_se),
        is.finite(B1), is.finite(SE1),
        is.finite(B2), is.finite(SE2),
        !is.na(EA_ref), !is.na(OA_ref),
        !is.na(EA1), !is.na(OA1),
        !is.na(EA2), !is.na(OA2)
      )
    
    if (nrow(dfu) == 0) {
      tab5_rv$df <- data.frame(Message = "After fetching, no SNPs had complete Exp1+Exp2+Outcome associations.")
      return()
    }
    
    # ---- 6) remove palindromic ambiguous SNPs using reference alleles (Outcome alleles)
    pal_ref <- is_palindromic(dfu$EA_ref, dfu$OA_ref)
    dfu <- dfu[!pal_ref, , drop = FALSE]
    
    if (nrow(dfu) == 0) {
      tab5_rv$df <- data.frame(Message = "All SNPs removed due to palindromic ambiguity (reference alleles).")
      return()
    }
    
    # ---- 7) align Exposure betas to the same reference allele pair (Outcome alleles)
    a1 <- align_to_ref(
      beta   = dfu$B1,
      ea     = dfu$EA1,
      oa     = dfu$OA1,
      ref_ea = dfu$EA_ref,
      ref_oa = dfu$OA_ref
    )
    
    a2 <- align_to_ref(
      beta   = dfu$B2,
      ea     = dfu$EA2,
      oa     = dfu$OA2,
      ref_ea = dfu$EA_ref,
      ref_oa = dfu$OA_ref
    )
    
    dfu$Exposure1_beta <- a1$beta
    dfu$Exposure2_beta <- a2$beta
    
    # any SNPs that could not be aligned become NA -> drop
    dfu <- dfu %>% filter(is.finite(Exposure1_beta), is.finite(Exposure2_beta))
    
    if (nrow(dfu) == 0) {
      tab5_rv$df <- data.frame(Message = "No SNPs could be aligned consistently across Exp1/Exp2 to Outcome alleles.")
      return()
    }
    
    # ---- 8) final MVMR table
    result_df <- dfu %>%
      transmute(
        SNP,
        Outcome_beta = as.numeric(Outcome_beta),
        Outcome_se   = as.numeric(Outcome_se),
        Exposure1_beta = as.numeric(Exposure1_beta),
        Exposure1_se   = as.numeric(SE1),
        Exposure2_beta = as.numeric(Exposure2_beta),
        Exposure2_se   = as.numeric(SE2)
      ) %>%
      distinct(SNP, .keep_all = TRUE)
    
    warn_msg <- ""
    if (nrow(result_df) < 3) {
      warn_msg <- paste0("⚠️ Only ", nrow(result_df), " SNPs after UNION+alignment. MVMR may be unstable.")
    }
    attr(result_df, "warn") <- warn_msg
    
    tab5_rv$df <- result_df
    session$sendCustomMessage("playSound", list())
  })
  
  output$tab5_warn <- renderText({
    df <- tab5_rv$df
    if (is.null(df)) return("")
    w <- attr(df, "warn")
    if (is.null(w) || is.na(w) || w == "") return("")
    w
  })
  # ===============================
  # TAB 5B: has SNP table flag
  # ===============================
  output$tab5_snp_ready <- reactive({
    
    df <- tab5_rv$df
    
    !is.null(df) &&
      !"Message" %in% names(df) &&
      nrow(df) > 0
  })
  outputOptions(output, "tab5_snp_ready", suspendWhenHidden = FALSE)
  
  output$tbl_tab5_snp_mvmr <- renderDT({
    df <- tab5_rv$df
    
    if (is.null(df) || nrow(df) == 0) {
      return(datatable(data.frame(Message = "No data yet. Click 'Build SNP Table'."),
                       rownames = FALSE, options = list(dom = "t")))
    }
    if ("Message" %in% names(df)) {
      return(datatable(df, rownames = FALSE, options = list(dom = "t")))
    }
    
    datatable(df, rownames = FALSE, options = list(pageLength = 12, scrollX = TRUE, dom = "tip"))
  })
  
  
  # --------------------------------------------------------------------
  # TAB 7: Mediation (2SMR)
  # --------------------------------------------------------------------
  
  # --- Tab7 mediator search module (only active if "new_search") ---
  # Reactive tracking second checkbox for mediator: filter to non‑overlapping cohorts
  tab7_only_non_overlap <- reactive({ isTRUE(input$tab7_non_overlap) })
  
  # Cohort overlap function for mediators: compares candidate mediators against both
  # Exposure 1 and Outcome cohorts. Returns "YES" if any cohort matches, "NO"
  # if none match, and "UNKNOWN" if cohort data is missing.
  mediator_overlap_fun <- reactive({
    # Build reference cohort set from exposure 1 and outcome GCSTs
    exp_ref <- tab7_ref_ids()
    exp_ref <- exp_ref[exp_ref != "" & !is.na(exp_ref)]
    out_ref <- unique(stringr::str_trim(as.character(outcome_mod$selected_gcst())))
    out_ref <- out_ref[out_ref != "" & !is.na(out_ref)]
    ref <- c(exp_ref, out_ref)
    ref <- ref[ref != "" & !is.na(ref)]
    if (length(ref) == 0) return(NULL)
    ref_cohorts <- unique(unlist(lapply(ref, function(gid) {
      rows <- which(DATA$merged$`STUDY ACCESSION` == gid)
      if (length(rows) == 0) return(character(0))
      co_strs <- DATA$merged$COHORT[rows]
      co_strs <- co_strs[!is.na(co_strs) & co_strs != ""]
      if (length(co_strs) == 0) return(character(0))
      unique(unlist(lapply(co_strs, parse_cohort)))
    })))
    # If no reference cohorts, return UNKNOWN for all
    if (length(ref_cohorts) == 0) {
      return(function(gcs) rep("UNKNOWN", length(gcs)))
    }
    function(gcs) {
      vapply(gcs, function(gid) {
        rows <- which(DATA$merged$`STUDY ACCESSION` == gid)
        if (length(rows) == 0) return("UNKNOWN")
        co_strs <- DATA$merged$COHORT[rows]
        co_strs <- co_strs[!is.na(co_strs) & co_strs != ""]
        if (length(co_strs) == 0) return("UNKNOWN")
        cand_coh <- unique(unlist(lapply(co_strs, parse_cohort)))
        # Treat presence of 'MULTIPLE' in candidate or reference as overlapping
        if ("MULTIPLE" %in% cand_coh || "MULTIPLE" %in% ref_cohorts) {
          return("YES")
        }
        if (length(intersect(ref_cohorts, cand_coh)) > 0) "YES" else "NO"
      }, character(1))
    }
  })
  
  tab7_search <- searchServer(
    id = "tab7_mediator_search",
    merged = DATA$merged,
    wide   = DATA$wide,
    snp_counts = DATA$snp_counts,
    status_fun_reactive = reactive({
      function(gcsts) {
        ref <- tab7_ref_ids()
        if (length(ref) == 0) return(rep("Unknown", length(gcsts)))
        
        vapply(gcsts, function(g) {
          ok <- is_independent_from_all(g, ref, DATA$merged)
          if (is.na(ok)) "Unknown" else if (ok) "Independent" else "Not independent"
        }, character(1))
      }
    }),
    # Use first checkbox to filter to same‑ancestry mediators (Not independent)
    filter_only_indep_reactive = reactive({ isTRUE(input$tab7_only_indep) }),
    overlap_fun_reactive = mediator_overlap_fun,
    filter_only_non_overlap_reactive = tab7_only_non_overlap,
    enable_selection = TRUE,
    # Restrict mediators: when the first checkbox is selected, return GCSTs that
    # share ancestry with the reference set (i.e., are not independent)
    restrict_accessions_reactive = reactive({
      if (!isTRUE(input$tab7_only_indep)) {
        return(character(0))
      }
      ref <- tab7_ref_ids()
      if (length(ref) == 0) return(character(0))
      all_gcst <- unique(as.character(DATA$merged$`STUDY ACCESSION`))
      all_gcst <- all_gcst[!is.na(all_gcst) & all_gcst != ""]
      ok <- vapply(all_gcst, function(g) {
        ans <- is_independent_from_all(g, ref, DATA$merged)
        !(isTRUE(ans))
      }, logical(1))
      all_gcst[ok]
    })
  )
  
  # --------------------------------------------------
  # Binary Mediator Notification (Tab 7)
  # --------------------------------------------------
  
  observeEvent(tab7_search$selected_gcst(), {
    
    ids <- tab7_search$selected_gcst()
    if (length(ids) == 0) return()
    
    mediator_id <- ids[1]
    
    row <- DATA$merged %>%
      dplyr::filter(`STUDY ACCESSION` == mediator_id)
    
    if (nrow(row) == 0) return()
    
    mediator_type <- row$`Data Type`[1]
    
    if (!is.na(mediator_type) && mediator_type == "Binary") {
      
      showNotification(
        paste0(
          "⚠ Exposure2 / Mediator ", mediator_id, " (",
          row$MAPPED_TRAIT[1],
          ") is binary. ",
          "X → M estimates are in log-odds scale. ",
          "Use exp(beta) for Odds Ratio interpretation."
        ),
        type = "warning",
        duration = 8
      )
      
    }
    
  })
  # ------------------------------------------------------------
  # AUTO-CLOSE mediator search modal when a mediator is selected
  # ------------------------------------------------------------
  observeEvent(tab7_search$selected_gcst(), {
    sel <- tab7_search$selected_gcst()
    
    if (length(sel) > 0 && !is.na(sel[1]) && sel[1] != "") {
      removeModal()
    }
  }, ignoreInit = TRUE)
  
  output$ui_tab7_mediator_picker <- renderUI({
    req(input$mediator_source)
    
    if (input$mediator_source == "from_tabA") {
      tagList(
        div(class="label", ""),
        verbatimTextOutput("txt_tab7_mediator_from_tab5")  # you can rename later
      )
    } else {
      div(class="label", "")
    }
  })
  
  output$txt_tab7_mediator_from_tab5 <- renderText({
    x <- exposure2_selected()[1]
    if (is.null(x) || is.na(x) || x == "") return("Exposure 2 not selected yet in Tab 5.")
    x
  })
  
  
  
  # ============================================================
  # Tab 7 (Mediation): Mediator → Outcome pair display
  # ============================================================
  
  tab7_my_pair_lhs <- reactive({
    req(mediator_gcst())
    mediator_gcst()
  })
  
  tab7_my_pair_rhs <- reactive({
    req(outcome_mod$selected_gcst())
    unique(str_trim(as.character(outcome_mod$selected_gcst())))[1]
  })
  
  output$tab7_my_pair_ui <- render_pair_box(
    tab7_my_pair_lhs,
    tab7_my_pair_rhs
  )
  
  observeEvent(mediator_gcst(), {
    
  }, ignoreInit = TRUE)
  
  # ============================================================
  # Tab 7B (Mediation): Exposure 1 → Mediator pair display
  # ============================================================
  
  tab7_x_pair_lhs <- reactive({
    req(exposure_mod$selected_gcst())
    unique(str_trim(as.character(exposure_mod$selected_gcst())))[1]
  })
  
  tab7_x_pair_rhs <- reactive({
    req(mediator_gcst())
    mediator_gcst()
  })
  
  output$tab7_x_pair_ui <- render_pair_box(
    tab7_x_pair_lhs,
    tab7_x_pair_rhs
  )
  
  output$ui_tab7_main_tabs <- renderUI({
    req(input$mediator_source)
    
    tab_C <- tabPanel(
      "Mediation effect",
      
      div(
        class = "card",
        
        # ================= HEADER =================
        div(
          class = "card-header",
          div(class = "card-title", "Calculate Mediation effect"),
          div(
            class = "card-subtitle",
            "Estimates indirect, direct, and total effects using IVW-based mediation with optional Sobel,bootstrap confidence intervals."
          )
        ),
        
        # ================= BODY =================
        div(
          class = "card-body",
          
          fluidRow(
            
            # ==================================================
            # LEFT COLUMN
            # ==================================================
            column(
              width = 3,
              
              div(
                class = "panel",
                div(class = "panel-title", "Options"),
                
                radioButtons(
                  "tab7_ci_method",
                  label = "Confidence interval method",
                  choices = c(
                    "Delta method (Sobel approximation)" = "delta",
                    "Bootstrap (robust)" = "bootstrap",
                    "Report both" = "both"
                  ),
                  selected = "delta"
                ),
                
                conditionalPanel(
                  condition = "input.tab7_ci_method == 'bootstrap' || input.tab7_ci_method == 'both'",
                  numericInput(
                    "tab7_boot_B",
                    "Bootstrap replicates (B)",
                    value = 1000,
                    min = 200,
                    step = 100
                  )
                )
                
              ),
              
              div(class="sep"),
              
              div(
                class = "panel",
                
                actionButton(
                  "run_mediation_2smr",
                  "⚡ Run Mediation Model",
                  class = "btn-primary w100"
                ),
                
                div(class="sep"),
                
                div(
                  class = "panel-title",
                  "Downloads"
                ),
                
                conditionalPanel(
                  condition = "output.tab7c_has_results === true",
                  
                  downloadButton(
                    "dl_tab7c_mediation_xlsx",
                    "⬇ Mediation results (Excel)",
                    class = "btn-ghost w100"
                  )
                  
                ),
                
                conditionalPanel(
                  condition = "output.tab7c_has_results === false",
                  tags$p(
                    style="opacity:0.6; font-style:italic;",
                    "Downloads will appear after running mediation."
                  )
                )
              )
            ),
            
            # ==================================================
            # RIGHT COLUMN
            # ==================================================
            column(
              width = 9,
              
              div(
                class = "panel",
                
                div(class = "panel-title", "Mediation Results & Statistics"),
                
                # ------------------------------------------------
                # Placeholder BEFORE running mediation
                # ------------------------------------------------
                conditionalPanel(
                  condition = "output.tab7c_has_results === false",
                  tags$p(
                    style = "margin-top:10px; opacity:0.7; font-style:italic;",
                    "Mediation Results & Statistics will appear here after running the Mediation analysis."
                  )
                ),
                
                # ------------------------------------------------
                #  Actual content AFTER mediation run
                # ------------------------------------------------
                conditionalPanel(
                  condition = "output.tab7c_has_results === true",
                  
                  div(
                    class="note",
                    "Two complementary mediation estimators:",
                    tags$ul(
                      tags$li("Difference method (MVMR-based decomposition)"),
                      tags$li("Product method (Two-sample MR-based mediation)")
                    )
                  ),
                  
                  div(class="sep"),
                  
                  # ---------- MVMR Difference ----------
                  div(
                    style="display:flex; align-items:center; justify-content:space-between;",
                    tags$h4("MVMR Mediation (Difference Method)"),
                    actionButton("info_mvmr_diff",
                                 label=NULL,
                                 icon=icon("info-circle"),
                                 class="btn-link")
                  ),
                  
                  div(class="note",
                      "Indirect effect calculated as Total − Direct from MVMR."
                  ),
                  
                  DTOutput("tbl_mediation_mvmr"),
                  
                  div(class="sep"),
                  
                  # ---------- Two-Sample Product ----------
                  div(
                    style="display:flex; align-items:center; justify-content:space-between;",
                    tags$h4("Two-Sample MR Mediation (Product Method)"),
                    actionButton("info_mvmr_prod",
                                 label=NULL,
                                 icon=icon("info-circle"),
                                 class="btn-link")
                  ),
                  
                  div(class="note",
                      "Indirect effect calculated as (X→M) × (M→Y|X)."
                  ),
                  
                  DTOutput("tbl_mediation_2smr"),
                  
                  div(class="sep"),
                  
                  # ---------- Diagram ----------
                  div(
                    class = "panel",
                    
                    div(class = "panel-title", "Causal Pathway Diagram"),
                    
                    div(
                      class = "note",
                      "This diagram illustrates the hypothesized mediation structure."
                    ),
                    
                    grVizOutput("plot_tab7_diagram", height = "300px")
                  )
                )
              )
            )
          )
        )
      )
    )
    
    tabsetPanel(
      id = "mediation_step",
      tab_C
    )
  })

  # -----------------------------
  # SAFE DT + Plot render helpers
  # -----------------------------
  .safe_dt <- function(x, msg = "No results yet.") {
    if (is.null(x)) {
      return(DT::datatable(data.frame(Message = msg), rownames = FALSE, options = list(dom="t")))
    }
    if (is.data.frame(x) && "Message" %in% names(x)) {
      return(DT::datatable(x, rownames = FALSE, options = list(dom="t")))
    }
    DT::datatable(x, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE, dom="tip"))
  }
  
  .safe_df <- function(obj, fallback = "No results yet.") {
    df <- as_df_or_msg(obj)
    if (is.null(df) || nrow(df) == 0) return(data.frame(Message = fallback))
    df
  }
  
  .safe_plot_print <- function(p) {
    if (is.null(p) || inherits(p, "error")) return(invisible())
    if (is.list(p) && length(p) >= 1) print(p[[1]]) else print(p)
  }
  
  
  .output_plot <- function(p) {
    if (is.null(p) || inherits(p,"error")) return(invisible())
    if (is.list(p) && length(p)>=1) print(p[[1]]) else print(p)
  }
  
  # ---- Helpers (local; no package loads needed)
  get_ivw_row <- function(mr_df) {
    d <- as.data.frame(mr_df)
    if (is.null(d) || nrow(d) == 0) return(NULL)
    if (!"method" %in% names(d)) return(d[1, , drop = FALSE])
    
    ix <- which(d$method %in% c("Inverse variance weighted", "mr_ivw", "IVW"))
    if (length(ix) == 0) d[1, , drop = FALSE] else d[ix[1], , drop = FALSE]
  }
  
  extract_bse <- function(mr_df) {
    row <- get_ivw_row(mr_df)
    if (is.null(row)) return(list(b = NA_real_, se = NA_real_, p = NA_real_))
    
    b  <- suppressWarnings(as.numeric(row$b[1]))
    se <- suppressWarnings(as.numeric(row$se[1]))
    p  <- suppressWarnings(as.numeric((row$pval %||% NA)[1]))
    
    if (!is.finite(b))  b  <- NA_real_
    if (!is.finite(se)) se <- NA_real_
    if (!is.finite(p))  p  <- NA_real_
    
    list(b = b, se = se, p = p)
  }
  
  # Extract MVMR IVW row by which exposure (1 = X, 2 = M)
  get_mvmr_bse <- function(ivw_df, which_exposure = 1) {
    d <- as.data.frame(ivw_df)
    if (is.null(d) || nrow(d) == 0) return(list(b=NA_real_, se=NA_real_, p=NA_real_))
    
    if (!"Exposure ID" %in% names(d) && !is.null(rownames(d))) {
      d <- cbind(`Exposure ID` = rownames(d), d)
      rownames(d) <- NULL
    }
    
    b_col  <- intersect(names(d), c("Estimate","b","Beta","beta"))[1]
    se_col <- intersect(names(d), c("StdError","SE","se","Std. Error"))[1]
    p_col  <- intersect(
      names(d),
      c("Pr(>|t|)","P_value","pval","p","P","P-value")
    )[1]
    
    if (is.na(b_col) || is.na(se_col)) return(list(b=NA_real_, se=NA_real_, p=NA_real_))
    if (nrow(d) < which_exposure) return(list(b=NA_real_, se=NA_real_, p=NA_real_))
    
    out <- list(
      b  = suppressWarnings(as.numeric(d[[b_col]][which_exposure])),
      se = suppressWarnings(as.numeric(d[[se_col]][which_exposure])),
      p  = if (!is.na(p_col)) suppressWarnings(as.numeric(d[[p_col]][which_exposure])) else NA_real_
    )
    
    if (!is.finite(out$b))  out$b  <- NA_real_
    if (!is.finite(out$se)) out$se <- NA_real_
    if (!is.finite(out$p))  out$p  <- NA_real_
    out
  }
  
  # Delta-method SE for product a*b assuming independence
  se_product_delta <- function(a, se_a, b, se_b) {
    sqrt((b^2) * (se_a^2) + (a^2) * (se_b^2))
  }
  
  # --- Storage reactives for Tab7 outputs ---
  tab7_rv <- reactiveValues(
    res_mvmr = data.frame(Message="Run mediation model."),
    res_2smr = data.frame(Message="Run mediation model."),
    interpretation = "Run mediation model to see interpretation."
  )
  # ===============================
  # TAB 7C: has-results flag
  # ===============================
  output$tab7c_has_results <- reactive({
    
    has_mvmr <- !is.null(tab7_rv$res_mvmr) &&
      !"Message" %in% names(tab7_rv$res_mvmr) &&
      nrow(tab7_rv$res_mvmr) > 0
    
    has_2smr <- !is.null(tab7_rv$res_2smr) &&
      !"Message" %in% names(tab7_rv$res_2smr) &&
      nrow(tab7_rv$res_2smr) > 0
    
    has_mvmr && has_2smr
  })
  outputOptions(output, "tab7c_has_results", suspendWhenHidden = FALSE)
  
  # -------------------------------------------------
  # Detect if Outcome is Binary (Tab 7)
  # -------------------------------------------------
  binary_outcome_tab7 <- reactive({
    
    out_ids <- outcome_mod$selected_gcst()
    if (length(out_ids) == 0) return(FALSE)
    
    row <- DATA$merged %>%
      dplyr::filter(`STUDY ACCESSION` == out_ids[1])
    
    if (nrow(row) == 0) return(FALSE)
    
    identical(row$`Data Type`[1], "Binary")
  })
  observeEvent(input$run_mediation_2smr, {
    
    ci_method <- input$tab7_ci_method
    use_delta     <- ci_method %in% c("delta", "both")
    use_bootstrap <- ci_method %in% c("bootstrap", "both")
    
    show_busy_modal("Running mediation...", 
                    "Reusing X→Y (Tab3) and MVMR results when possible…")
    on.exit(removeModal(), add = TRUE)
    
    # ---------------------------------------------------------
    # 0) Extract selections
    # ---------------------------------------------------------
    exp_gcst <- unique(str_trim(as.character(exposure_mod$selected_gcst())))[1]
    out_gcst <- unique(str_trim(as.character(outcome_mod$selected_gcst())))[1]
    med_gcst <- unique(str_trim(as.character(mediator_gcst())))[1]
    
    if (any(is.na(c(exp_gcst, out_gcst, med_gcst))) ||
        any(c(exp_gcst, out_gcst, med_gcst) == "")) {
      
      tab7_rv$res_mvmr <- data.frame(Message="Need Exposure 1, Outcome, and Mediator.")
      tab7_rv$res_2smr <- data.frame(Message="Need Exposure 1, Outcome, and Mediator.")
      tab7_rv$steiger  <- data.frame(Message="Missing required selections.")
      tab7_rv$interpretation <- "Missing required selections."
      return()
    }
    
    # ---------------------------------------------------------
    # 1) TOTAL effect X → Y  (Tab 3)
    # ---------------------------------------------------------
    key_xy <- paste0(exp_gcst, "  →  ", out_gcst)
    res_xy <- mr_rv$results_by_pair[[key_xy]]
    
    if (is.null(res_xy) || !is.null(res_xy$error)) {
      
      tab7_rv$res_mvmr <- data.frame(Message="Run Tab 3 MR first.")
      tab7_rv$res_2smr <- data.frame(Message="Run Tab 3 MR first.")
      tab7_rv$steiger  <- data.frame(Message="No Steiger results.")
      tab7_rv$interpretation <- "Run Tab 3 MR before mediation."
      return()
    }
    
    beta_total <- extract_bse(as_df_or_msg(res_xy$mr_or))
    
    # ---------------------------------------------------------
    # 2) X → M (Tab 5 IVW)
    # ---------------------------------------------------------
    if (is.null(tab5_mr_rv$meta) ||
        !identical(str_trim(tab5_mr_rv$meta$exp2_gcst), str_trim(exp_gcst)) ||
        !identical(str_trim(tab5_mr_rv$meta$out_gcst),  str_trim(med_gcst))) {
      
      tab7_rv$res_mvmr <- data.frame(Message="Run Tab 5 MR for X→M first.")
      tab7_rv$res_2smr <- data.frame(Message="Run Tab 5 MR for X→M first.")
      tab7_rv$steiger  <- data.frame(Message="No Steiger results.")
      tab7_rv$interpretation <- "Tab 5 must match selected mediator."
      return()
    }
    
    beta_xm <- tab5_mr_rv$ivw_summary
    
    # ---------------------------------------------------------
    # 3) DIRECT + M EFFECT (Tab 6 MVMR)
    # ---------------------------------------------------------
    if (is.null(mvmr_rv$res_ivw) ||
        is.data.frame(mvmr_rv$res_ivw) && "Message" %in% names(mvmr_rv$res_ivw)) {
      
      tab7_rv$res_mvmr <- data.frame(Message="Run Tab 6 MVMR first.")
      tab7_rv$res_2smr <- data.frame(Message="Run Tab 6 MVMR first.")
      tab7_rv$steiger  <- data.frame(Message="No Steiger results.")
      tab7_rv$interpretation <- "MVMR results required."
      return()
    }
    
    if (!identical(str_trim(mvmr_rv$used$exp2_gcst), str_trim(med_gcst))) {
      
      tab7_rv$res_mvmr <- data.frame(Message="Mediator differs from last MVMR run.")
      tab7_rv$res_2smr <- data.frame(Message="Rebuild SNP table and rerun MVMR.")
      tab7_rv$steiger  <- data.frame(Message="No Steiger results.")
      tab7_rv$interpretation <- "Mediator mismatch with MVMR."
      return()
    }
    
    beta_direct_x <- get_mvmr_bse(mvmr_rv$res_ivw, 1)
    beta_my_adj   <- get_mvmr_bse(mvmr_rv$res_ivw, 2)
    
    # ---------------------------------------------------------
    # Sanity check
    # ---------------------------------------------------------
    if (!all(is.finite(c(beta_total$b, beta_total$se,
                         beta_xm$b, beta_xm$se,
                         beta_direct_x$b, beta_direct_x$se,
                         beta_my_adj$b, beta_my_adj$se)))) {
      
      tab7_rv$res_mvmr <- data.frame(Message="Could not extract valid effect sizes.")
      tab7_rv$res_2smr <- data.frame(Message="Invalid MR estimates.")
      tab7_rv$steiger  <- data.frame(Message="No Steiger results.")
      tab7_rv$interpretation <- "Extraction error."
      return()
    }
    
    # ---------------------------------------------------------
    # Pathway support
    # ---------------------------------------------------------
    pathway_supported <- (beta_xm$p < 0.05) && (beta_my_adj$p < 0.05)
    
    pathway_note <- if (pathway_supported)
      "Supported (X→M and M→Y|X significant)"
    else
      "⚠️ Weak/unsupported pathway (X→M and/or M→Y|X not significant)"
    
    # ---------------------------------------------------------
    # DIFFERENCE METHOD
    # ---------------------------------------------------------
    ind_diff <- beta_total$b - beta_direct_x$b
    if (use_delta) {
      
      se_diff  <- sqrt(beta_total$se^2 + beta_direct_x$se^2)
      p_diff   <- 2 * pnorm(-abs(ind_diff / se_diff))
      
      ci_diff_l <- ind_diff - 1.96 * se_diff
      ci_diff_u <- ind_diff + 1.96 * se_diff
      
    } else {
      
      p_diff <- NA
      ci_diff_l <- NA
      ci_diff_u <- NA
    }
    
    if (use_bootstrap) {
      
      B <- input$tab7_boot_B
      
      boot_vals <- replicate(B, {
        total_b  <- rnorm(1, beta_total$b, beta_total$se)
        direct_b <- rnorm(1, beta_direct_x$b, beta_direct_x$se)
        total_b - direct_b
      })
      
      # Percentile CI
      ci_boot_l <- quantile(boot_vals, 0.025, na.rm = TRUE)
      ci_boot_u <- quantile(boot_vals, 0.975, na.rm = TRUE)
      
      # Percentile two-sided p-value
      p_boot_diff <- 2 * min(
        mean(boot_vals <= 0, na.rm = TRUE),
        mean(boot_vals >= 0, na.rm = TRUE)
      )
      
    } else {
      p_boot_diff <- NA_real_
    }
    
    prop_diff <- if (beta_total$b != 0) ind_diff / beta_total$b else NA_real_
    
    # Safe CI
    if (use_delta) {
      
      ci_vals <- c(
        ind_diff - 1.96 * se_diff,
        ind_diff + 1.96 * se_diff
      )
      
      ci_diff_l <- min(ci_vals)
      ci_diff_u <- max(ci_vals)
      
    } else {
      
      ci_diff_l <- NA
      ci_diff_u <- NA
    }
    
    
    # ---- Proportion CI (Difference method) ----
    if (use_bootstrap && is.finite(beta_total$b) && beta_total$b != 0) {
      
      boot_prop <- boot_vals / beta_total$b
      
      ci_prop_diff_l <- quantile(boot_prop, 0.025, na.rm = TRUE)
      ci_prop_diff_u <- quantile(boot_prop, 0.975, na.rm = TRUE)
      
    } else if (use_delta && is.finite(beta_total$b) && beta_total$b != 0) {
      
      se_prop_diff <- se_diff / abs(beta_total$b)
      
      ci_prop_vals <- c(
        prop_diff - 1.96 * se_prop_diff,
        prop_diff + 1.96 * se_prop_diff
      )
      
      ci_prop_diff_l <- min(ci_prop_vals)
      ci_prop_diff_u <- max(ci_prop_vals)
      
    } else {
      
      ci_prop_diff_l <- NA_real_
      ci_prop_diff_u <- NA_real_
      
    }
    
    # ---------------------------------------------------------
    # PRODUCT METHOD
    # ---------------------------------------------------------
    ind_prod <- beta_xm$b * beta_my_adj$b
    if (use_delta) {
      
      se_prod  <- se_product_delta(beta_xm$b, beta_xm$se,
                                   beta_my_adj$b, beta_my_adj$se)
      
      p_prod   <- 2 * pnorm(-abs(ind_prod / se_prod))
      
      ci_prod_l <- ind_prod - 1.96 * se_prod
      ci_prod_u <- ind_prod + 1.96 * se_prod
      
    } else {
      
      p_prod <- NA
      ci_prod_l <- NA
      ci_prod_u <- NA
    }
    if (use_bootstrap) {
      
      B <- input$tab7_boot_B
      
      boot_prod <- replicate(B, {
        a <- rnorm(1, beta_xm$b, beta_xm$se)
        b <- rnorm(1, beta_my_adj$b, beta_my_adj$se)
        a * b
      })
      
      # Percentile CI
      ci_prod_boot_l <- quantile(boot_prod, 0.025, na.rm = TRUE)
      ci_prod_boot_u <- quantile(boot_prod, 0.975, na.rm = TRUE)
      
      # Percentile two-sided p-value
      p_boot_prod <- 2 * min(
        mean(boot_prod <= 0, na.rm = TRUE),
        mean(boot_prod >= 0, na.rm = TRUE)
      )
      
    } else {
      p_boot_prod <- NA_real_
    }
    
    prop_prod <- if (beta_total$b != 0) ind_prod / beta_total$b else NA_real_
    
    
    # ---- Proportion CI (Product method) ----
    if (use_bootstrap && is.finite(beta_total$b) && beta_total$b != 0) {
      
      boot_prop_prod <- boot_prod / beta_total$b
      
      ci_prop_prod_l <- quantile(boot_prop_prod, 0.025, na.rm = TRUE)
      ci_prop_prod_u <- quantile(boot_prop_prod, 0.975, na.rm = TRUE)
      
    } else if (use_delta && is.finite(beta_total$b) && beta_total$b != 0) {
      
      se_prop_prod <- se_prod / abs(beta_total$b)
      
      ci_prop_vals2 <- c(
        prop_prod - 1.96 * se_prop_prod,
        prop_prod + 1.96 * se_prop_prod
      )
      
      ci_prop_prod_l <- min(ci_prop_vals2)
      ci_prop_prod_u <- max(ci_prop_vals2)
      
    } else {
      
      ci_prop_prod_l <- NA_real_
      ci_prop_prod_u <- NA_real_
      
    }
    
    # ---------------------------------------------------------
    # Interpretation
    # ---------------------------------------------------------
    sig_diff <- if (ci_method == "delta") {
      p_diff < 0.05
    } else if (ci_method == "bootstrap") {
      p_boot_diff < 0.05
    } else {
      (p_diff < 0.05) || (p_boot_diff < 0.05)
    }
    
    sig_prod <- if (ci_method == "delta") {
      p_prod < 0.05
    } else if (ci_method == "bootstrap") {
      p_boot_prod < 0.05
    } else {
      (p_prod < 0.05) || (p_boot_prod < 0.05)
    }
    
    tab7_rv$interpretation <- if (!pathway_supported) {
      "⚠️ Pathway not fully supported."
    } else if (sig_diff && sig_prod) {
      "✅ Strong mediation evidence."
    } else if (sig_diff || sig_prod) {
      "ℹ️ Partial mediation evidence."
    } else {
      "❌ No significant mediation detected."
    }
    
    # -----------------------------
    # Format Indirect_DIFF
    # -----------------------------
    if (ci_method == "delta") {
      
      Indirect_DIFF_txt <- paste0(
        round(ind_diff,6),
        "\n(", round(ci_diff_l,6),
        ", ", round(ci_diff_u,6), ")"
      )
      
    } else if (ci_method == "bootstrap") {
      
      Indirect_DIFF_txt <- paste0(
        round(ind_diff,6),
        "\nBootstrap CI: (",
        round(ci_boot_l,6), ", ",
        round(ci_boot_u,6), ")"
      )
      
    } else {  # both
      
      Indirect_DIFF_txt <- paste0(
        round(ind_diff,6),
        "\nDelta: (", round(ci_diff_l,6), ", ", round(ci_diff_u,6), ")",
        "\nBootstrap: (", round(ci_boot_l,6), ", ", round(ci_boot_u,6), ")"
      )
    }
    # ---------------------------------------------------------
    # SAVE RESULTS
    # ---------------------------------------------------------
    tab7_rv$res_mvmr <- data.frame(
      Exposure = get_trait_label(exp_gcst, DATA$merged),
      Mediator = get_trait_label(med_gcst, DATA$merged),
      Outcome  = get_trait_label(out_gcst, DATA$merged),
      Pathway_support = pathway_note,
      Total_X_to_Y = round(beta_total$b, 6),
      Direct_X_to_Y_adjM = round(beta_direct_x$b, 6),
      Indirect_DIFF = Indirect_DIFF_txt,
      P_DIFF = if (ci_method == "delta") {
        signif(p_diff,4)
      } else if (ci_method == "bootstrap") {
        signif(p_boot_diff,4)
      } else {
        paste0(
          "Delta: ", signif(p_diff,4),
          "\nBoot: ", signif(p_boot_diff,4)
        )
      },
      Prop_DIFF = paste0(round(prop_diff,6),
                         "\n(", round(ci_prop_diff_l,6),
                         ", ", round(ci_prop_diff_u,6), ")"),
      stringsAsFactors = FALSE
    )
    # -----------------------------
    # Format Indirect_PROD
    # -----------------------------
    if (ci_method == "delta") {
      
      Indirect_PROD_txt <- paste0(
        round(ind_prod,6),
        "\n(", round(ci_prod_l,6),
        ", ", round(ci_prod_u,6), ")"
      )
      
    } else if (ci_method == "bootstrap") {
      
      Indirect_PROD_txt <- paste0(
        round(ind_prod,6),
        "\nBootstrap CI: (",
        round(ci_prod_boot_l,6), ", ",
        round(ci_prod_boot_u,6), ")"
      )
      
    } else {
      
      Indirect_PROD_txt <- paste0(
        round(ind_prod,6),
        "\nDelta: (", round(ci_prod_l,6), ", ", round(ci_prod_u,6), ")",
        "\nBootstrap: (", round(ci_prod_boot_l,6), ", ", round(ci_prod_boot_u,6), ")"
      )
    }
    
    tab7_rv$res_2smr <- data.frame(
      Exposure = get_trait_label(exp_gcst, DATA$merged),
      Mediator = get_trait_label(med_gcst, DATA$merged),
      Outcome  = get_trait_label(out_gcst, DATA$merged),
      Pathway_support = pathway_note,
      X_to_M = round(beta_xm$b,6),
      M_to_Y_adjX = round(beta_my_adj$b,6),
      Indirect_PROD = Indirect_PROD_txt,
      P_PROD = if (ci_method == "delta") {
        signif(p_prod,4)
      } else if (ci_method == "bootstrap") {
        signif(p_boot_prod,4)
      } else {
        paste0(
          "Delta: ", signif(p_prod,4),
          "\nBoot: ", signif(p_boot_prod,4)
        )
      },
      Prop_PROD = paste0(round(prop_prod,6),
                         "\n(", round(ci_prop_prod_l,6),
                         ", ", round(ci_prop_prod_u,6), ")"),
      stringsAsFactors = FALSE
    )
    # ---------------------------------------------------------
    # Data Type–Dependent Mediation Warning
    # ---------------------------------------------------------
    
    out_row <- DATA$merged %>%
      dplyr::filter(`STUDY ACCESSION` == out_gcst)
    
    med_row <- DATA$merged %>%
      dplyr::filter(`STUDY ACCESSION` == med_gcst)
    
    if (nrow(out_row) > 0 && nrow(med_row) > 0) {
      
      outcome_type  <- out_row$`Data Type`[1]
      mediator_type <- med_row$`Data Type`[1]
      
      msg <- NULL
      
      # ---------------------------------------------------------
      # Case 1: Outcome binary & Mediator continuous
      # ---------------------------------------------------------
      if (outcome_type == "Binary" && mediator_type == "Continuous") {
        
        msg <- paste0(
          "⚠ Outcome is binary.\n\n",
          "• MVMR Difference Method may be inappropriate due to non-collapsibility of log-odds.\n",
          "• Product method is more appropriate in this setting.\n",
          "• Indirect effect is on log-odds scale.\n",
          "• Use exp(beta) for Odds Ratio interpretation."
        )
        
        # ---------------------------------------------------------
        # Case 2: Mediator binary & Outcome continuous
        # ---------------------------------------------------------
      } else if (mediator_type == "Binary" && outcome_type == "Continuous") {
        
        msg <- paste0(
          "⚠ Mediator is binary.\n\n",
          "• X → M path is estimated on log-odds scale.\n",
          "• Indirect effect involves log-odds for mediator.\n",
          "• Consider interpreting mediator effects using Odds Ratios."
        )
        
        # ---------------------------------------------------------
        # Case 3: Both mediator and outcome binary
        # ---------------------------------------------------------
      } else if (outcome_type == "Binary" && mediator_type == "Binary") {
        
        msg <- paste0(
          "⚠ Both Mediator and Outcome are binary.\n\n",
          "• MVMR Difference Method is not ideal due to log-odds non-collapsibility.\n",
          "• Product method may be preferable, but interpretation remains complex.\n",
          "• Indirect effect is on log-odds scale.\n",
          "• Interpretation requires caution due to non-linear modeling structure."
        )
        
      }
      
      if (!is.null(msg)) {
        
        showNotification(
          msg,
          type = "warning",
          duration = 15
        )
        
      }
    }
  })
  
  # -------------------------------
  # TAB 7C: Mediation results download
  # -------------------------------
  output$dl_tab7c_mediation_xlsx <- downloadHandler(
    filename = function() {
      paste0(
        "TAB7_Mediation_Results_",
        Sys.Date(),
        ".xlsx"
      )
    },
    content = function(file) {
      
      # Create workbook
      wb <- openxlsx::createWorkbook()
      
      # ---------------------------
      # Sheet 1: MVMR Difference
      # ---------------------------
      df_mvmr <- tab7_rv$res_mvmr
      
      if (is.null(df_mvmr) || nrow(df_mvmr) == 0) {
        df_mvmr <- data.frame(Message = "No MVMR mediation results.")
      }
      
      openxlsx::addWorksheet(wb, "MVMR_Difference")
      openxlsx::writeData(wb, "MVMR_Difference", df_mvmr)
      
      
      # ---------------------------
      # Sheet 2: TwoSMR Product
      # ---------------------------
      df_2smr <- tab7_rv$res_2smr
      
      if (is.null(df_2smr) || nrow(df_2smr) == 0) {
        df_2smr <- data.frame(Message = "No TwoSMR mediation results.")
      }
      
      openxlsx::addWorksheet(wb, "TwoSMR_Product")
      openxlsx::writeData(wb, "TwoSMR_Product", df_2smr)
      
      # Save file
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  output$tbl_mediation_mvmr <- renderDT({
    
    df <- tab7_rv$res_mvmr
    
    if (binary_outcome_tab7()) {
      
      return(
        datatable(
          data.frame(
            Message = "⚠ Difference Method not appropriate for binary outcomes due to log-odds non-collapsibility. Use Product Method instead."
          ),
          rownames = FALSE,
          options = list(dom = "t")
        )
      )
    }
    
    if (is.null(df) || "Message" %in% names(df)) {
      return(datatable(df, rownames = FALSE, options = list(dom="t")))
    }
    
    colnames(df) <- c(
      "Exposure",
      "Mediator",
      "Outcome",
      "Pathway Support",
      "Total Effect (X → Y)",
      "Direct Effect (X → Y | M)",
      "Indirect Effect (Difference)",
      "P-value (Difference)",
      "Mediated Proportion"
    )
    
    datatable(
      df,
      rownames = FALSE,
      options = list(scrollX = TRUE, pageLength = 5, dom = "tip")
    )
  })
  
  output$tbl_mediation_2smr <- renderDT({
    
    df <- tab7_rv$res_2smr
    
    if (is.null(df) || "Message" %in% names(df)) {
      return(datatable(df, rownames = FALSE, options = list(dom="t")))
    }
    
    colnames(df) <- c(
      "Exposure",
      "Mediator",
      "Outcome",
      "Pathway Support",
      "X → M Effect",
      "M → Y Effect (adjusted for X)",
      "Indirect Effect (Product)",
      "P-value (Product)",
      "Mediated Proportion"
    )
    
    datatable(
      df,
      rownames = FALSE,
      options = list(scrollX = TRUE, pageLength = 5, dom = "tip")
    )
  })
  
  # --------------------------------------------------------------------
  # TAB 6: Multivariable MR (MVMR)
  # --------------------------------------------------------------------
  mvmr_rv <- reactiveValues(
    status = "Build SNP table, then click Run MVMR.",
    res_ivw = NULL,
    res_strength = NULL,
    res_pleio = NULL,
    used = list()
  )
  output$mvmr_status <- renderText({ mvmr_rv$status })
  
  # ===============================
  # TAB 6: has-results flag
  # ===============================
  output$mvmr_has_results <- reactive({
    !is.null(mvmr_rv$res_ivw) &&
      !"Message" %in% names(mvmr_rv$res_ivw)
  })
  outputOptions(output, "mvmr_has_results", suspendWhenHidden = FALSE)
  
  observeEvent(input$run_mvmr, {
    show_busy_modal("Performing MVMR...", "Running MVMR now… (No OpenGWAS)")
    on.exit({ removeModal() }, add = TRUE)
    
    if (!requireNamespace("MVMR", quietly = TRUE)) {
      mvmr_rv$status <- "Package MVMR is not installed. Install with: remotes::install_github('WSpiller/MVMR')"
      return()
    }
    
    df <- tab5_rv$df
    if (is.null(df) || nrow(df) == 0) {
      mvmr_rv$status <- "Cannot run MVMR: Build Tab 5 SNP table first."
      mvmr_rv$res_ivw <- mvmr_rv$res_strength <- mvmr_rv$res_pleio <- NULL
      return()
    }
    if ("Message" %in% names(df)) {
      mvmr_rv$status <- paste0("Cannot run MVMR: ", (df$Message %||% "Build Tab 5 SNP table first.")[1])
      mvmr_rv$res_ivw <- mvmr_rv$res_strength <- mvmr_rv$res_pleio <- NULL
      return()
    }
    
    # Record which GCSTs were used for THIS MVMR run
    mvmr_rv$used <- list(
      exp1_gcst = unique(str_trim(as.character(exposure_mod$selected_gcst())))[1],
      exp2_gcst = unique(str_trim(as.character(exp2_for_mvmr())))[1],
      out_gcst  = unique(str_trim(as.character(outcome_mod$selected_gcst())))[1]
    )
    
    if (nrow(df) < 3) {
      mvmr_rv$status <- paste0(
        "⚠️ Running MVMR with very few SNPs (n=", nrow(df), "). ",
        "Results may be unstable/under-identified. Diagnostics may fail, but outputs will be shown if available."
      )
    } else {
      mvmr_rv$status <- "Performing MVMR..."
    }
    df2 <- df %>%
      dplyr::mutate(
        Outcome_beta   = as.numeric(Outcome_beta),
        Outcome_se     = as.numeric(Outcome_se),
        Exposure1_beta = as.numeric(Exposure1_beta),
        Exposure1_se   = as.numeric(Exposure1_se),
        Exposure2_beta = as.numeric(Exposure2_beta),
        Exposure2_se   = as.numeric(Exposure2_se)
      ) %>%
      dplyr::filter(
        is.finite(Outcome_beta), is.finite(Outcome_se),
        is.finite(Exposure1_beta), is.finite(Exposure1_se),
        is.finite(Exposure2_beta), is.finite(Exposure2_se),
        !is.na(SNP), SNP != ""
      )
    
    if (nrow(df2) < 3) {
      mvmr_rv$status <- paste0("Cannot run MVMR: only ", nrow(df2), " valid SNP rows after filtering.")
      mvmr_rv$res_ivw <- mvmr_rv$res_strength <- mvmr_rv$res_pleio <- NULL
      return()
    }
    
    F.data <- tryCatch(
      MVMR::format_mvmr(
        BXGs   = as.matrix(df2[, c("Exposure1_beta", "Exposure2_beta"), drop = FALSE]),
        BYG    = df2$Outcome_beta,
        seBXGs = as.matrix(df2[, c("Exposure1_se", "Exposure2_se"), drop = FALSE]),
        seBYG  = df2$Outcome_se,
        RSID   = df2$SNP
      ),
      error = function(e) e
    )
    
    if (inherits(F.data, "error")) {
      mvmr_rv$status <- paste0("format_mvmr failed: ", conditionMessage(F.data))
      mvmr_rv$res_ivw <- NULL
      mvmr_rv$res_strength <- NULL
      mvmr_rv$res_pleio <- NULL
      return()
    }
    
    sres <- tryCatch(MVMR::strength_mvmr(r_input = F.data, gencov = 0), error = function(e) e)
    pres <- tryCatch(MVMR::pleiotropy_mvmr(r_input = F.data, gencov = 0), error = function(e) e)
    ivw  <- tryCatch(MVMR::ivw_mvmr(r_input = F.data), error = function(e) e)
    
    # -------------------------
    # Instrument Strength (Conditional F only)
    # -------------------------
    
    if (inherits(sres, "error")) {
      
      strength_processed <- data.frame(
        Message = paste0("strength_mvmr failed: ", conditionMessage(sres)),
        stringsAsFactors = FALSE
      )
      
    } else {
      
      F_vals <- as.numeric(sres)
      names_vals <- names(sres)
      if (is.null(names_vals)) {
        names_vals <- paste0("Exposure", seq_along(F_vals))
      }
      
      strength_processed <- data.frame(
        Exposure = names_vals,
        `Conditional F-statistic (Sanderson-Windmeijer)` = round(F_vals, 4),
        stringsAsFactors = FALSE,
        check.names = FALSE
        
      )
      colnames(strength_processed)[2] <- 
        "Conditional F-statistic (Sanderson-Windmeijer)"
      
      strength_processed$`Instrument Strength` <- ifelse(
        F_vals > 10,
        "Strong (F > 10)",
        "⚠ Weak (F ≤ 10)"
      )
    }
    
    # -------------------------
    # Pleiotropy (ROBUST)  
    # -------------------------
    if (inherits(pres, "error")) {
      pleio_processed <- data.frame(
        Message = paste0("pleiotropy_mvmr failed: ", conditionMessage(pres)),
        stringsAsFactors = FALSE
      )
    } else {
      q_val <- NA_real_
      p_val <- NA_real_
      
      if (is.list(pres)) {
        q_val <- pres$Qstat %||% pres$Q %||% pres$Q_stat %||% pres$q %||% NA
        p_val <- pres$Q_pval %||% pres$Qpval %||% pres$pval %||% pres$p %||% NA
        
      } else if (is.matrix(pres) || is.data.frame(pres)) {
        tmp <- as.data.frame(pres)
        q_col <- intersect(names(tmp), c("Qstat","Q","Q_stat","Q-statistic","Q_statistic"))[1]
        p_col <- intersect(names(tmp), c("Q_pval","Qpval","pval","p.value","p_value","p","Pr(>|t|)","Pr(>|z|)"))[1]
        if (!is.na(q_col)) q_val <- suppressWarnings(as.numeric(tmp[[q_col]][1]))
        if (!is.na(p_col)) p_val <- suppressWarnings(as.numeric(tmp[[p_col]][1]))
      }
      
      pleio_processed <- data.frame(
        Method = "IVW",
        Qstat  = suppressWarnings(as.numeric(q_val)),
        Q_pval = suppressWarnings(as.numeric(p_val)),
        stringsAsFactors = FALSE
      )
      
      pleio_processed <- add_significance_column(pleio_processed)
      pleio_processed <- round_df_display(pleio_processed, digits = 4)
    }
    
    # -------------------------
    # IVW (ROBUST)
    # -------------------------
    if (inherits(ivw, "error")) {
      ivw_processed <- data.frame(
        Message = paste0("ivw_mvmr failed: ", conditionMessage(ivw)),
        stringsAsFactors = FALSE
      )
    } else {
      
      if (is.data.frame(ivw)) {
        ivw_processed <- ivw
      } else if (is.matrix(ivw)) {
        ivw_processed <- as.data.frame(ivw)
      } else if (is.list(ivw)) {
        ivw_processed <- as.data.frame(ivw)
      } else {
        ivw_processed <- data.frame(Estimate = as.numeric(ivw), stringsAsFactors = FALSE)
      }
      
      if (!"Exposure ID" %in% names(ivw_processed)) {
        if (!is.null(rownames(ivw_processed)) && any(rownames(ivw_processed) != "")) {
          ivw_processed <- cbind(`Exposure ID` = rownames(ivw_processed), ivw_processed)
          rownames(ivw_processed) <- NULL
        } else {
          ivw_processed <- cbind(
            `Exposure ID` = paste0("Exposure", seq_len(nrow(ivw_processed))),
            ivw_processed
          )
        }
      }
      
      ivw_processed <- round_df_display(ivw_processed, digits = 6)
    }
    
    # Save results
    mvmr_rv$res_strength <- strength_processed
    mvmr_rv$res_pleio    <- pleio_processed
    mvmr_rv$res_ivw      <- ivw_processed
    mvmr_rv$status <- "Finished (MVMR complete)"
    
    exp1 <- unique(str_trim(as.character(exposure_mod$selected_gcst())))
    exp2 <- unique(str_trim(as.character(mediator_gcst())))
    out1 <- unique(str_trim(as.character(outcome_mod$selected_gcst())))
    exp1_gcst <- exp1[1]
    exp2_gcst <- exp2[1]
    out_gcst  <- out1[1]
    mvmr_rv$used <- list(
      exp1_gcst = exp1_gcst,
      exp2_gcst = exp2_gcst,
      out_gcst  = out_gcst,
      n_snps    = nrow(df)
    )
    
    mvmr_rv$status <- paste0(
      (if (nrow(df) < 3) "Finished (with warning) — " else "Finished — "),
      "No-OpenGWAS MVMR. SNPs used: ", nrow(df)
    )
    # --------------------------------------------
    # Binary Outcome Warning After MVMR
    # --------------------------------------------
    
    out_ids <- outcome_mod$selected_gcst()
    
    if (length(out_ids) > 0) {
      
      out_id <- out_ids[1]
      
      row <- DATA$merged %>%
        dplyr::filter(`STUDY ACCESSION` == out_id)
      
      if (nrow(row) > 0) {
        
        outcome_type <- row$`Data Type`[1]
        
        if (!is.na(outcome_type) && outcome_type == "Binary") {
          
          showNotification(
            paste0(
              "⚠ Outcome ", out_id, " (",
              row$MAPPED_TRAIT[1],
              ") is binary. Direct effects from MVMR are in log-odds scale. ",
              "Use exp(beta) for Odds Ratio interpretation."
            ),
            type = "warning",
            duration = 10
          )
          
        }
        
      }
      
    }
  })
  
  output$tbl_mvmr_pleio <- renderDT({
    .safe_dt(mvmr_rv$res_pleio, msg = "No pleiotropy results yet.")
  })
  
  output$tbl_mvmr_ivw <- renderDT({
    
    df <- mvmr_rv$res_ivw
    
    if (is.null(df)) {
      return(datatable(
        data.frame(Message = "No results yet."),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }
    
    # -------------------------------------------------
    # Get GCST IDs
    # -------------------------------------------------
    exp1_gcst <- unique(str_trim(as.character(exposure_mod$selected_gcst())))[1]
    med_gcst  <- unique(str_trim(as.character(exp2_for_mvmr())))[1]
    out_gcst  <- unique(str_trim(as.character(outcome_mod$selected_gcst())))[1]
    
    # -------------------------------------------------
    # Get trait names from DATA$merged
    # -------------------------------------------------
    get_trait <- function(gcst) {
      x <- DATA$merged %>%
        dplyr::filter(`STUDY ACCESSION` == gcst) %>%
        dplyr::pull(MAPPED_TRAIT)
      
      if (length(x) == 0 || is.na(x[1])) return(gcst)
      x[1]
    }
    
    exp1_trait <- get_trait(exp1_gcst)
    med_trait  <- get_trait(med_gcst)
    out_trait  <- get_trait(out_gcst)
    
    # -------------------------------------------------
    # Convert rownames if needed
    # -------------------------------------------------
    if (!"Exposure ID" %in% names(df) && !is.null(rownames(df))) {
      df <- cbind(`Exposure ID` = rownames(df), df)
      rownames(df) <- NULL
    }
    
    # Detect columns
    b_col  <- intersect(names(df), c("Estimate","b","Beta","beta"))[1]
    se_col <- intersect(names(df), c("Std. Error","StdError","SE","se"))[1]
    # ---- Auto-detect p-value column ----
    # ---- Robust p-value column detection ----
    p_col <- names(df)[
      grepl("p", names(df), ignore.case = TRUE) &
        !grepl("exp", names(df), ignore.case = TRUE) &
        !grepl("prop", names(df), ignore.case = TRUE)
    ][1]
    
    # Force numeric conversion (important)
    if (!is.na(p_col)) {
      df[[p_col]] <- suppressWarnings(as.numeric(df[[p_col]]))
    }
    t_col  <- intersect(names(df), c("t value","t","t_value"))[1]
    
    if (is.na(b_col) || is.na(se_col)) {
      return(datatable(df, rownames = FALSE))
    }
    
    # -------------------------------------------------
    # Dynamic effect labels
    # -------------------------------------------------
    df$Effect <- ifelse(
      tolower(df$`Exposure ID`) %in% c("exposure1","x","exp1"),
      paste0("Direct effect of ", exp1_trait,
             " on ", out_trait,
             " (adjusted for ", med_trait, ")"),
      
      paste0("Effect of ", med_trait,
             " on ", out_trait,
             " (adjusted for ", exp1_trait, ")")
    )
    
    # -------------------------------------------------
    # Clean table
    # -------------------------------------------------
    clean_df <- data.frame(
      `Causal Effect` = df$Effect,
      `Beta (Estimate)` = round(as.numeric(df[[b_col]]), 6),
      `Standard Error (SE)` = round(as.numeric(df[[se_col]]), 6),
      `t-statistic` = if (!is.na(t_col)) round(as.numeric(df[[t_col]]), 4) else NA,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    # Always create numeric P column safely
    if (!is.na(p_col)) {
      clean_df$`P value` <- as.numeric(df[[p_col]])
    } else {
      clean_df$`P value` <- NA_real_
    }
    
    # Add Significance column
    p_vals <- clean_df$`P value`
    
    clean_df$Significance <- ifelse(
      is.finite(p_vals) & p_vals < 0.001, "*** (p<0.001)",
      ifelse(
        is.finite(p_vals) & p_vals < 0.01, "** (p<0.01)",
        ifelse(
          is.finite(p_vals) & p_vals < 0.05, "* (p<0.05)",
          ifelse(
            is.finite(p_vals) & p_vals < 0.1, ". (p<0.1)",
            "Not significant"
          )
        )
      )
    )
    
    colnames(clean_df) <- c(
      "Causal Effect",
      "Beta (Estimate)",
      "Standard Error (SE)",
      "t-statistic",
      "P value",
      "Significance"
    )
    
    datatable(
      clean_df,
      rownames = FALSE,
      options = list(pageLength = 5, scrollX = TRUE, dom = "tip")
    )
  })
  
  output$tbl_mvmr_strength <- renderDT({
    
    req(mvmr_rv$res_strength)
    
    df <- mvmr_rv$res_strength
    
    # If error message
    if ("Message" %in% names(df)) {
      return(datatable(df, rownames = FALSE, options = list(dom = "t")))
    }
    
    # --------------------------------------
    # Get GCSTs used in THIS MVMR run
    # --------------------------------------
    exp1_gcst <- mvmr_rv$used$exp1_gcst
    med_gcst  <- mvmr_rv$used$exp2_gcst
    
    # --------------------------------------
    # Get trait names dynamically
    # --------------------------------------
    get_trait <- function(gcst) {
      x <- DATA$merged %>%
        dplyr::filter(`STUDY ACCESSION` == gcst) %>%
        dplyr::pull(MAPPED_TRAIT)
      
      if (length(x) == 0 || is.na(x[1])) return(gcst)
      x[1]
    }
    
    exp1_trait <- get_trait(exp1_gcst)
    med_trait  <- get_trait(med_gcst)
    
    # --------------------------------------
    # Extract F statistics safely
    # --------------------------------------
    f_vals <- suppressWarnings(as.numeric(df[[2]]))
    if (length(f_vals) < 2) {
      return(datatable(
        data.frame(Message = "Strength results did not contain two F-statistics."),
        rownames = FALSE, options = list(dom = "t")
      ))
    }
    
    # Build EXACTLY 3 columns BEFORE datatable()
    clean_df <- data.frame(
      Effect = c(
        paste0("Conditional instruments for ", exp1_trait),
        paste0("Conditional instruments for ", med_trait)
      ),
      `Conditional F-statistic` = round(f_vals[1:2], 3),
      `Instrument Strength` = ifelse(f_vals[1:2] > 10, "Strong (F > 10)", "⚠ Weak (F ≤ 10)"),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    #  Ensure EXACTLY these 3 columns (and in the right order)
    clean_df <- clean_df[, c("Effect", "Conditional F-statistic", "Instrument Strength"), drop = FALSE]
    
    datatable(
      clean_df,
      rownames = FALSE,
      options = list(
        pageLength = 5,
        scrollX = TRUE,
        dom = "tip"
      )
    )
  })
  
  output$tbl_mvmr_pleio <- renderDT({
    
    req(mvmr_rv$res_pleio)
    
    df <- mvmr_rv$res_pleio
    
    if ("Message" %in% names(df)) {
      return(datatable(df, options=list(dom="t"), rownames=FALSE))
    }
    
    df$Qstat <- round(df$Qstat, 3)
    df$Q_pval <- signif(df$Q_pval, 4)
    
    df$Interpretation <- ifelse(
      df$Q_pval < 0.05,
      "Evidence of excess heterogeneity (possible pleiotropy)",
      "No strong evidence of heterogeneity"
    )
    
    clean_df <- df[, c("Method","Qstat","Q_pval","Interpretation")]
    
    colnames(clean_df) <- c(
      "Method",
      "Q-statistic",
      "Q p-value",
      "Interpretation"
    )
    
    datatable(
      clean_df,
      rownames = FALSE,
      options = list(pageLength = 5, scrollX = TRUE, dom = "tip")
    )
  })
  
  # --------------------------------------------------------------------
  # Info modals for MR tables (mr results, heterogeneity, pleiotropy, direction)
  # The body of these showModal calls is unchanged from app.R; include if needed.
  # --------------------------------------------------------------------
  
  observeEvent(input$info_mr_results, {
    showModal(modalDialog(
      title = "MR Results table columns",
      HTML(paste0(
        "<table style='width:100%; border-collapse:collapse;'>",
        "<thead>",
        "<tr>",
        "<th style='text-align:left; padding:4px;'>Column name</th>",
        "<th style='text-align:left; padding:4px;'>Description</th>",
        "</tr>",
        "</thead>",
        "<tbody>",
        
        "<tr><td style='padding:4px;'><strong>Exposure (Trait | GCST ID)</strong></td>",
        "<td style='padding:4px;'>Unique identifier of the exposure GWAS dataset corresponding to the exposure trait used in the Mendelian Randomization analysis.</td></tr>",
        
        "<tr><td style='padding:4px;'><strong>Outcome (Trait | GCST ID)</strong></td>",
        "<td style='padding:4px;'>Unique identifier of the outcome GWAS dataset corresponding to the outcome trait analyzed in the Mendelian Randomization analysis.</td></tr>",
        
        "<tr><td style='padding:4px;'><strong>MR Method</strong></td>",
        "<td style='padding:4px;'>Mendelian Randomization estimation method used to compute the causal effect (e.g., IVW, MR&#8209;Egger, Weighted Median, Mode-based estimator).</td></tr>",
        
        "<tr><td style='padding:4px;'><strong>Number of SNPs</strong></td>",
        "<td style='padding:4px;'>Number of single nucleotide polymorphisms (SNPs) used as instrumental variables after filtering and harmonization.</td></tr>",
        
        "<tr><td style='padding:4px;'><strong>Beta (Effect)</strong></td>",
        "<td style='padding:4px;'>stimated causal effect size of the exposure on the outcome per unit increase in the exposure.For continuous outcomes, this represents the linear change in the outcome.For binary outcomes, this represents the change in the log-odds of the outcome (log odds ratio).</td></tr>",
        
        "<tr><td style='padding:4px;'><strong>Standard Error</strong></td>",
        "<td style='padding:4px;'>Standard error of the estimated causal effect, representing statistical uncertainty around the effect estimate.</td></tr>",
        
        "<tr><td style='padding:4px;'><strong>Odds Ratio</strong></td>",
        "<td style='padding:4px;'>Odds ratio corresponding to the causal effect estimate. Reported when the outcome is binary and interpreted as the change in odds of the outcome per unit increase in the exposure.</td></tr>",
       
        "<tr><td style='padding:4px;'><strong>P-value</strong></td>",
        "<td style='padding:4px;'>P-value testing the null hypothesis that the causal effect estimate is equal to zero.</td></tr>",
        
        "<tr><td style='padding:4px;'><strong>Significance</strong></td>",
        "<td style='padding:4px;'>Indicator of statistical significance based on a predefined threshold, reported as significant (***p < 0.001, **p < 0.01, *p < 0.05, and .p < 0.1) and not significant (p ≥ 0.05).</td></tr>",
        
        "</tbody>",
        "</table>"
      )),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  observeEvent(input$info_het, {
    showModal(modalDialog(
      title = "Heterogeneity statistics table columns",
      HTML(paste0(
        "<table style='width:100%; border-collapse:collapse;'>",
        "<thead><tr><th style='text-align:left; padding:4px;'>Column name</th><th style='text-align:left; padding:4px;'>Description</th></tr></thead><tbody>",
        "<tr><td style='padding:4px;'><strong>Exposure (Trait | GCST ID)</strong></td><td style='padding:4px;'>Unique identifier of the exposure GWAS dataset used in the Mendelian Randomization analysis.</td></tr>",
        "<tr><td style='padding:4px;'><strong>Outcome (Trait | GCST ID)</strong></td><td style='padding:4px;'>Unique identifier of the outcome GWAS dataset used in the analysis.</td></tr>",
        "<tr><td style='padding:4px;'><strong>MR Method</strong></td><td style='padding:4px;'>Mendelian Randomization estimation method applied (e.g., IVW, MR&#8209;Egger, Weighted Median, Mode‑based estimator).</td></tr>",
        "<tr><td style='padding:4px;'><strong>Cochran's Q Statistic</strong></td><td style='padding:4px;'>Cochran’s Q statistic used to assess heterogeneity among SNP‑specific causal estimates.</td></tr>",
        "<tr><td style='padding:4px;'><strong>Degrees of Freedom</strong></td><td style='padding:4px;'>Degrees of freedom for the heterogeneity test, typically equal to nsnp − 1.</td></tr>",
        "<tr><td style='padding:4px;'><strong>Q-test P-value</strong></td><td style='padding:4px;'>P‑value associated with Cochran’s Q statistic, testing the null hypothesis of no heterogeneity across SNPs.</td></tr>",
        "<tr><td style='padding:4px;'><strong>Heterogeneity Significance</strong></td><td style='padding:4px;'>Indicator of statistically significant heterogeneity based on the Q‑test. Significant heterogeneity suggests potential pleiotropy, invalid instruments, or model misspecification.Reported as significant (***p < 0.001, **p < 0.01, *p < 0.05, and .p < 0.1) and not significant (p ≥ 0.05).</td></tr>",
        "</tbody></table>"
      )),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  observeEvent(input$info_pleio, {
    showModal(modalDialog(
      title = "Horizontal pleiotropy table columns",
      HTML(paste0(
        "<table style='width:100%; border-collapse:collapse;'>",
        "<thead><tr><th style='text-align:left; padding:4px;'>Column name</th><th style='text-align:left; padding:4px;'>Description</th></tr></thead><tbody>",
        "<tr><td style='padding:4px;'><strong>Exposure (Trait | GCST ID)</strong></td><td style='padding:4px;'>Unique identifier of the exposure GWAS dataset used in the Mendelian Randomization analysis.</td></tr>",
        "<tr><td style='padding:4px;'><strong>Outcome (Trait | GCST ID)</strong></td><td style='padding:4px;'>Unique identifier of the outcome GWAS dataset used in the analysis.</td></tr>",
        "<tr><td style='padding:4px;'><strong>MR-Egger Intercept</strong></td><td style='padding:4px;'>Estimated intercept from MR‑Egger regression, representing the average direct (horizontal pleiotropic) effect of the genetic instruments on the outcome.</td></tr>",
        "<tr><td style='padding:4px;'><strong>Standard Error</strong></td><td style='padding:4px;'>Standard error of the MR‑Egger intercept estimate, reflecting uncertainty in the estimated pleiotropic effect.</td></tr>",
        "<tr><td style='padding:4px;'><strong>P-value</strong></td><td style='padding:4px;'>P‑value testing the null hypothesis that the MR‑Egger intercept equals zero (H&#8320;: intercept = 0).</td></tr>",
        "<tr><td style='padding:4px;'><strong>Pleiotropy Significance</strong></td><td style='padding:4px;'>Indicator of statistically significant horizontal pleiotropy, typically defined as pval &lt; 0.05. A significant value suggests evidence of directional (unbalanced) pleiotropy.</td></tr>",
        "</tbody></table>"
      )),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  observeEvent(input$info_dir, {
    showModal(modalDialog(
      title = "Causal direction test table columns",
      HTML(paste0(
        "<table style='width:100%; border-collapse:collapse;'>",
        "<thead><tr><th style='text-align:left; padding:4px;'>Column name</th><th style='text-align:left; padding:4px;'>Description</th></tr></thead><tbody>",
        "<tr><td style='padding:4px;'><strong>Exposure (Trait | GCST ID)</strong></td><td style='padding:4px;'>Unique identifier of the exposure GWAS dataset used in the Mendelian Randomization analysis.</td></tr>",
        "<tr><td style='padding:4px;'><strong>Outcome (Trait | GCST ID)</strong></td><td style='padding:4px;'>Unique identifier of the outcome GWAS dataset used in the analysis.</td></tr>",
        "<tr><td style='padding:4px;'><strong>R2 (Exposure)</strong></td><td style='padding:4px;'>Proportion of variance (R²) in the exposure explained by the genetic instruments (SNPs).</td></tr>",
        "<tr><td style='padding:4px;'><strong>R2 (Outcome)</strong></td><td style='padding:4px;'>Proportion of variance (R²) in the outcome explained by the same genetic instruments.</td></tr>",
        "<tr><td style='padding:4px;'><strong>Direction Supported?</strong></td><td style='padding:4px;'>Logical indicator showing whether the inferred causal direction (exposure → outcome) is supported (TRUE if SNPs explain more variance in the exposure than in the outcome).</td></tr>",
        "<tr><td style='padding:4px;'><strong>Steiger P-value</strong></td><td style='padding:4px;'>P‑value from the Steiger directionality test assessing whether the variance explained in the exposure is significantly greater than in the outcome.</td></tr>",
        "<tr><td style='padding:4px;'><strong>Direction Significance</strong></td><td style='padding:4px;'>Indicator of statistically significant support for the proposed causal direction. Reported as significant (***p < 0.001, **p < 0.01, *p < 0.05, and .p < 0.1) and not significant (p ≥ 0.05).</td></tr>",
        "</tbody></table>"
      )),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  observeEvent(input$info_mvmr_diff, {
    
    showModal(modalDialog(
      title = "MVMR Mediation (Difference Method)",
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close"),
      
      HTML("
    <p><strong>Indirect effect calculated as Total − Direct from MVMR.</strong><br>
    </p>

    <table style='width:100%; border-collapse:collapse;'>
      <thead>
        <tr>
          <th style='text-align:left; padding:8px; border-bottom:1px solid #ccc;'>Column</th>
          <th style='text-align:left; padding:8px; border-bottom:1px solid #ccc;'>Description</th>
        </tr>
      </thead>
      <tbody>

        <tr>
          <td style='padding:8px; border-bottom:1px solid #eee;'><strong>Exposure</strong></td>
          <td style='padding:8px; border-bottom:1px solid #eee;'>Unique identifier of the exposure GWAS dataset corresponding to the exposure trait used in the Mendelian Randomization analysis.</td>
        </tr>

        <tr>
          <td style='padding:8px; border-bottom:1px solid #eee;'><strong>Mediator</strong></td>
          <td style='padding:8px; border-bottom:1px solid #eee;'>Unique identifier of the Mediator GWAS dataset corresponding to the Mediator trait analyzed in the Mendelian Randomization analysis.</td>
        </tr>

        <tr>
          <td style='padding:8px; border-bottom:1px solid #eee;'><strong>Outcome</strong></td>
          <td style='padding:8px; border-bottom:1px solid #eee;'>Unique identifier of the outcome GWAS dataset corresponding to the outcome trait analyzed in the Mendelian Randomization analysis.</td>
        </tr>

        <tr>
          <td style='padding:8px; border-bottom:1px solid #eee;'><strong>Pathway Support</strong></td>
          <td style='padding:8px; border-bottom:1px solid #eee;'>Indicates whether both X→M and M→Y|X paths are statistically significant.</td>
        </tr>

        <tr>
          <td style='padding:8px; border-bottom:1px solid #eee;'><strong>Total Effect (X → Y)</strong></td>
          <td style='padding:8px; border-bottom:1px solid #eee;'>Total causal effect of exposure on outcome estimated using univariable MR (Tab 3).</td>
        </tr>

        <tr>
          <td style='padding:8px; border-bottom:1px solid #eee;'><strong>Direct Effect (X → Y | M)</strong></td>
          <td style='padding:8px; border-bottom:1px solid #eee;'>Direct causal effect of exposure on outcome adjusted for mediator using MVMR.</td>
        </tr>

        <tr>
          <td style='padding:8px; border-bottom:1px solid #eee;'><strong>Indirect Effect (Difference)</strong></td>
          <td style='padding:8px; border-bottom:1px solid #eee;'>Indirect (mediated) effect calculated as Total − Direct with 95% confidence interval.</td>
        </tr>

        <tr>
          <td style='padding:8px; border-bottom:1px solid #eee;'><strong>P-value (Difference)</strong></td>
          <td style='padding:8px; border-bottom:1px solid #eee;'>P-value testing whether the indirect (mediated) effect differs from zero.</td>
        </tr>

        <tr>
          <td style='padding:8px;'><strong>Mediated Proportion</strong></td>
          <td style='padding:8px;'>Proportion of the total effect mediated through the mediator with 95% confidence interval.Mediated Proportion = Indirect Effect (Difference) / Total Effect (X → Y)</td>
        </tr>

      </tbody>
    </table>
    ")
      
    ))
    
  })
  observeEvent(input$info_mvmr_prod, {
    
    showModal(modalDialog(
      title = "Two-Sample MR Mediation (Product Method)",
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close"),
      
      HTML("
    <p><strong>Indirect effect calculated as (X→M) × (M→Y|X).</strong><br>
    Confidence intervals derived using the delta method.</p>

    <table style='width:100%; border-collapse:collapse;'>
      <thead>
        <tr>
          <th style='text-align:left; padding:8px; border-bottom:1px solid #ccc;'>Column</th>
          <th style='text-align:left; padding:8px; border-bottom:1px solid #ccc;'>Description</th>
        </tr>
      </thead>
      <tbody>

        <tr>
          <td style='padding:8px; border-bottom:1px solid #eee;'><strong>Exposure</strong></td>
          <td style='padding:8px; border-bottom:1px solid #eee;'>Unique identifier of the exposure GWAS dataset corresponding to the exposure trait used in the Mendelian Randomization analysis.</td>
        </tr>

        <tr>
          <td style='padding:8px; border-bottom:1px solid #eee;'><strong>Mediator</strong></td>
          <td style='padding:8px; border-bottom:1px solid #eee;'>Unique identifier of the Mediator GWAS dataset corresponding to the Mediator trait analyzed in the Mendelian Randomization analysis.</td>
        </tr>

        <tr>
          <td style='padding:8px; border-bottom:1px solid #eee;'><strong>Outcome</strong></td>
          <td style='padding:8px; border-bottom:1px solid #eee;'>Unique identifier of the outcome GWAS dataset corresponding to the outcome trait analyzed in the Mendelian Randomization analysis.</td>
        </tr>

        <tr>
          <td style='padding:8px; border-bottom:1px solid #eee;'><strong>Pathway Support</strong></td>
          <td style='padding:8px; border-bottom:1px solid #eee;'>Indicates whether X→M and M→Y|X are statistically significant.</td>
        </tr>

        <tr>
          <td style='padding:8px; border-bottom:1px solid #eee;'><strong>X → M Effect</strong></td>
          <td style='padding:8px; border-bottom:1px solid #eee;'>Effect of exposure on mediator estimated from two-sample MR (Tab (C) Run Exposure 1 → Exposure 2 / Mediator MR).</td>
        </tr>

        <tr>
          <td style='padding:8px; border-bottom:1px solid #eee;'><strong>M → Y Effect (adjusted for X)</strong></td>
          <td style='padding:8px; border-bottom:1px solid #eee;'>Effect of mediator on outcome adjusted for exposure estimated using MVMR.</td>
        </tr>

        <tr>
          <td style='padding:8px; border-bottom:1px solid #eee;'><strong>Indirect Effect (Product)</strong></td>
          <td style='padding:8px; border-bottom:1px solid #eee;'>Indirect (mediated) effect calculated as (X→M) × (M→Y|X) with 95% confidence interval.</td>
        </tr>

        <tr>
          <td style='padding:8px; border-bottom:1px solid #eee;'><strong>P-value (Product)</strong></td>
          <td style='padding:8px; border-bottom:1px solid #eee;'>P-value for the product-based indirect effect.</td>
        </tr>

        <tr>
          <td style='padding:8px;'><strong>Mediated Proportion</strong></td>
          <td style='padding:8px;'>Proportion of total effect mediated using the product method with 95% confidence interval.Mediated Proportion = Indirect Effect (Product) / Total Effect (X → Y)</td>
        </tr>

      </tbody>
    </table>
    ")
      
    ))
    
  })
  # Similar modals for heterogeneity, pleiotropy, and direction omitted for brevity
  # (copy from original app.R if you need them)
  
  # --------------------------------------------------------------------
  # MVMR PDF download
  # --------------------------------------------------------------------
  output$dl_mvmr_full_xlsx <- downloadHandler(
    filename = function() {
      
      exp1 <- mvmr_rv$used$exp1_gcst
      exp2 <- mvmr_rv$used$exp2_gcst
      out  <- mvmr_rv$used$out_gcst
      
      paste0(
        "MVMR_",
        exp1, "_",
        exp2, "_",
        out, "_",
        Sys.Date(),
        ".xlsx"
      )
    },
    
    content = function(file) {
      
      wb <- openxlsx::createWorkbook()
      
      # Sheet 1 — IVW
      df_ivw <- mvmr_rv$res_ivw
      if (is.null(df_ivw)) df_ivw <- data.frame(Message="No IVW results.")
      
      openxlsx::addWorksheet(wb, "IVW_Estimates")
      openxlsx::writeData(wb, "IVW_Estimates", df_ivw)
      
      # Sheet 2 — Instrument Strength
      df_strength <- mvmr_rv$res_strength
      if (is.null(df_strength)) df_strength <- data.frame(Message="No strength results.")
      
      openxlsx::addWorksheet(wb, "Instrument_Strength")
      openxlsx::writeData(wb, "Instrument_Strength", df_strength)
      
      # Sheet 3 — Pleiotropy
      df_pleio <- mvmr_rv$res_pleio
      if (is.null(df_pleio)) df_pleio <- data.frame(Message="No pleiotropy results.")
      
      openxlsx::addWorksheet(wb, "Pleiotropy_Test")
      openxlsx::writeData(wb, "Pleiotropy_Test", df_pleio)
      
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
}
