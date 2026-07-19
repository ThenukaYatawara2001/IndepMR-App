# R/05_search_module.R
# ============================================================
# IndepMR — Search Module (UI + Server)
#
# Displays GWAS study search results with ancestry-based
# independence classification derived from EXTRACTED_TERMS.
#
# COLUMN SEMANTICS in the study details table
# ─────────────────────────────────────────────────────────────
#   Ancestry        "Same" / "Different" / "Unknown"
#                   Whether the candidate shares a main ancestry
#                   group with the selected exposure(s).
#
#   Sub-group       Fine-grained within-ancestry classification
#                   (only shown when Ancestry == "Same"):
#
#     "Independent"   — same main group, different known
#                       sub-populations (e.g. Scandinavian vs
#                       Iberian within European).
#     "Overlapping"   — same sub-population cluster.
#     "Broad only"    — only continental label available; finer
#                       resolution not possible.
#     "—"             — not applicable.
#
# ============================================================

library(shiny)
library(DT)
library(dplyr)
library(stringr)

# ─────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────

searchUI <- function(id,
                     title           = "Search",
                     show_selection_box = FALSE,
                     selection_title    = "Selected GCST IDs",
                     extra_ui           = NULL) {
  ns <- NS(id)
  tagList(
    div(
      class = "card",
      div(
        class = "card-header",
        div(class = "card-title", title),
        div(class = "card-subtitle",
            "Search by publication information, GCST ID, or disease/trait keyword. ",
            "First select a PubMed ID from results, then choose study accession(s).")
      ),
      div(
        class = "card-body",
        if (!is.null(extra_ui)) extra_ui,
        fluidRow(
          column(4,
            div(class = "label", "Search by"),
            selectInput(ns("search_type"), label = NULL,
                        choices  = c("PubMed ID", "First Author", "Journal Name",
                                     "Study Title", "Study Accession (GCST ID)",
                                     "Disease / Trait Keyword"),
                        selected = "Disease / Trait Keyword")),
          column(6,
            div(class = "label", "Search term"),
            textInput(ns("query"), label = NULL,
                      placeholder = "e.g. heart failure, kidney disease, Smith, GCST90012345")),
          column(2,
            div(class = "label", "\u00a0"),
            actionButton(ns("go"), "Search", class = "btn-primary w100"))
        ),
        if (show_selection_box) div(class = "sep"),
        if (show_selection_box) div(
          class = "panel",
          div(class = "panel-title", selection_title),
          verbatimTextOutput(ns("selected_box"))
        ),
        div(class = "sep"),
        div(class = "panel",
            div(class = "panel-title", "Results (Select a PubMed ID)"),
            DT::DTOutput(ns("tbl_main"))),
        div(class = "sep"),
        div(class = "panel",
            div(class = "panel-title", "Study Details (Select Study Accession(s))"),
            DT::DTOutput(ns("tbl_details")),
            div(class = "sep"),
            div(class = "label", "Message"),
            uiOutput(ns("message_box")))
      )
    )
  )
}

# ─────────────────────────────────────────────────────────────
# Server
# ─────────────────────────────────────────────────────────────

searchServer <- function(id,
                         merged,
                         wide,
                         snp_counts                       = NULL,
                         status_fun_reactive              = NULL,
                         filter_only_indep_reactive       = NULL,
                         subgroup_indep_fun_reactive      = NULL,
                         filter_only_non_overlap_reactive = NULL,
                         enable_selection                 = TRUE,
                         restrict_accessions_reactive     = NULL,
                         # legacy alias — ignored, kept for compatibility
                         overlap_fun_reactive             = NULL) {

  moduleServer(id, function(input, output, session) {

    rv <- reactiveValues(
      main          = NULL,
      details       = data.frame(),
      filtered_rows = merged[0, ],
      selected_gcst = character(0)
    )

    # ── Label helpers ────────────────────────────────────────

    .ancestry_label <- function(independence_codes) {
      dplyr::case_when(
        independence_codes == "Independent"    ~ "Different",
        independence_codes == "Not independent"~ "Same",
        TRUE                                   ~ "Unknown"
      )
    }

    .subgroup_label <- function(codes) {
      dplyr::case_when(
        codes == "SUBGROUP_INDEPENDENT"  ~ "Independent",
        codes == "SUBGROUP_OVERLAPPING"  ~ "Overlapping",
        codes == "UMBRELLA_ONLY"         ~ "Broad only",
        codes == "OTHER_PRESENT"         ~ "Broad / unknown",
        TRUE                             ~ "\u2014"
      )
    }

    # ── Apply independence annotation + optional filtering ───

    .annotate_and_filter <- function(df) {
      if (!"STUDY ACCESSION" %in% names(df)) return(df)

      # Level 1: ancestry status
      if (!is.null(status_fun_reactive)) {
        sf <- status_fun_reactive()
        if (!is.null(sf)) {
          df$`_indep_raw` <- sf(df$`STUDY ACCESSION`)
          df$Ancestry     <- .ancestry_label(df$`_indep_raw`)
        }
      }

      # Level 2: within-ancestry sub-group status
      if (!is.null(subgroup_indep_fun_reactive)) {
        sgf <- subgroup_indep_fun_reactive()
        if (!is.null(sgf)) {
          df$`_subgroup_raw` <- sgf(df$`STUDY ACCESSION`)
          df$`Sub-group`     <- .subgroup_label(df$`_subgroup_raw`)
        }
      }

      # Filter: first checkbox (same ancestry only)
      if (!is.null(filter_only_indep_reactive) &&
          isTRUE(filter_only_indep_reactive())) {
        if ("_indep_raw" %in% names(df))
          df <- df %>% dplyr::filter(`_indep_raw` == "Not independent")

        # Filter: second checkbox (sub-group independent only)
        if (!is.null(filter_only_non_overlap_reactive) &&
            isTRUE(filter_only_non_overlap_reactive())) {
          if ("_subgroup_raw" %in% names(df))
            df <- df %>%
              dplyr::filter(`_subgroup_raw` == "SUBGROUP_INDEPENDENT")
        }
      }

      df
    }

    # ── Restrict to allowed accessions ─────────────────────

    .apply_restriction <- function(df) {
      if (is.null(restrict_accessions_reactive)) return(df)
      allowed <- unique(stringr::str_trim(
        as.character(restrict_accessions_reactive())))
      allowed <- allowed[allowed != "" & !is.na(allowed)]
      if (length(allowed) == 0) return(df)
      df %>% dplyr::filter(`STUDY ACCESSION` %in% allowed)
    }

    # ── Search ──────────────────────────────────────────────

    observeEvent(input$go, {
      nid <- showNotification("\U0001f50d Searching\u2026 please wait",
                              type = "message", duration = NULL)
      on.exit(removeNotification(nid), add = TRUE)

      type <- input$search_type
      q    <- stringr::str_trim(input$query)

      if (is.na(q) || q == "") {
        rv$main          <- data.frame(Message = "Please enter a search term.")
        rv$details       <- data.frame()
        rv$filtered_rows <- merged[0, ]
        rv$selected_gcst <- character(0)
        return()
      }

      df_hit <- .apply_restriction(merged)

      if (type != "Disease / Trait Keyword") {
        df_hit <- switch(type,
          "PubMed ID" =
            df_hit %>% dplyr::filter(
              stringr::str_detect(PUBMEDID, stringr::fixed(q))),
          "First Author" =
            df_hit %>% dplyr::filter(
              stringr::str_detect(
                stringr::str_to_lower(`FIRST AUTHOR`), stringr::fixed(tolower(q)))),
          "Journal Name" =
            df_hit %>% dplyr::filter(
              stringr::str_detect(
                stringr::str_to_lower(JOURNAL), stringr::fixed(tolower(q)))),
          "Study Title" =
            df_hit %>% dplyr::filter(
              stringr::str_detect(
                stringr::str_to_lower(STUDY), stringr::fixed(tolower(q)))),
          "Study Accession (GCST ID)" =
            df_hit %>% dplyr::filter(
              stringr::str_detect(
                stringr::str_to_lower(`STUDY ACCESSION`), stringr::fixed(tolower(q)))),
          df_hit
        )

        df_hit <- .annotate_and_filter(df_hit)

        if (nrow(df_hit) == 0) {
          rv$main          <- data.frame(Message = "No results found.")
          rv$details       <- data.frame()
          rv$filtered_rows <- merged[0, ]
          rv$selected_gcst <- character(0)
          return()
        }

        rv$filtered_rows <- df_hit
        rv$main          <- publication_summary(df_hit)
        rv$details       <- data.frame()
        rv$selected_gcst <- character(0)

      } else {
        accessions <- wide_keyword_to_accessions(wide, q)
        if (!is.null(restrict_accessions_reactive)) {
          allowed <- unique(stringr::str_trim(
            as.character(restrict_accessions_reactive())))
          allowed <- allowed[allowed != "" & !is.na(allowed)]
          if (length(allowed) > 0)
            accessions <- intersect(accessions, allowed)
        }

        if (length(accessions) == 0) {
          rv$main          <- data.frame(Message = "No studies matched this keyword.")
          rv$details       <- data.frame()
          rv$filtered_rows <- merged[0, ]
          rv$selected_gcst <- character(0)
          return()
        }

        df_hit <- fetch_by_accessions(merged, accessions) %>%
          .annotate_and_filter()

        if (nrow(df_hit) == 0) {
          rv$main          <- data.frame(Message = "No results found after filters.")
          rv$details       <- data.frame()
          rv$filtered_rows <- merged[0, ]
          rv$selected_gcst <- character(0)
          return()
        }

        rv$filtered_rows <- df_hit
        rv$main          <- publication_summary(df_hit)
        rv$details       <- data.frame()
        rv$selected_gcst <- character(0)
      }
    })

    # ── Render publication summary ──────────────────────────

    output$tbl_main <- DT::renderDT({
      req(rv$main)
      df <- rv$main
      if ("Message" %in% names(df))
        return(DT::datatable(df, rownames = FALSE,
                             options = list(dom = "t")))
      DT::datatable(df, selection = "single", rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE, dom = "tip"))
    })

    # ── Publication selected → populate details ─────────────

    observeEvent(input$tbl_main_rows_selected, {
      req(rv$filtered_rows)
      sel <- input$tbl_main_rows_selected
      if (is.null(sel) || length(sel) == 0) return()
      if (is.null(rv$main) || "Message" %in% names(rv$main)) return()

      pubmed_id <- rv$main$PUBMEDID[sel]
      if (is.null(pubmed_id) || is.na(pubmed_id)) return()

      det <- publication_details(rv$filtered_rows, pubmed_id,
                                 snp_counts = snp_counts)

      # Join Data Type if available
      if ("Study Accession" %in% names(det) &&
          "Data Type" %in% names(merged)) {
        det <- det %>%
          dplyr::left_join(
            merged %>% dplyr::select(`STUDY ACCESSION`, `Data Type`),
            by = c("Study Accession" = "STUDY ACCESSION"))
      }

      # Ancestry level-1
      if (!is.null(status_fun_reactive)) {
        sf <- status_fun_reactive()
        if (!is.null(sf) && "Study Accession" %in% names(det)) {
          det$`_indep_raw` <- sf(det$`Study Accession`)
          det$Ancestry     <- .ancestry_label(det$`_indep_raw`)
        }
      }

      # Within-ancestry level-2
      if (!is.null(subgroup_indep_fun_reactive)) {
        sgf <- subgroup_indep_fun_reactive()
        if (!is.null(sgf) && "Study Accession" %in% names(det)) {
          det$`_subgroup_raw` <- sgf(det$`Study Accession`)
          det$`Sub-group`     <- .subgroup_label(det$`_subgroup_raw`)
        }
      }

      # Apply filters
      if (!is.null(filter_only_indep_reactive) &&
          isTRUE(filter_only_indep_reactive())) {
        if ("_indep_raw" %in% names(det))
          det <- det %>% dplyr::filter(`_indep_raw` == "Not independent")
        if (!is.null(filter_only_non_overlap_reactive) &&
            isTRUE(filter_only_non_overlap_reactive()) &&
            "_subgroup_raw" %in% names(det))
          det <- det %>%
            dplyr::filter(`_subgroup_raw` == "SUBGROUP_INDEPENDENT")
      }

      rv$details       <- det
      rv$selected_gcst <- character(0)
    })

    # ── Render details table ────────────────────────────────

    output$tbl_details <- DT::renderDT({
      df <- rv$details

      # Drop internal/helper columns
      drop <- c("_indep_raw", "_subgroup_raw", "COHORT")
      df   <- df[, !names(df) %in% drop, drop = FALSE]

      if (nrow(df) == 0)
        return(DT::datatable(
          data.frame(Message = "Select a publication above to view its study accessions."),
          rownames = FALSE, options = list(dom = "t")))

      # Reorder: Data Type right after Study Accession
      if ("Data Type" %in% names(df) && "Study Accession" %in% names(df))
        df <- df %>% dplyr::relocate(`Data Type`, .after = `Study Accession`)

      # Reorder: Ancestry + Sub-group near front
      ancestry_cols <- intersect(c("Ancestry", "Sub-group"), names(df))
      if (length(ancestry_cols) > 0 && "Study Accession" %in% names(df)) {
        rest <- setdiff(names(df), c("Study Accession", ancestry_cols))
        df   <- df[, c("Study Accession", ancestry_cols, rest), drop = FALSE]
      }

      dt <- DT::datatable(
        df,
        selection = if (enable_selection) "multiple" else "single",
        rownames  = FALSE,
        options   = list(pageLength = 10, scrollX = TRUE, dom = "tip")
      )

      # Colour Ancestry column
      if ("Ancestry" %in% names(df)) {
        dt <- dt %>% DT::formatStyle(
          "Ancestry",
          fontWeight = "bold",
          color = DT::styleEqual(
            c("Same", "Different", "Unknown"),
            c("#34d399", "#fb7185", "#fbbf24")
          )
        )
      }

      # Colour Sub-group column
      if ("Sub-group" %in% names(df)) {
        dt <- dt %>% DT::formatStyle(
          "Sub-group",
          color = DT::styleEqual(
            c("Independent", "Overlapping", "Broad only",
              "Broad / unknown", "\u2014"),
            c("#34d399", "#fb7185", "#fbbf24", "#f97316", "#94a3b8")
          )
        )
      }

      dt
    })

    # ── Accessions selected ─────────────────────────────────

    observeEvent(input$tbl_details_rows_selected, {
      if (!enable_selection) return()
      df  <- rv$details
      sel <- input$tbl_details_rows_selected
      if (is.null(sel) || length(sel) == 0 || nrow(df) == 0) {
        rv$selected_gcst <- character(0)
        return()
      }
      col <- if ("Study Accession" %in% names(df)) "Study Accession" else
             if ("STUDY ACCESSION" %in% names(df)) "STUDY ACCESSION" else NULL
      if (is.null(col)) { rv$selected_gcst <- character(0); return() }
      ids  <- unique(stringr::str_trim(as.character(df[[col]][sel])))
      rv$selected_gcst <- ids[ids != "" & !is.na(ids)]
    })

    # ── Selection display ───────────────────────────────────

    output$selected_box <- renderText({
      if (!enable_selection) return("")
      ids <- rv$selected_gcst
      if (length(ids) == 0) return("None selected yet.")
      paste(ids, collapse = ", ")
    })

    # ── Contextual message box (HTML) ───────────────────────

    output$message_box <- renderUI({
      same_on <- !is.null(filter_only_indep_reactive) &&
                 isTRUE(filter_only_indep_reactive())
      sub_on  <- same_on &&
                 !is.null(filter_only_non_overlap_reactive) &&
                 isTRUE(filter_only_non_overlap_reactive())

      df <- rv$details
      if (nrow(df) == 0 || "Message" %in% names(df)) return(NULL)
      if (!same_on) return(NULL)

      if ("_subgroup_raw" %in% names(df)) {
        has_umbrella  <- any(df$`_subgroup_raw` == "UMBRELLA_ONLY", na.rm = TRUE)
        has_sub_indep <- any(df$`_subgroup_raw` == "SUBGROUP_INDEPENDENT", na.rm = TRUE)

        if (!sub_on && has_umbrella && !has_sub_indep) {
          return(div(class = "ancestry-notice ancestry-info",
            HTML(paste0(
              "\u2139\ufe0f  These studies share the same main ancestry group as the exposure, ",
              "but only a broad continental label is available (e.g. <em>European</em>). ",
              "Sub-population independence cannot be verified at this level. ",
              "The sub-group filter is therefore hidden."
            ))))
        }

        if (sub_on && has_sub_indep) {
          return(div(class = "ancestry-notice ancestry-info",
            HTML(paste0(
              "\u2139\ufe0f  These studies come from distinct sub-populations within the ",
              "same main ancestry group (e.g. Scandinavian vs. Iberian within European). ",
              "Sample overlap is unlikely; however, LD patterns, haplotype blocks, and ",
              "allele frequencies may differ between sub-populations. ",
              "Verify that your genetic instruments replicate well across these groups."
            ))))
        }
      }

      NULL
    })

    # ── Return selected GCSTs ───────────────────────────────

    list(selected_gcst = reactive(rv$selected_gcst))
  })
}
