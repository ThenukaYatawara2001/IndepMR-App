# R/05_search_module.R
# ============================================================
# IndepMR — Search module (UI + Server)
# Safe independence labeling
# Restriction + selection logic unchanged
# ============================================================

searchUI <- function(id, title = "Search", show_selection_box = FALSE,
                     selection_title = "Selected GCST IDs", extra_ui = NULL) {
  ns <- NS(id)
  tagList(
    div(
      class = "card",
      div(
        class = "card-header",
        div(class = "card-title", title),
        div(class = "card-subtitle",
            "Search by publication information, GCST ID, or disease/trait keyword. First, select a PubMed ID from results, then choose study accession(s).")
      ),
      div(
        class = "card-body",
        if (!is.null(extra_ui)) extra_ui,
        fluidRow(
          column(
            4,
            div(class = "label", "Search by"),
            selectInput(
              ns("search_type"), label = NULL,
              choices = c(
                "PubMed ID",
                "First Author",
                "Journal Name",
                "Study Title",
                "Study Accession (GCST ID)",
                "Disease / Trait Keyword"
              ),
              selected = "Disease / Trait Keyword"
            )
          ),
          column(
            6,
            div(class = "label", "Search term"),
            textInput(
              ns("query"), label = NULL,
              placeholder = "e.g. heart failure, kidney disease, Smith, GCST90012345, 29892016"
            )
          ),
          column(
            2,
            div(class = "label", " "),
            actionButton(ns("go"), "Search", class = "btn-primary w100")
          )
        ),
        
        if (show_selection_box) div(class = "sep"),
        if (show_selection_box) div(
          class = "panel",
          div(class = "panel-title", selection_title),
          verbatimTextOutput(ns("selected_box"))
        ),
        
        div(class = "sep"),
        div(class = "panel",
            div(class = "panel-title", "Results (Select a PUBMEDID)"),
            DTOutput(ns("tbl_main"))
        ),
        
        div(class = "sep"),
        div(class = "panel",
            div(class = "panel-title", "Study Details (Select Study Accession(s))"),
            DTOutput(ns("tbl_details"))
        )
      )
    )
  )
}

searchServer <- function(id, merged, wide,
                         snp_counts = NULL,
                         status_fun_reactive = NULL,
                         filter_only_indep_reactive = NULL,
                         # Optional: compute cohort overlap between selected reference IDs
                         # and candidate GCSTs.  Should return a function taking a vector
                         # of GCST IDs and returning "YES", "NO" or "UNKNOWN".
                         overlap_fun_reactive = NULL,
                         # When TRUE (and the first filter is TRUE), keep only rows
                         # without cohort overlap ("NO").
                         filter_only_non_overlap_reactive = NULL,
                         enable_selection = TRUE,
                         restrict_accessions_reactive = NULL) {
  
  moduleServer(id, function(input, output, session) {
    
    rv <- reactiveValues(
      main = NULL,
      details = data.frame(),
      filtered_rows = merged[0, ],
      selected_gcst = character(0)
    )
    
    apply_independence <- function(df_in) {
      # Annotate each row with independence status and cohort overlap (if requested).
      df2 <- df_in
      # Compute independence status using the provided status function, if any.
      if (!is.null(status_fun_reactive)) {
        sf <- status_fun_reactive()
        if (!is.null(sf) && ("STUDY ACCESSION" %in% names(df2))) {
          df2$`Independence` <- sf(df2$`STUDY ACCESSION`)
          df2$`Ancestry` <- dplyr::recode(
            df2$`Independence`,
            "Independent" = "Different",
            "Not independent" = "Same"
          )
        }
      }
      # Compute cohort overlap status using the provided overlap function, if any.
      if (!is.null(overlap_fun_reactive)) {
        of <- overlap_fun_reactive()
        if (!is.null(of) && ("STUDY ACCESSION" %in% names(df2))) {
          df2$`Independence overlap` <- of(df2$`STUDY ACCESSION`)
        }
      }
      # Apply filtering logic: when filter_only_indep_reactive() is TRUE, keep
      # only rows labelled "Not independent" (i.e. same ancestry).  If
      # filter_only_non_overlap_reactive() is TRUE, further restrict to
      # non-overlapping cohorts (overlap == "NO").
      if (!is.null(filter_only_indep_reactive) && isTRUE(filter_only_indep_reactive())) {
        if ("Independence" %in% names(df2)) {
          df2 <- df2 %>% dplyr::filter(`Independence` == "Not independent")
        }
        if (!is.null(filter_only_non_overlap_reactive) && isTRUE(filter_only_non_overlap_reactive())) {
          if ("Independence overlap" %in% names(df2)) {
            df2 <- df2 %>% dplyr::filter(`Independence overlap` == "NO")
          }
        }
      }
      df2
    }
    
    apply_restriction <- function(df_in) {
      if (is.null(restrict_accessions_reactive)) return(df_in)
      allowed <- restrict_accessions_reactive()
      allowed <- unique(stringr::str_trim(as.character(allowed)))
      allowed <- allowed[allowed != "" & !is.na(allowed)]
      if (length(allowed) == 0) return(df_in)
      df_in %>% dplyr::filter(`STUDY ACCESSION` %in% allowed)
    }
    
    observeEvent(input$go, {
      
      notif_id <- showNotification(
        "🔍 Searching… please wait",
        type = "message",
        duration = NULL
      )
      
      on.exit(removeNotification(notif_id), add = TRUE)
      
      type <- input$search_type
      q <- stringr::str_trim(input$query)
      
      if (is.na(q) || q == "") {
        rv$main <- data.frame(Message = "Please enter a search term.")
        rv$details <- data.frame()
        rv$filtered_rows <- merged[0, ]
        rv$selected_gcst <- character(0)
        return()
      }
      
      if (type %in% c("PubMed ID","First Author","Journal Name","Study Title","Study Accession (GCST ID)")) {
        
        df_hit <- merged %>% apply_restriction()
        
        if (type == "PubMed ID") {
          df_hit <- df_hit %>% dplyr::filter(stringr::str_detect(PUBMEDID, stringr::fixed(q)))
        } else if (type == "First Author") {
          df_hit <- df_hit %>% dplyr::filter(stringr::str_detect(stringr::str_to_lower(`FIRST AUTHOR`), stringr::fixed(stringr::str_to_lower(q))))
        } else if (type == "Journal Name") {
          df_hit <- df_hit %>% dplyr::filter(stringr::str_detect(stringr::str_to_lower(JOURNAL), stringr::fixed(stringr::str_to_lower(q))))
        } else if (type == "Study Title") {
          df_hit <- df_hit %>% dplyr::filter(stringr::str_detect(stringr::str_to_lower(STUDY), stringr::fixed(stringr::str_to_lower(q))))
        } else if (type == "Study Accession (GCST ID)") {
          df_hit <- df_hit %>% dplyr::filter(stringr::str_detect(stringr::str_to_lower(`STUDY ACCESSION`), stringr::fixed(stringr::str_to_lower(q))))
        }
        
        df_hit <- apply_independence(df_hit)
        
        if (nrow(df_hit) == 0) {
          rv$main <- data.frame(Message = "No results found.")
          rv$details <- data.frame()
          rv$filtered_rows <- merged[0, ]
          rv$selected_gcst <- character(0)
          return()
        }
        
        rv$filtered_rows <- df_hit
        rv$main <- publication_summary(df_hit)
        rv$details <- data.frame()
        rv$selected_gcst <- character(0)
        
      } else {
        
        accessions <- wide_keyword_to_accessions(wide, q)
        
        if (!is.null(restrict_accessions_reactive)) {
          allowed <- restrict_accessions_reactive()
          allowed <- unique(stringr::str_trim(as.character(allowed)))
          allowed <- allowed[allowed != "" & !is.na(allowed)]
          if (length(allowed) > 0) accessions <- intersect(accessions, allowed)
        }
        
        if (length(accessions) == 0) {
          rv$main <- data.frame(Message = "No studies matched this keyword in the wide ontology table.")
          rv$details <- data.frame()
          rv$filtered_rows <- merged[0, ]
          rv$selected_gcst <- character(0)
          return()
        }
        
        df_hit <- fetch_by_accessions(merged, accessions) %>% apply_independence()
        
        if (nrow(df_hit) == 0) {
          rv$main <- data.frame(Message = "Keyword matched studies, but no results found in merged table.")
          rv$details <- data.frame()
          rv$filtered_rows <- merged[0, ]
          rv$selected_gcst <- character(0)
          return()
        }
        
        rv$filtered_rows <- df_hit
        rv$main <- publication_summary(df_hit)
        rv$details <- data.frame()
        rv$selected_gcst <- character(0)
      }
      
    })
    
    output$tbl_main <- renderDT({
      req(rv$main)
      df <- rv$main
      if ("Message" %in% names(df)) {
        return(DT::datatable(df, rownames = FALSE, options = list(dom = "t")))
      }
      DT::datatable(df, selection = "single", rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE, dom = "tip"))
    })
    
    observeEvent(input$tbl_main_rows_selected, {
      req(rv$filtered_rows)
      sel <- input$tbl_main_rows_selected
      if (is.null(sel) || length(sel) == 0) return()
      if (is.null(rv$main) || "Message" %in% names(rv$main)) return()
      
      pubmed_id <- rv$main$PUBMEDID[sel]
      if (is.null(pubmed_id) || is.na(pubmed_id)) return()
      
      rv$details <- publication_details(rv$filtered_rows, pubmed_id, snp_counts = snp_counts)
      if ("Study Accession" %in% names(rv$details)) {
        
        rv$details <- rv$details %>%
          dplyr::left_join(
            merged %>%
              dplyr::select(`STUDY ACCESSION`, `Data Type`),
            by = c("Study Accession" = "STUDY ACCESSION")
          )
        
      }
      # Annotate study details with independence status and cohort overlap, then
      # apply any requested filtering.  We compute independence for each
      # accession using the provided status function (if available), compute
      # cohort overlap via the overlap function (if provided), and then filter
      # according to the user's checkboxes.  The independence filter is
      # inverted: when the first checkbox is active, we retain only "Not
      # independent" (same‑ancestry) entries.  When the second checkbox is
      # active in addition, we further restrict to studies whose cohorts do
      # not overlap with the reference cohort set (overlap == "NO").
      #
      # Compute independence status on the friendly "Study Accession" column.
      if (!is.null(status_fun_reactive)) {
        sf <- status_fun_reactive()
        if (!is.null(sf) && ("Study Accession" %in% names(rv$details))) {
          rv$details$`Independence` <- sf(rv$details$`Study Accession`)
          rv$details$`Ancestry` <- dplyr::recode(
            rv$details$`Independence`,
            "Independent" = "Different",
            "Not independent" = "Same"
          )
        }
      }
      # Compute cohort overlap if an overlap function is supplied.
      if (!is.null(overlap_fun_reactive)) {
        of <- overlap_fun_reactive()
        if (!is.null(of) && ("Study Accession" %in% names(rv$details))) {
          rv$details$`Independence overlap` <- of(rv$details$`Study Accession`)
        }
      }
      # Apply filtering: when the first checkbox is TRUE, retain only
      # "Not independent" rows.  When the second checkbox is also TRUE,
      # further retain only those rows with no cohort overlap ("NO").
      if (!is.null(filter_only_indep_reactive) && isTRUE(filter_only_indep_reactive())) {
        if ("Independence" %in% names(rv$details)) {
          rv$details <- rv$details %>% dplyr::filter(`Independence` == "Not independent")
        }
        if (!is.null(filter_only_non_overlap_reactive) && isTRUE(filter_only_non_overlap_reactive())) {
          if ("Independence overlap" %in% names(rv$details)) {
            rv$details <- rv$details %>% dplyr::filter(`Independence overlap` == "NO")
          }
        }
      }
      
      rv$selected_gcst <- character(0)
    })
    
    output$tbl_details <- renderDT({
      df <- rv$details
      
      # Remove internal Independence column
      if ("Independence" %in% names(df)) {
        df <- df %>% dplyr::select(-Independence)
      }
      
      # Rename Independence overlap → Cohort Overlap (DISPLAY ONLY)
      if ("Independence overlap" %in% names(df)) {
        colnames(df)[colnames(df) == "Independence overlap"] <- "Cohort Overlap"
      }
      
      if (nrow(df) == 0) {
        return(DT::datatable(
          data.frame(Message = "Select a publication above to view its Study Accession(s)."),
          rownames = FALSE,
          options = list(dom = "t")
        ))
      }
      
      if ("Data Type" %in% names(df)) {
        df <- df %>%
          dplyr::relocate(`Data Type`, .after = `Study Accession`)
      }
      
      DT::datatable(df,
                    selection = if (enable_selection) "multiple" else "single",
                    rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE, dom = "tip"))
    })
    
    observeEvent(input$tbl_details_rows_selected, {
      if (!enable_selection) return()
      df <- rv$details
      sel <- input$tbl_details_rows_selected
      if (is.null(sel) || length(sel) == 0 || nrow(df) == 0) {
        rv$selected_gcst <- character(0)
        return()
      }
      
      if (!("Study Accession" %in% names(df))) {
        rv$selected_gcst <- character(0)
        return()
      }
      
      ids <- df$`Study Accession`[sel]
      ids <- unique(stringr::str_trim(as.character(ids)))
      ids <- ids[ids != "" & !is.na(ids)]
      rv$selected_gcst <- ids
    })
    
    output$selected_box <- renderText({
      if (!enable_selection) return("")
      ids <- rv$selected_gcst
      if (length(ids) == 0) return("None selected yet.")
      paste(ids, collapse = ", ")
    })
    
    list(selected_gcst = reactive(rv$selected_gcst))
  })
}