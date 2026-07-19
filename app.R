# app.R
# ============================================================
# IndepMR — Application entry point
#
# Execution order:
#   1. library(shiny) + library(DT)
#   2. source(global.R)  → loads all R/ modules and DATA
#   3. Define UI helper functions + build ui object
#   4. source(R/server_original.R)  → defines server function
#   5. shinyApp(ui, server)
# ============================================================

library(shiny)
library(DT)

# Project root — single definition
APP_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Load all modules, packages, and data
source(file.path(APP_DIR, "global.R"), local = FALSE)

# Register static assets (logo etc.)
addResourcePath("assets", APP_DIR)

# ============================================================
# Custom CSS
# ============================================================

custom_css <- "
:root {
  --bg1:    #0b1220;
  --card:   #0f1b33;
  --panel:  #0b1730;
  --text:   #eaf0ff;
  --muted:  #a9b4c7;
  --line:   rgba(255,255,255,0.10);
  --accent: #4da3ff;
  --accent2:#6d72ff;
  --ok:     #34d399;
  --bad:    #fb7185;
  --unk:    #fbbf24;
}
body {
  background: radial-gradient(1200px 750px at 18% 12%,#163a7a 0%,var(--bg1) 45%,#070b14 100%);
  color: var(--text);
  font-family: -apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,Arial,sans-serif;
}
.container-fluid { max-width:1350px; }
.app-header {
  margin-top:12px; border-radius:18px; padding:18px;
  background:linear-gradient(135deg,rgba(77,163,255,0.22),rgba(109,114,255,0.14));
  border:1px solid rgba(255,255,255,0.12);
  box-shadow:0 14px 40px rgba(0,0,0,0.35);
}
.app-title  { font-size:26px; font-weight:900; margin:0; }
.app-subtitle { margin-top:6px; color:var(--muted); font-size:13px; line-height:1.35; }
.card {
  margin-top:14px; border-radius:18px;
  background:rgba(15,27,51,0.86);
  border:1px solid rgba(255,255,255,0.12);
  box-shadow:0 12px 36px rgba(0,0,0,0.32); overflow:hidden;
}
.card-header {
  padding:14px 16px;
  background:linear-gradient(180deg,rgba(255,255,255,0.06),rgba(255,255,255,0.02));
  border-bottom:1px solid var(--line);
}
.card-title    { font-size:18px; font-weight:900; margin:0; }
.card-subtitle { margin-top:6px; color:var(--muted); font-size:12px; }
.card-body     { padding:14px 16px 16px 16px; }
.label         { color:var(--muted); font-size:12px; margin-bottom:6px; }
.sep           { height:14px; }
.panel {
  border-radius:16px; background:rgba(11,23,48,0.62);
  border:1px solid rgba(255,255,255,0.10); padding:12px;
}
.panel-title { font-size:13px; font-weight:900; margin-bottom:10px; }
pre {
  color:var(--text); background:rgba(0,0,0,0.18);
  border:1px solid rgba(255,255,255,0.10); border-radius:12px;
  padding:10px; white-space:pre-wrap;
}
.form-control,.selectize-input {
  background:rgba(8,16,34,0.70)!important; color:var(--text)!important;
  border:1px solid rgba(255,255,255,0.14)!important; border-radius:12px!important;
}
.selectize-dropdown,.selectize-dropdown-content {
  background:#081022!important; color:var(--text)!important;
  border:1px solid rgba(255,255,255,0.14)!important;
}
.btn-primary {
  border:0!important; border-radius:12px!important; height:40px!important;
  background:linear-gradient(90deg,rgba(77,163,255,0.95),rgba(109,114,255,0.95))!important;
  color:white!important; font-weight:800!important;
}
.btn-ghost {
  border-radius:12px!important; height:40px!important;
  background:rgba(255,255,255,0.06)!important; color:var(--text)!important;
  font-weight:800!important; border:1px solid rgba(255,255,255,0.12)!important;
}
.w100 { width:100%; }
.nav-tabs { border-bottom:1px solid rgba(255,255,255,0.10); margin-top:12px; }
.nav-tabs > li > a {
  border-radius:14px 14px 0 0!important;
  background:rgba(255,255,255,0.05)!important;
  border:1px solid rgba(255,255,255,0.10)!important;
  color:var(--muted)!important; font-weight:800!important; margin-right:8px!important;
}
.nav-tabs > li.active > a {
  color:var(--text)!important;
  background:linear-gradient(180deg,rgba(77,163,255,0.20),rgba(255,255,255,0.05))!important;
  border-color:rgba(77,163,255,0.35)!important;
}
table.dataTable { border-collapse:separate!important; border-spacing:0 8px!important; }
table.dataTable thead th {
  background:rgba(255,255,255,0.06)!important; color:var(--text)!important; border:none!important;
}
table.dataTable tbody tr {
  background:rgba(255,255,255,0.04)!important;
  border:1px solid rgba(255,255,255,0.10)!important; border-radius:12px!important;
}
table.dataTable tbody td { border:none!important; color:var(--text)!important; }
.dataTables_info,.dataTables_paginate,.dataTables_length,.dataTables_filter {
  color:var(--muted)!important;
}
.note { color:var(--muted); font-size:12px; line-height:1.35; }
.spinWrap { display:flex; align-items:center; gap:12px; }
.spinner {
  width:22px; height:22px;
  border:3px solid rgba(255,255,255,0.25);
  border-top:3px solid rgba(77,163,255,0.95);
  border-radius:50%; animation:spin 0.9s linear infinite;
}
@keyframes spin { to { transform:rotate(360deg); } }
.modal-content {
  background:rgba(15,27,51,0.96)!important; color:var(--text)!important;
  border:1px solid rgba(255,255,255,0.14)!important; border-radius:16px!important;
}
.header-card {
  display:flex; align-items:center;
  border:1px solid rgba(255,255,255,0.35); border-radius:16px;
  overflow:hidden; margin-top:16px;
}
.header-card .logo-col {
  background:#ffffff; padding:20px; display:flex;
  align-items:center; justify-content:center; min-width:180px;
}
.header-card .text-col {
  background:rgba(11,23,48,0.62); color:var(--text); padding:20px;
  flex:1; display:flex; flex-direction:column; justify-content:center;
}
.header-title   { font-size:22px; font-weight:900; margin:0; line-height:1.2; }
.header-subtitle { margin-top:8px; font-size:12px; line-height:1.35; color:var(--muted); }
.tab7-box    { border:2px solid #1f1f1f; border-radius:8px; padding:14px; background:#fff; min-height:240px; }
.tab7-runbox { display:flex; justify-content:center; align-items:center; min-height:240px; }
.btn-runmodel { font-size:22px; font-weight:700; padding:18px 38px; border-radius:14px; }
/* ── Ancestry context message boxes ────────────────────────── */
.ancestry-notice {
  margin-top:10px; padding:10px 14px;
  border-radius:8px; font-size:0.88rem; line-height:1.55;
}
.ancestry-notice.ancestry-warn {
  background:rgba(251,113,133,0.12); border-left:3px solid #fb7185; color:#fca5a5;
}
.ancestry-notice.ancestry-ok {
  background:rgba(52,211,153,0.10); border-left:3px solid #34d399; color:#6ee7b7;
}
.ancestry-notice.ancestry-info {
  background:rgba(251,191,36,0.10); border-left:3px solid #fbbf24; color:#fcd34d;
}
"

# ============================================================
# JavaScript helpers
# ============================================================

custom_scroll_script <- function() {
  tags$script(HTML(
    "Shiny.addCustomMessageHandler('scrollTo',function(msg){
       var el=document.getElementById(msg.id);
       if(el) el.scrollIntoView({behavior:'smooth',block:'start'});
     });"
  ))
}

# ============================================================
# UI builder functions
# ============================================================

create_header_card <- function() {
  div(
    class = "app-header",
    div(
      class = "header-card",
      div(
        class = "logo-col",
        tags$img(src = "assets/logo.png", style = "height:120px;width:auto;")
      ),
      div(
        class = "text-col",
        h1(class = "header-title",
           "An Ontology-Driven, Ancestry-Aware Automated Framework for Two-Sample and Multivariable Mendelian Randomization with Mediation Analysis Using Public GWAS Summary Data"),
        div(class = "header-subtitle",
            "IndepMR integrates GWAS Catalog and OpenGWAS resources to automatically identify ancestry-independent studies for valid causal inference. Supports 2SMR, MVMR, and mediation analysis.")
      )
    ),
    div(class = "sep"),
    div(
      class = "panel",
      style = "max-width:720px;margin:0 auto;text-align:center;padding:18px;",
      div(class = "panel-title", style = "font-size:18px;", "OpenGWAS Authentication"),
      div(class = "note", style = "font-size:13px;margin-bottom:14px;",
          "Paste a valid OpenGWAS JWT token before running any 2SMR or MVMR analysis."),
      passwordInput("opengwas_jwt", "OpenGWAS JWT Token",
                    placeholder = "Paste your OpenGWAS JWT token here", width = "100%"),
      tags$a(href = "https://api.opengwas.io/profile/", target = "_blank",
             class = "btn btn-ghost", "Get OpenGWAS Token")
    )
  )
}

# ── Tab 1: Select Exposure 1 ──────────────────────────────────
create_select_exposure_tab <- function() {
  tabPanel(
    "Select Exposure 1",
    div(
      searchUI("exposure", "Select Exposure 1",
               show_selection_box = TRUE,
               selection_title    = "Selected Exposure 1 GCST IDs"),
      div(class = "sep"),
      tags$a(href = "https://www.ebi.ac.uk/gwas/", target = "_blank",
             class = "btn btn-ghost", "More study details (GWAS Catalog)")
    )
  )
}

# ── Tab 2: Select Outcome ────────────────────────────────────
create_select_outcome_tab <- function() {
  tabPanel(
    "Select Outcome",
    div(
      searchUI(
        "outcome", "Select Outcome",
        show_selection_box = TRUE,
        selection_title    = "Selected Outcome GCST IDs",
        extra_ui = div(
          class = "panel", style = "margin-bottom:8px;",
          div(class = "panel-title", "Outcome independence options (vs Exposure 1)"),

          # Always shown: filter to same-ancestry outcomes
          checkboxInput("show_only_independent",
                        "Show only same-ancestry Outcome studies (recommended)",
                        value = FALSE),

          # Conditionally rendered by server — only when:
          #   (a) first checkbox is ON, AND
          #   (b) exposure has specific subgroup/member-level ancestry resolution
          uiOutput("outcome_subgroup_checkbox_ui"),

          # Contextual guidance message (server-rendered)
          uiOutput("outcome_ancestry_msg")
        )
      ),
      div(class = "sep"),
      tags$a(href = "https://www.ebi.ac.uk/gwas/", target = "_blank",
             class = "btn btn-ghost", "More study details (GWAS Catalog)")
    )
  )
}

# ── Tab 3: Two-sample MR ─────────────────────────────────────
create_two_sample_mr_tab <- function() {
  tabPanel(
    "Two-sample MR",
    div(
      class = "card",
      div(class = "card-header",
          div(class = "card-title", "Two-sample MR"),
          div(class = "card-subtitle",
              "Adjust instrument selection parameters, choose MR methods, and evaluate causal effects with sensitivity diagnostics.")),
      div(
        class = "card-body",
        fluidRow(
          # ── Left column: settings + downloads ──────────────
          column(
            4,
            div(class = "panel",
                div(class = "panel-title", "Instrument settings"),
                checkboxInput("do_clump", "Use LD clumping (OpenGWAS)", value = TRUE),
                numericInput("p_thresh",  "Instrument p-value (p1)", value = 5e-8, min = 0, max = 1, step = 1e-8),
                sliderInput("clump_r2", "Clump r\u00b2", min = 0.0001, max = 0.2, value = 0.01, step = 0.0005),
                numericInput("clump_kb", "Clump window (kb)", value = 10000, min = 250, max = 50000, step = 250)),
            div(class = "sep"),
            div(class = "panel",
                div(class = "panel-title", "Outcome extraction options"),
                checkboxInput("use_proxies", "Use proxies if SNP missing", value = TRUE),
                sliderInput("proxy_rsq",  "Min proxy r\u00b2", min = 0.1, max = 1.0, value = 0.8, step = 0.05),
                checkboxInput("allow_palindromes", "Allow palindromic SNPs", value = TRUE),
                sliderInput("maf_thresh", "MAF threshold", min = 0.0, max = 0.5, value = 0.30, step = 0.01),
                radioButtons("harm_action", "Harmonise action",
                             choices = c("Positive strand" = "pos",
                                         "Align (recommended)" = "align",
                                         "Exclude palindromes" = "exclude"),
                             selected = "align")),
            div(class = "sep"),
            div(class = "panel",
                div(class = "panel-title", "MR methods"),
                checkboxGroupInput("mr_methods", "Methods",
                                   choices  = mr_method_choices(),
                                   selected = c("mr_wald_ratio","mr_ivw",
                                                "mr_egger_regression",
                                                "mr_weighted_median","mr_weighted_mode")),
                actionButton("run_mr", "\u26a1 Run MR", class = "btn-primary w100"),
                div(class = "sep"),
                div(class = "panel-title", "Downloads"),
                conditionalPanel(
                  condition = "output.tab3_has_results === true",
                  downloadButton("dl_tab3_mr_xlsx",   "\u2b07 MR Results (Excel)",           class = "btn-ghost w100"), div(class="sep"),
                  downloadButton("dl_tab3_sens_xlsx",  "\u2b07 Sensitivity Tables (Excel)",   class = "btn-ghost w100"), div(class="sep"),
                  downloadButton("dl_plot_scatter",    "\u2b07 Method comparison plot (PDF)",  class = "btn-ghost w100"), div(class="sep"),
                  downloadButton("dl_plot_forest",     "\u2b07 Single SNP analysis plot (PDF)",class = "btn-ghost w100"), div(class="sep"),
                  downloadButton("dl_plot_loo",        "\u2b07 Leave-one-out plot (PDF)",      class = "btn-ghost w100"), div(class="sep"),
                  downloadButton("dl_plot_funnel",     "\u2b07 Funnel plot (PDF)",             class = "btn-ghost w100"), div(class="sep"),
                  downloadButton("dl_tab3_snp_xlsx",   "\u2b07 SNP Effects Table (Excel)",     class = "btn-ghost w100")),
                conditionalPanel(
                  condition = "output.tab3_has_results === false",
                  tags$p(style="opacity:0.6;font-style:italic;",
                         "Downloads appear after MR results are available.")))
          ),
          # ── Right column: status + results ─────────────────
          column(
            8,
            div(class = "panel",
                div(class = "panel-title", "Status"),
                verbatimTextOutput("mr_status"),
                tags$hr(),
                verbatimTextOutput("opengwas_countdown")),
            div(class = "sep"),
            div(class = "panel",
                div(class = "panel-title", "Choose pair"),
                uiOutput("pair_picker_ui")),
            div(class = "sep"),
            div(
              class = "panel",
              div(class = "panel-title", "MR Results & Statistics"),
              conditionalPanel(
                condition = "output.tab3_has_results === false",
                tags$p(style="opacity:0.7;font-style:italic;",
                       "MR Results & Statistics will appear here after running the analysis.")),
              conditionalPanel(
                condition = "output.tab3_has_results === true",
                div(style="display:flex;align-items:center;justify-content:space-between;",
                    tags$h4("MR Results"),
                    actionButton("info_mr_results", label=NULL, icon=icon("info-circle"), class="btn-link")),
                div(class="note","The MR Results table summarizes causal effect estimates. The null hypothesis (H\u2080) states no causal effect (\u03b2\u202f=\u202f0; OR\u202f=\u202f1). A small p-value (< 0.05) supports rejection of H\u2080."),
                DT::DTOutput("tbl_mr_results"),
                div(style="display:flex;align-items:center;justify-content:space-between;",
                    tags$h4("Heterogeneity statistics"),
                    actionButton("info_het", label=NULL, icon=icon("info-circle"), class="btn-link")),
                div(class="note","Cochran\u2019s Q tests whether all instruments estimate a common causal effect. A significant Q p-value (< 0.05) suggests heterogeneity, possibly from horizontal pleiotropy."),
                DT::DTOutput("tbl_het"),
                div(style="display:flex;align-items:center;justify-content:space-between;",
                    tags$h4("Horizontal pleiotropy"),
                    actionButton("info_pleio", label=NULL, icon=icon("info-circle"), class="btn-link")),
                div(class="note","The MR-Egger intercept estimates average directional pleiotropy. A significant intercept (p < 0.05) suggests instruments affect the outcome through pathways other than the exposure."),
                verbatimTextOutput("pleio_summary"),
                DT::DTOutput("tbl_pleio"),
                div(style="display:flex;align-items:center;justify-content:space-between;",
                    tags$h4("Causal direction test (Steiger)"),
                    actionButton("info_dir", label=NULL, icon=icon("info-circle"), class="btn-link")),
                div(class="note","Tests whether instruments explain more variance in the exposure than the outcome. A small p-value supports the specified direction (Exposure \u2192 Outcome)."),
                uiOutput("steiger_binary_msg"),
                DT::DTOutput("tbl_dir"))),
            div(class = "sep"),
            div(
              class = "panel",
              div(class = "panel-title", "Plots"),
              conditionalPanel(
                condition = "output.tab3_has_results === false",
                tags$p(style="opacity:0.7;font-style:italic;margin-top:8px;",
                       "Plots appear after running the MR analysis.")),
              conditionalPanel(
                condition = "output.tab3_has_results === true",
                tags$h4("Method comparison plot"),
                plotOutput("plot_scatter", height="380px"),
                div(class="note","SNP effects on the outcome vs. exposure. Slope = causal estimate per method."),
                tags$h4("Single SNP analysis"),
                plotOutput("plot_forest",  height="380px"),
                div(class="note","Wald ratio estimates per SNP shown as forest plot."),
                tags$h4("Leave-one-out analysis"),
                plotOutput("plot_loo",     height="380px"),
                div(class="note","IVW estimate recalculated excluding one SNP at a time."),
                tags$h4("Funnel plot"),
                plotOutput("plot_funnel",  height="380px"),
                div(class="note","Asymmetry suggests directional pleiotropy.")))
          )
        )
      )
    )
  )
}

# ── Tab 4: Multivariable MR ──────────────────────────────────
create_multivariable_mr_tab <- function() {
  tabPanel(
    "Multivariable MR",
    div(
      class = "card",
      div(class = "card-header",
          div(class = "card-title", "Multivariable Mendelian Randomization Workflow"),
          div(class = "card-subtitle",
              "Workflow: Identify repeated-SNP studies \u2192 Select Exposure 2 / Mediator \u2192 Run Exposure\u202f1\u2192Mediator MR \u2192 Build SNP table \u2192 Run MVMR.")),
      div(
        class = "card-body",
        tabsetPanel(

          # ── A) Repeated SNP Independent Studies ──────────
          tabPanel(
            "A) Repeated SNP Independent Studies",
            div(
              class = "card",
              div(class = "card-header",
                  div(class = "card-title", "Repeated SNP studies (from Tab 3 SNP table)"),
                  div(class = "card-subtitle",
                      "Studies sharing SNP instruments with the selected Exposure 1 \u2192 Outcome pair. Use the filters below to restrict to same-ancestry candidates.")),
              div(
                class = "card-body",
                div(
                  class = "panel", style = "margin-bottom:8px;",
                  div(class = "panel-title", "Candidate study filter options"),

                  # First checkbox: same-ancestry filter
                  checkboxInput("rep_only_indep",
                                "Show only same-ancestry repeated-SNP studies (recommended)",
                                value = FALSE),

                  # Second checkbox: subgroup-independent filter
                  # Shown only when first checkbox is ON and exposure has specific subgroup resolution
                  # (output$has_cohort_exp reflects subgroup resolution, not cohort presence)
                  conditionalPanel(
                    condition = "input.rep_only_indep == true && output.has_cohort_exp == 'TRUE'",
                    checkboxInput("rep_non_overlap",
                                  "Within same ancestry, show only sub-group independent studies (vs. both Exposure\u202f1 and Outcome)",
                                  value = FALSE)),

                  # Contextual message (server-rendered)
                  uiOutput("rep_ancestry_msg")
                ),
                div(class = "sep"),
                DT::DTOutput("tbl_repeated")
              )
            )
          ),

          # ── B) Select Exposure 2 / Mediator ──────────────
          tabPanel(
            "B) Select Exposure 2 / Mediator",
            div(class = "sep"),
            div(
              class = "card",
              div(class = "card-header",
                  div(class = "card-title", "Select Exposure 2 / Mediator"),
                  div(class = "card-subtitle",
                      "Choose a mediator from repeated-SNP candidates (Tab A) or search for a new one. The mediator must share ancestry with both Exposure 1 and Outcome for valid MVMR.")),
              div(
                class = "card-body",
                div(class = "panel",
                    div(class = "panel-title", "Auto-filled Exposure 1 (from Tab 1)"),
                    verbatimTextOutput("exp1_auto_box")),
                div(class = "sep"),
                div(
                  class = "card-body",
                  radioButtons("mediator_source", "Mediator source",
                               choices  = c("Use Exposure 2 (from Tab A: Repeated SNP Studies)" = "from_tabA",
                                            "Search a New Exposure 2 / Mediator"                 = "new_search"),
                               selected = "from_tabA"),

                  # New mediator search
                  conditionalPanel(
                    condition = "input.mediator_source == 'new_search'",
                    div(class = "sep"),
                    searchUI(
                      "tab7_mediator_search",
                      "Search & Select NEW Mediator",
                      show_selection_box = TRUE,
                      selection_title    = "Selected Mediator GCST",
                      extra_ui = div(
                        class = "panel", style = "margin-bottom:8px;",
                        div(class = "panel-title",
                            "Mediator independence options (vs Exposure 1 & Outcome)"),

                        # First checkbox: same-ancestry mediator filter
                        checkboxInput(
                          "tab7_only_indep",
                          paste0("Show only same-ancestry Mediator studies (recommended)",
                                 " \u2014 required for valid MVMR estimation"),
                          value = FALSE),

                        # Second checkbox: subgroup-independent filter
                        # Rendered by server — appears only when first is ON and ancestry
                        # has specific subgroup resolution
                        uiOutput("tab7_subgroup_checkbox_ui"),

                        # Contextual message (server-rendered)
                        uiOutput("tab7_ancestry_msg")
                      )
                    )
                  ),

                  div(class = "sep"),

                  # Pick from Tab A
                  conditionalPanel(
                    condition = "input.mediator_source == 'from_tabA'",
                    div(class = "panel",
                        div(class = "panel-title", "Pick Exposure 2 from repeated-SNP studies"),
                        uiOutput("exp2_pick_ui"),
                        div(class = "sep"),
                        div(class = "label", "Currently selected Exposure 2"),
                        verbatimTextOutput("exp2_pick_box"))
                  ),
                  div(class = "sep")
                )
              )
            )
          ),

          # ── C) Run Exposure 1 → Mediator MR ──────────────
          tabPanel(
            "C) Run Exposure 1 \u2192 Exposure 2 / Mediator MR",
            div(
              class = "card",
              div(class = "card-header",
                  div(class = "card-title",
                      "Two-sample MR: Exposure 1 \u2192 Exposure 2 / Mediator"),
                  div(class = "card-subtitle",
                      "Estimate the causal effect of Exposure 1 on the Mediator (X \u2192 M pathway).")),
              div(
                class = "card-body",
                fluidRow(
                  column(
                    4,
                    div(class = "panel",
                        div(class = "panel-title", "Instrument settings"),
                        checkboxInput("tab5_do_clump", "Use LD clumping (OpenGWAS)", value = TRUE),
                        numericInput("tab5_p_thresh",  "Instrument p-value (p1)", value = 5e-8, min = 0, max = 1, step = 1e-8),
                        sliderInput("tab5_clump_r2", "Clump r\u00b2", min = 0.0001, max = 0.2, value = 0.01, step = 0.0005),
                        numericInput("tab5_clump_kb", "Clump window (kb)", value = 10000, min = 250, max = 50000, step = 250)),
                    div(class = "sep"),
                    div(class = "panel",
                        div(class = "panel-title", "Outcome extraction options"),
                        checkboxInput("tab5_use_proxies", "Use proxies if SNP missing", value = TRUE),
                        sliderInput("tab5_proxy_rsq",  "Min proxy r\u00b2", min = 0.1, max = 1.0, value = 0.8, step = 0.05),
                        checkboxInput("tab5_allow_palindromes", "Allow palindromic SNPs", value = TRUE),
                        sliderInput("tab5_maf_thresh", "MAF threshold", min = 0.0, max = 0.5, value = 0.30, step = 0.01),
                        radioButtons("tab5_harm_action", "Harmonise action",
                                     choices  = c("Positive strand" = "pos",
                                                  "Align (recommended)" = "align",
                                                  "Exclude palindromes" = "exclude"),
                                     selected = "align")),
                    div(class = "sep"),
                    div(class = "panel",
                        div(class = "panel-title", "MR methods"),
                        checkboxGroupInput("tab5_mr_methods", "Methods",
                                           choices  = mr_method_choices(),
                                           selected = c("mr_wald_ratio","mr_ivw",
                                                        "mr_egger_regression",
                                                        "mr_weighted_median","mr_weighted_mode")),
                        actionButton("run_mr_tab5",
                                     "\u26a1 Run MR (Exposure\u202f1\u2192Mediator)",
                                     class = "btn-primary w100"),
                        div(class = "sep"),
                        div(class = "panel-title", "Downloads"),
                        conditionalPanel(
                          condition = "output.tab5_has_results === true",
                          downloadButton("dl_tab5_mr_xlsx",       "\u2b07 MR Results (Excel)",            class="btn-ghost w100"), div(class="sep"),
                          downloadButton("dl_tab5_sens_xlsx",      "\u2b07 Sensitivity Tables (Excel)",    class="btn-ghost w100"), div(class="sep"),
                          downloadButton("dl_tab5_plot_scatter",   "\u2b07 Method comparison plot (PDF)",  class="btn-ghost w100"), div(class="sep"),
                          downloadButton("dl_tab5_plot_forest",    "\u2b07 Single SNP analysis plot (PDF)",class="btn-ghost w100"), div(class="sep"),
                          downloadButton("dl_tab5_plot_loo",       "\u2b07 Leave-one-out plot (PDF)",      class="btn-ghost w100"), div(class="sep"),
                          downloadButton("dl_tab5_plot_funnel",    "\u2b07 Funnel plot (PDF)",             class="btn-ghost w100"), div(class="sep"),
                          downloadButton("dl_tab5_snp_table_xlsx", "\u2b07 SNP Effects Table (Excel)",     class="btn-ghost w100")),
                        conditionalPanel(
                          condition = "output.tab5_has_results === false",
                          tags$p(style="opacity:0.6;font-style:italic;",
                                 "Downloads appear after running the MR analysis.")))
                  ),
                  column(
                    8,
                    div(class = "sep"),
                    div(class = "panel",
                        div(class = "panel-title", "Status"),
                        verbatimTextOutput("tab5_mr_status"),
                        tags$hr(),
                        verbatimTextOutput("opengwas_countdown")),
                    div(class = "panel",
                        div(class = "panel-title", "Choose pair (X \u2192 M)"),
                        uiOutput("tab5_pair_ui"),
                        conditionalPanel(
                          condition = "output.tab5_has_results === false",
                          tags$p(style="margin-top:10px;opacity:0.7;font-style:italic;",
                                 "No MR results yet."))),
                    div(class = "sep"),
                    div(
                      class = "panel",
                      div(class = "panel-title", "MR Results & Statistics"),
                      conditionalPanel(
                        condition = "output.tab5_has_results === false",
                        tags$p(style="opacity:0.7;font-style:italic;",
                               "MR Results & Statistics will appear here after running the analysis.")),
                      conditionalPanel(
                        condition = "output.tab5_has_results === true",
                        div(style="display:flex;align-items:center;justify-content:space-between;",
                            tags$h4("MR Results"),
                            actionButton("info_tab5_mr", label=NULL, icon=icon("info-circle"), class="btn-link")),
                        div(class="note","Causal effect of Exposure\u202f1 on Mediator (X\u202f\u2192\u202fM). For binary mediators, effects are in log-odds scale."),
                        DT::DTOutput("tbl_tab5_mr_results"),
                        div(class = "sep"),
                        div(style="display:flex;align-items:center;justify-content:space-between;",
                            tags$h4("Heterogeneity statistics"),
                            actionButton("info_tab5_het", label=NULL, icon=icon("info-circle"), class="btn-link")),
                        div(class="note","Cochran\u2019s Q test for instrument heterogeneity."),
                        DT::DTOutput("tbl_tab5_het"),
                        div(class = "sep"),
                        div(style="display:flex;align-items:center;justify-content:space-between;",
                            tags$h4("Horizontal pleiotropy"),
                            actionButton("info_tab5_pleio", label=NULL, icon=icon("info-circle"), class="btn-link")),
                        div(class="note","MR-Egger intercept test for directional pleiotropy."),
                        DT::DTOutput("tbl_tab5_pleio"),
                        div(class = "sep"),
                        div(style="display:flex;align-items:center;justify-content:space-between;",
                            tags$h4("Causal direction test (Steiger)"),
                            actionButton("info_tab5_dir", label=NULL, icon=icon("info-circle"), class="btn-link")),
                        div(class="note","Tests whether instruments explain more variance in Exposure\u202f1 than in the Mediator."),
                        uiOutput("steiger_binary_msg_tab5"),
                        DT::DTOutput("tbl_tab5_dir"))),
                    div(class = "sep"),
                    div(
                      class = "panel",
                      div(class = "panel-title", "Plots"),
                      conditionalPanel(
                        condition = "output.tab5_has_results === true",
                        tags$h4("Method comparison plot"),   plotOutput("tab5_plot_scatter", height="360px"),
                        tags$h4("Single SNP analysis"),      plotOutput("tab5_plot_forest",  height="360px"),
                        tags$h4("Leave-one-out analysis"),   plotOutput("tab5_plot_loo",     height="360px"),
                        tags$h4("Funnel plot"),              plotOutput("tab5_plot_funnel",  height="360px")),
                      conditionalPanel(
                        condition = "output.tab5_has_results === false",
                        tags$p(style="opacity:0.7;font-style:italic;",
                               "Plots appear after running the MR analysis.")))
                  )
                )
              )
            )
          ),

          # ── D) Build SNP Table ─────────────────────────────
          tabPanel(
            "D) Build SNP Table",
            div(
              class = "card",
              div(class = "card-header",
                  div(class = "card-title", "Build SNP Table for MVMR"),
                  div(class = "card-subtitle",
                      "Combines SNP instruments from Exposure 1 and Exposure 2 / Mediator into a single MVMR SNP effects table.")),
              div(
                class = "card-body",
                div(
                  class = "panel",
                  actionButton("build_tab5", "\U0001f9e9 Build SNP Table", class = "btn-primary w100"),
                  div(class = "sep"),
                  conditionalPanel(
                    condition = "output.tab5_snp_ready === true",
                    downloadButton("dl_tab5_snp_mvmr_xlsx", "\u2b07 Download SNP Table (Excel)",
                                   class = "btn-ghost w100")),
                  conditionalPanel(
                    condition = "output.tab5_snp_ready === false",
                    tags$p(style="opacity:0.6;font-style:italic;",
                           "Download available after building SNP table.")),
                  div(class = "sep"),
                  verbatimTextOutput("tab5_warn"),
                  DT::DTOutput("tbl_tab5_snp_mvmr")
                )
              )
            )
          ),

          # ── E) Run MVMR ────────────────────────────────────
          tabPanel(
            "E) Run MVMR",
            div(
              class = "card",
              div(class = "card-header",
                  div(class = "card-title", "Multivariable MR"),
                  div(class = "card-subtitle",
                      "MVMR estimates the direct causal effects of Exposure\u202f1 and Exposure\u202f2\u202f/\u202fMediator on the Outcome, adjusted for each other.")),
              div(
                class = "card-body",
                fluidRow(
                  column(
                    4,
                    div(class = "panel",
                        actionButton("run_mvmr", "\u26a1 Run MVMR", class = "btn-primary w100")),
                    div(class = "sep"),
                    div(class = "panel",
                        div(class = "panel-title", "Downloads"),
                        conditionalPanel(
                          condition = "output.mvmr_has_results === true",
                          downloadButton("dl_mvmr_full_xlsx", "\u2b07 Download Full MVMR Results",
                                         class = "btn-ghost w100")),
                        conditionalPanel(
                          condition = "output.mvmr_has_results === false",
                          tags$p(style="opacity:0.7;font-style:italic;",
                                 "Download appears after MVMR results are available.")))
                  ),
                  column(
                    8,
                    div(class = "panel",
                        div(class = "panel-title", "Status"),
                        verbatimTextOutput("mvmr_status")),
                    div(class = "sep"),
                    div(
                      class = "panel",
                      div(class = "panel-title", "MVMR Results & Statistics"),
                      conditionalPanel(
                        condition = "output.mvmr_has_results === false",
                        tags$p(style="opacity:0.7;font-style:italic;",
                               "MVMR Results will appear here after running the analysis.")),
                      conditionalPanel(
                        condition = "output.mvmr_has_results === true",
                        tags$h4("MVMR IVW result"),
                        div(class="note","Direct causal effects adjusted for other exposures. H\u2080: \u03b2\u202f=\u202f0."),
                        DT::DTOutput("tbl_mvmr_ivw"),
                        div(class = "sep"),
                        tags$h4("Instrument strength"),
                        div(class="note","Conditional F-statistics: F\u202f>\u202f10 suggests low weak-instrument bias."),
                        DT::DTOutput("tbl_mvmr_strength"),
                        div(class = "sep"),
                        tags$h4("Test for horizontal pleiotropy"),
                        div(class="note","Heterogeneity Q-statistic. Significant result (p\u202f<\u202f0.05) may indicate pleiotropy."),
                        DT::DTOutput("tbl_mvmr_pleio")))
                  )
                )
              )
            )
          )
        )
      )
    )
  )
}

# ── Tab 5: Mediation Analysis ────────────────────────────────
create_mediation_tab <- function() {
  tabPanel(
    "Mediation Analysis",
    div(
      class = "card",
      div(class = "card-header",
          div(class = "card-title", "Mediation Analysis"),
          div(class = "card-subtitle",
              "Decomposes the total causal effect (X\u202f\u2192\u202fY) into a direct effect (X\u202f\u2192\u202fY\u202f|\u202fM) and an indirect effect mediated through M, using both difference and product methods.")),
      div(
        class = "card-body",
        fluidRow(
          column(12, uiOutput("ui_tab7_main_tabs"))
        )
      )
    )
  )
}

# ============================================================
# Assemble the UI
# ============================================================

build_ui <- function() {
  fluidPage(
    tags$head(
      tags$style(HTML(custom_css)),
      custom_scroll_script(),
      tags$script(HTML(
        "Shiny.addCustomMessageHandler('playSound',function(msg){
           var a=new Audio('https://actions.google.com/sounds/v1/alarms/beep_short.ogg');
           a.play();
         });"
      ))
    ),
    # Show a banner if data failed to load
    if (!is.null(DATA$load_error)) {
      div(style = "background:#7f1d1d;color:#fca5a5;padding:12px 18px;border-radius:8px;margin:10px 0;font-weight:700;",
          paste0("\u26a0\ufe0f Data load error: ", DATA$load_error,
                 " \u2014 Please check your data/ folder and restart."))
    },
    create_header_card(),
    tabsetPanel(
      id = "main_tabs",
      create_select_exposure_tab(),
      create_select_outcome_tab(),
      create_two_sample_mr_tab(),
      create_multivariable_mr_tab(),
      create_mediation_tab()
    )
  )
}

# Build and register the UI
ui <- build_ui()

# Source server (defines `server` function)
source(file.path(APP_DIR, "R", "server_original.R"), local = FALSE)

shinyApp(ui, server)
