# IndepMR — Main Application

library(shiny)

##----------------------------------------------------------------------------
##  Resource Path and Configuration
##----------------------------------------------------------------------------

# Path to the directory containing image assets.  Register this path under
# the `assets` alias so that images can be referenced via `src="assets/..."`.
#img_dir <- "/Users/thenukayatawara/Desktop/IndepMR"
#addResourcePath("assets", img_dir)
tags$img(src = "logo.png")
##----------------------------------------------------------------------------
##  Custom Stylesheet and JavaScript
##----------------------------------------------------------------------------

# Consolidate all of the custom CSS rules used in the app into a single
# character string.  This makes it easy to maintain and reason about the
# overall look and feel of the interface.
custom_css <- "
:root {
  --bg1: #0b1220;
  --card: #0f1b33;
  --panel: #0b1730;
  --text: #eaf0ff;
  --muted: #a9b4c7;
  --line: rgba(255, 255, 255, 0.10);
  --accent: #4da3ff;
  --accent2: #6d72ff;
  --ok: #34d399;
  --bad: #fb7185;
  --unk: #fbbf24;
}
body {
  background: radial-gradient(1200px 750px at 18% 12%, #163a7a 0%, var(--bg1) 45%, #070b14 100%);
  color: var(--text);
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Arial, sans-serif;
}
.container-fluid { max-width: 1350px; }
.app-header {
  margin-top: 12px;
  border-radius: 18px;
  padding: 18px;
  background: linear-gradient(135deg, rgba(77, 163, 255, 0.22), rgba(109, 114, 255, 0.14));
  border: 1px solid rgba(255, 255, 255, 0.12);
  box-shadow: 0 14px 40px rgba(0, 0, 0, 0.35);
}
.app-title { font-size: 26px; font-weight: 900; margin: 0; }
.app-subtitle { margin-top: 6px; color: var(--muted); font-size: 13px; line-height: 1.35; }
.card {
  margin-top: 14px;
  border-radius: 18px;
  background: rgba(15, 27, 51, 0.86);
  border: 1px solid rgba(255, 255, 255, 0.12);
  box-shadow: 0 12px 36px rgba(0, 0, 0, 0.32);
  overflow: hidden;
}
.card-header {
  padding: 14px 16px;
  background: linear-gradient(180deg, rgba(255, 255, 255, 0.06), rgba(255, 255, 255, 0.02));
  border-bottom: 1px solid var(--line);
}
.card-title { font-size: 18px; font-weight: 900; margin: 0; }
.card-subtitle { margin-top: 6px; color: var(--muted); font-size: 12px; }
.card-body { padding: 14px 16px 16px 16px; }
.label { color: var(--muted); font-size: 12px; margin-bottom: 6px; }
.sep { height: 14px; }
.panel {
  border-radius: 16px;
  background: rgba(11, 23, 48, 0.62);
  border: 1px solid rgba(255, 255, 255, 0.10);
  padding: 12px;
}
.panel-title { font-size: 13px; font-weight: 900; margin-bottom: 10px; }
pre {
  color: var(--text);
  background: rgba(0, 0, 0, 0.18);
  border: 1px solid rgba(255, 255, 255, 0.10);
  border-radius: 12px;
  padding: 10px;
  white-space: pre-wrap;
}
.form-control, .selectize-input {
  background: rgba(8, 16, 34, 0.70) !important;
  color: var(--text) !important;
  border: 1px solid rgba(255, 255, 255, 0.14) !important;
  border-radius: 12px !important;
}
.selectize-dropdown, .selectize-dropdown-content {
  background: #081022 !important;
  color: var(--text) !important;
  border: 1px solid rgba(255, 255, 255, 0.14) !important;
}
.btn-primary {
  border: 0 !important;
  border-radius: 12px !important;
  height: 40px !important;
  background: linear-gradient(90deg, rgba(77, 163, 255, 0.95), rgba(109, 114, 255, 0.95)) !important;
  color: white !important;
  font-weight: 800 !important;
}
.btn-ghost {
  border-radius: 12px !important;
  height: 40px !important;
  background: rgba(255, 255, 255, 0.06) !important;
  color: var(--text) !important;
  font-weight: 800 !important;
  border: 1px solid rgba(255, 255, 255, 0.12) !important;
}
.w100 { width: 100%; }
.nav-tabs {
  border-bottom: 1px solid rgba(255, 255, 255, 0.10);
  margin-top: 12px;
}
.nav-tabs > li > a {
  border-radius: 14px 14px 0 0 !important;
  background: rgba(255, 255, 255, 0.05) !important;
  border: 1px solid rgba(255, 255, 255, 0.10) !important;
  color: var(--muted) !important;
  font-weight: 800 !important;
  margin-right: 8px !important;
}
.nav-tabs > li.active > a {
  color: var(--text) !important;
  background: linear-gradient(180deg, rgba(77, 163, 255, 0.20), rgba(255, 255, 255, 0.05)) !important;
  border-color: rgba(77, 163, 255, 0.35) !important;
}
table.dataTable {
  border-collapse: separate !important;
  border-spacing: 0 8px !important;
}
table.dataTable thead th {
  background: rgba(255, 255, 255, 0.06) !important;
  color: var(--text) !important;
  border: none !important;
}
table.dataTable tbody tr {
  background: rgba(255, 255, 255, 0.04) !important;
  border: 1px solid rgba(255, 255, 255, 0.10) !important;
  border-radius: 12px !important;
}
table.dataTable tbody td {
  border: none !important;
  color: var(--text) !important;
}
.dataTables_info, .dataTables_paginate, .dataTables_length, .dataTables_filter {
  color: var(--muted) !important;
}
.note {
  color: var(--muted);
  font-size: 12px;
  line-height: 1.35;
}
.spinWrap {
  display: flex;
  align-items: center;
  gap: 12px;
}
.spinner {
  width: 22px;
  height: 22px;
  border: 3px solid rgba(255, 255, 255, 0.25);
  border-top: 3px solid rgba(77, 163, 255, 0.95);
  border-radius: 50%;
  animation: spin 0.9s linear infinite;
}
@keyframes spin { to { transform: rotate(360deg); } }
.modal-content {
  background: rgba(15, 27, 51, 0.96) !important;
  color: var(--text) !important;
  border: 1px solid rgba(255, 255, 255, 0.14) !important;
  border-radius: 16px !important;
}
/* Custom header styling */
.header-card {
  display: flex;
  align-items: center;
  border: 1px solid rgba(255, 255, 255, 0.35);
  border-radius: 16px;
  overflow: hidden;
  margin-top: 16px;
}
.header-card .logo-col {
  background: #ffffff;
  padding: 20px;
  display: flex;
  align-items: center;
  justify-content: center;
  min-width: 180px;
}
.header-card .text-col {
  background: rgba(11, 23, 48, 0.62);
  color: var(--text);
  padding: 20px;
  flex: 1;
  display: flex;
  flex-direction: column;
  justify-content: center;
}
.header-title {
  font-size: 22px;
  font-weight: 900;
  margin: 0;
  line-height: 1.2;
}
.header-subtitle {
  margin-top: 8px;
  font-size: 12px;
  line-height: 1.35;
  color: var(--muted);
}
.tab7-box {
  border: 2px solid #1f1f1f;
  border-radius: 8px;
  padding: 14px;
  background: #ffffff;
  min-height: 240px;
}
.tab7-runbox {
  display: flex;
  justify-content: center;
  align-items: center;
  min-height: 240px;
}
.btn-runmodel {
  font-size: 22px;
  font-weight: 700;
  padding: 18px 38px;
  border-radius: 14px;
}
"

# JavaScript helper for smooth scrolling to an element when triggered from the
# server.  This is registered as a custom message handler via
# `session$sendCustomMessage('scrollTo', list(id = ...))`.
custom_scroll_script <- function() {
  tags$script(HTML(
    "
    Shiny.addCustomMessageHandler('scrollTo', function(message) {
      var el = document.getElementById(message.id);
      if (el) el.scrollIntoView({ behavior: 'smooth', block: 'start' });
    });
    "
  ))
}

##----------------------------------------------------------------------------
##  UI Helper Functions
##----------------------------------------------------------------------------

# Header card with the logo and a concise description of the application.
create_header_card <- function() {
  div(
    class = "app-header",
    div(
      class = "header-card",
      div(
        class = "logo-col",
        tags$img(
          src = "assets/logo.png",
          style = "height: 120px; width: auto;"
        )
      ),
      div(
        class = "text-col",
        h1(class = "header-title",
           "An Ontology-Driven and Ancestry-Aware Automated Framework for Two-Sample, Multivariable, and Mediation Mendelian Randomization Using Public GWAS Summary Data"
        ),
        div(
          class = "header-subtitle",
          "IndepMR is an automated framework for Two-Sample MR, Multivariable MR, and Mediation Analysis using public GWAS summary data. It integrates GWAS Catalog and OpenGWAS resources to automatically identify ancestry-independent studies and avoid sample overlap, enabling valid causal inference without manual data curation.The platform supports multiple MR methods, sensitivity diagnostics, MVMR for direct effect estimation, and automated mediation analysis using both product and difference approaches. IndepMR provides a streamlined, reproducible solution for estimating total, direct, and indirect causal effects without requiring individual-level data."
        )
      )
    )
  )
}

# Tab 1: Select Exposure 1
create_select_exposure_tab <- function() {
  tabPanel(
    "Select Exposure 1",
    div(
      searchUI(
        "exposure",
        "Select Exposure 1",
        show_selection_box = TRUE,
        selection_title = "Selected Exposure 1 GCST IDs"
      ),
      
      div(class = "sep"),
      
      tags$a(
        href = "https://www.ebi.ac.uk/gwas/",
        target = "_blank",
        class = "btn btn-ghost",
        "More study details (GWAS Catalog)"
      )
    )
  )
}

# Tab 2: Select Outcome
create_select_outcome_tab <- function() {
  tabPanel(
    "Select Outcome",
    div(
      searchUI(
        "outcome",
        "Select Outcome",
        show_selection_box = TRUE,
        selection_title = "Selected Outcome GCST IDs",
        extra_ui = div(
          class = "panel",
          div(class = "panel-title", "Outcome independence options (vs Exposure 1)"),
          # First checkbox: restrict to same‑ancestry (not independent) outcome studies
          checkboxInput(
            "show_only_independent",
            "Show only same-ancestry Outcome studies (recommended)",
            value = TRUE
          ),
          # Second checkbox appears only when the first is checked AND cohort
          # information exists for both exposure and outcome.  It restricts to
          # non‑overlapping cohorts.
          conditionalPanel(
            condition = "input.show_only_independent && output.has_cohort_exp",
            checkboxInput(
              "show_non_overlap_outcome",
              "Within same ancestry, show only non-overlapping cohorts (to Avoid Sample Overlap)",
              value = FALSE
            )
          )
        )
      ),
      div(class = "sep"),
      tags$a(
        href = "https://www.ebi.ac.uk/gwas/",
        target = "_blank",
        class = "btn btn-ghost",
        "More study details (GWAS Catalog)"
      )
    )
  )
}

# Tab 3: Two-sample Mendelian Randomization
create_two_sample_mr_tab <- function() {
  tabPanel(
    "Two-sample MR",
    div(
      class = "card",
      div(
        class = "card-header",
        div(class = "card-title", "Two-sample MR"),
        div(
          class = "card-subtitle",
          "Users can adjust instrument selection parameters (p-value threshold and LD clumping settings), choose MR estimation methods, and evaluate causal effects using comprehensive sensitivity analyses and diagnostic plots."
        )
      ),
      div(
        class = "card-body",
        fluidRow(
          column(
            4,
            # Instrument settings
            div(
              class = "panel",
              div(class = "panel-title", "Instrument settings"),
              checkboxInput("do_clump", "Use LD clumping (OpenGWAS)", value = TRUE),
              numericInput("p_thresh", "Instrument p-value (p1)", value = 5e-8, min = 0, max = 1, step = 1e-8),
              sliderInput("clump_r2", "Clump r²", min = 0.0001, max = 0.2, value = 0.01, step = 0.0005),
              numericInput("clump_kb", "Clump window (kb)", value = 10000, min = 250, max = 50000, step = 250)
            ),
            div(class = "sep"),
            # Outcome extraction options
            div(
              class = "panel",
              div(class = "panel-title", "Outcome extraction options"),
              checkboxInput("use_proxies", "Use proxies if SNP missing", value = TRUE),
              sliderInput("proxy_rsq", "Min proxy r²", min = 0.1, max = 1.0, value = 0.8, step = 0.05),
              checkboxInput("allow_palindromes", "Allow palindromic SNPs", value = TRUE),
              sliderInput("maf_thresh", "MAF threshold", min = 0.0, max = 0.5, value = 0.30, step = 0.01),
              radioButtons(
                "harm_action",
                "Harmonise action",
                choices = c("Positive strand" = "pos", "Align (recommended)" = "align", "Exclude palindromes" = "exclude"),
                selected = "align"
              )
            ),
            div(class = "sep"),
            # MR methods and run button
            div(
              class = "panel",
              div(class = "panel-title", "MR methods"),
              checkboxGroupInput(
                "mr_methods",
                "Methods",
                choices = mr_method_choices(),
                selected = c(
                  "mr_wald_ratio",
                  "mr_ivw",
                  "mr_egger_regression",
                  "mr_weighted_median",
                  "mr_weighted_mode"
                )
              ),
              actionButton("run_mr", "⚡ Run MR", class = "btn-primary w100"),
              div(class = "sep"),
              div(class = "panel-title", "Downloads (current selected pair)"),
              conditionalPanel(
                condition = "output.tab3_has_results === true",
                downloadButton("dl_tab3_mr_xlsx", "⬇ MR Results (Excel)", class = "btn-ghost w100"),
                div(class = "sep"),
                downloadButton("dl_tab3_sens_xlsx", "⬇ Sensitivity Tables (Excel)", class = "btn-ghost w100"),
                div(class = "sep"),
                downloadButton("dl_plot_scatter", "⬇ Method comparison plot (PDF)", class = "btn-ghost w100"),
                div(class = "sep"),
                downloadButton("dl_plot_forest", "⬇ Single SNP analysis plot (PDF)", class = "btn-ghost w100"),
                div(class = "sep"),
                downloadButton("dl_plot_loo", "⬇ Leave-one-out analysis plot (PDF)", class = "btn-ghost w100"),
                div(class = "sep"),
                downloadButton("dl_plot_funnel", "⬇ Funnel plot (PDF)", class = "btn-ghost w100"),
                div(class = "sep"),
                downloadButton("dl_tab3_snp_xlsx", "⬇ SNP Effects Table (Excel)", class = "btn-ghost w100")
              ),
              conditionalPanel(
                condition = "output.tab3_has_results === false",
                tags$p(
                  style = "opacity:0.6; font-style:italic;",
                  "Downloads will appear after MR results are available."
                )
              )
            )
          ),
          column(
            8,
            # Status and pair picker
            div(
              class = "panel",
              div(class = "panel-title", "Status"),
              verbatimTextOutput("mr_status"),
              tags$hr(),
              verbatimTextOutput("opengwas_countdown")
            ),
            div(class = "sep"),
            div(class = "panel", div(class = "panel-title", "Choose pair"), uiOutput("pair_picker_ui")),
            div(class = "sep"),
            # MR results and stats
            div(
              class = "panel",
              div(class = "panel-title", "MR Results & Statistics"),
              
              # 🔹 Message before MR run
              conditionalPanel(
                condition = "output.tab3_has_results === false",
                tags$p(
                  style = "opacity:0.7; font-style:italic;",
                  "MR Results & Statistics will appear here after running the MR analysis."
                )
              ),
              
              # 🔹 Actual tables after MR run
              conditionalPanel(
                condition = "output.tab3_has_results === true",
                
                # MR Results
                div(
                  style = "display:flex; align-items:center; justify-content:space-between;",
                  tags$h4("MR Results"),
                  actionButton("info_mr_results", label = NULL, icon = icon("info-circle"), class = "btn-link")
                ),
                div(
                  class = "note",
                  "The MR Results table summarizes causal effect estimates obtained using different Mendelian Randomization methods. 
  The null hypothesis (H₀) states that there is no causal effect of the exposure on the outcome (β = 0; OR = 1). 
  The p-value evaluates the probability of observing the estimated effect under H₀. 
  A small p-value (e.g., < 0.05) provides evidence to reject H₀, supporting a potential causal association."
                ),
                DTOutput("tbl_mr_results"),
                
                # Heterogeneity
                div(
                  style = "display:flex; align-items:center; justify-content:space-between;",
                  tags$h4("Heterogeneity statistics"),
                  actionButton("info_het", label = NULL, icon = icon("info-circle"), class = "btn-link")
                ),
                div(
                  class = "note",
                  "Heterogeneity refers to the variability in SNP-specific causal estimates. 
  Cochran’s Q statistic tests the null hypothesis that all genetic instruments 
  estimate a common causal effect. A small Q-test p-value (e.g., < 0.05) 
  indicates statistically significant heterogeneity, which may suggest 
  horizontal pleiotropy, invalid instruments, or model misspecification."
                ),
                DTOutput("tbl_het"),
                
                # Horizontal pleiotropy
                div(
                  style = "display:flex; align-items:center; justify-content:space-between;",
                  tags$h4("Horizontal pleiotropy"),
                  actionButton("info_pleio", label = NULL, icon = icon("info-circle"), class = "btn-link")
                ),
                div(
                  class = "note",
                  "Horizontal pleiotropy occurs when genetic instruments (SNPs) affect the outcome through pathways other than the exposure of interest, thereby violating the exclusion restriction assumption of Mendelian Randomization. 
  If pleiotropic effects are directionally unbalanced, they can bias causal effect estimates. 
  The MR-Egger regression intercept provides an estimate of the average directional pleiotropic effect across instruments. 
  A statistically significant intercept (p < 0.05) suggests the presence of directional horizontal pleiotropy."
                ),
                verbatimTextOutput("pleio_summary"),
                DTOutput("tbl_pleio"),
                
                # Causal direction
                div(
                  style = "display:flex; align-items:center; justify-content:space-between;",
                  tags$h4("Causal direction test"),
                  actionButton("info_dir", label = NULL, icon = icon("info-circle"), class = "btn-link")
                ),
                div(
                  class = "note",
                  "The causal direction (Steiger) test evaluates whether genetic instruments explain more variance in the exposure than in the outcome. 
  The null hypothesis states that the instruments explain equal or greater variance in the outcome, suggesting possible reverse causation. 
  A small p-value supports the specified direction (Exposure → Outcome). 
  For binary outcomes, variance is estimated on the liability scale."
                ),
                DTOutput("tbl_dir")
              )
            ),
            div(class = "sep"),
            # Plots panel
            div(
              class = "panel",
              div(class = "panel-title", "Plots"),
              
              # Message shown when no results
              conditionalPanel(
                condition = "output.tab3_has_results === false",
                tags$p(
                  style = "opacity:0.7; font-style:italic; margin-top:8px;",
                  "Plots will appear here after running the MR analysis."
                )
              ),
              
              # Actual plots shown when results exist
              conditionalPanel(
                condition = "output.tab3_has_results === true",
                
                tags$h4("Method comparison plot"),
                plotOutput("plot_scatter", height = "380px"),
                div(class = "note",
                    "SNP effects on the outcome are plotted against SNP effects on the exposure (all SNPs with negative effects on the exposure are shown to be positive, with the sign of the effect on the outcome flipped). The slope of the line represents the causal association, and each method has a different line. The Egger estimate is the only line which doesn't automatically pass through the origin."
                ),
                
                tags$h4("Single SNP analysis"),
                plotOutput("plot_forest", height = "380px"),
                div(class = "note",
                    "The causal effect of exposure on outcome is estimated using each SNP singly using the Wald ratio, and represented in a forest plot. The MR estimate using all SNPs using the MR Egger and IVW methods are also shown. Formal estimates of heterogeneity are shown in the tables above."
                ),
                
                tags$h4("Leave-one-out analysis"),
                plotOutput("plot_loo", height = "380px"),
                div(class = "note",
                    "Leave-one-out sensitivity analysis is performed to ascertain if an association is being disproportionately influenced by a single SNP. Each black point in the forest plot represents the MR analysis (using IVW) excluding that particular SNP. The overall analysis including all SNPs is also shown for comparison."
                ),
                
                tags$h4("Funnel plot"),
                plotOutput("plot_funnel", height = "380px"),
                div(class = "note",
                    "Funnel plot to assess heterogeneity. Less precise estimates (lower values on the y-axis) 'funnel' in as they increase in precision. Larger spread suggests higher heterogeneity, which may be due to horizontal pleiotropy. Asymmetry in the funnel plot indicates directional horizontal pleiotropy which can bias many MR methods. MR Egger regression should guard against this."
                )
              )
            )
          )
        )
      )
    )
  )
}

# Tab 4: Multivariable MR workflow
create_multivariable_mr_tab <- function() {
  tabPanel(
    "Multivariable MR",
    div(
      class = "card",
      div(
        class = "card-header",
        div(class = "card-title", "Multivariable Mendelian Randomization Workflow"),
        div(class = "card-subtitle", "Step-by-step workflow: Identify repeated SNP studies \u2192 Select Exposure 2 or Mediator \u2192 Run Exposure 1 to Exposure 2 / Mediator MR \u2192 Build SNP table \u2192 Run MVMR.")
      ),
      div(
        class = "card-body",
        tabsetPanel(
          # A) Repeated SNP Independent Studies
          tabPanel(
            "A) Repeated SNP Independent Studies",
            div(
              class = "card",
              div(
                class = "card-header",
                div(class = "card-title", "Repeated SNP studies (from Tab 3 SNP table)"),
                div(
                  class = "card-subtitle",
                  "Select Exposure 2 / Mediator from the table below or proceed to Tab B) Select Exposure 2 / Mediator.The Top Category and Sub Category tables below will filter the Repeated SNP studies table based on the selection."
                ),
              ),
              div(
                class = "card-body",
                div(
                  class = "panel",
                  # Independence filter: show only same‑ancestry repeated‑SNP studies
                  checkboxInput(
                    "rep_only_indep",
                    "Show only same-ancestry repeated-SNP studies (recommended)",
                    TRUE
                  ),
                  # Second checkbox: when the first is checked and cohort info exists,
                  # allow filtering to non‑overlapping cohorts
                  conditionalPanel(
                    condition = "input.rep_only_indep && output.has_cohort_union",
                    checkboxInput(
                      "rep_non_overlap",
                      "Within same ancestry, show only non-overlapping cohorts (to Avoid Sample Overlap)",
                      value = FALSE
                    )
                  )
                ),
                
                div(class = "sep"),
                DTOutput("tbl_repeated"),
                div(class = "sep"),
                fluidRow(
                  column(6, DTOutput("tbl_cat_counts")),
                  column(6, DTOutput("tbl_subcat_counts"))
                ),
                actionButton("clear_rep_filters", "Clear filters in Top Category Section & Sub Category Section", class = "btn-ghost w100")
              )
            )
          ),
          # B) Select Exposure 2
          tabPanel(
            "B) Select Exposure 2 / Mediator",
            div(class = "sep"),
            div(
              class = "card",
              div(
                class = "card-header",
                div(class = "card-title", "Select Exposure 2"),
                div(class = "card-subtitle", "Exposure 1 is from Tab 1. Exposure 2 / Mediator is from repeated-SNP candidates or New Mediator. Next tabs can run Exposure 1 \u2192 Exposure 2 / Mediator MR and build the SNP table for MVMR.")
              ),
              div(
                class = "card-body",
                # Auto-filled Exposure 1 box
                div(
                  class = "panel",
                  div(class = "panel-title", "Auto-filled Exposure 1 (from Tab 1)"),
                  verbatimTextOutput("exp1_auto_box")
                ),
                div(class = "sep"),
                # Mediator source selection
                div(
                  class = "card-body",
                  radioButtons(
                    "mediator_source",
                    "Mediator source",
                    choices = c(
                      "Use Exposure 2 (from Tab A: Repeated SNP Studies)" = "from_tabA",
                      "Search a New Exposure 2 / Mediator" = "new_search"
                    ),
                    selected = "from_tabA"
                  ),
                  conditionalPanel(
                    condition = "input.mediator_source == 'new_search'",
                    
                    div(class = "sep"),
                    searchUI(
                      "tab7_mediator_search",
                      "Search & Select NEW mediator",
                      show_selection_box = TRUE,
                      selection_title = "Selected Mediator GCST",
                      extra_ui =   div(
                        class = "panel",
                        div(class = "panel-title", "Mediator independence filter"),
                        # First checkbox: show only same‑ancestry mediators (shared ancestry with exposure and outcome)
                        checkboxInput(
                          "tab7_only_indep",
                          "Show only same-ancestry mediators (recommended)",
                          value = TRUE
                        ),
                        # Second checkbox: show only non‑overlapping cohorts within same ancestry
                        # Appears only when cohort information exists for the selected exposure/outcome pair
                        conditionalPanel(
                          condition = "input.tab7_only_indep && output.has_cohort_union",
                          checkboxInput(
                            "tab7_non_overlap",
                            "Within same ancestry, show only non-overlapping cohorts (to Avoid Sample Overlap)",
                            value = FALSE
                          )
                        )
                      )
                    )
                  ),
                  div(class = "sep"),
                  # Exposure 2 picker from repeated-SNP studies
                  conditionalPanel(
                    condition = "input.mediator_source == 'from_tabA'",
                    div(
                      class = "panel",
                      div(class = "panel-title", "Pick Exposure 2 from repeated-SNP studies"),
                      uiOutput("exp2_pick_ui"),
                      div(class = "sep"),
                      div(class = "label", "Currently selected Exposure 2"),
                      verbatimTextOutput("exp2_pick_box")
                    )
                  ),
                  div(class = "sep"),
                  
                  
                )
              )
            )
          ),
          # C) Run Exposure 1 → Mediator MR
          tabPanel(
            "C) Run Exposure 1 \u2192 Exposure 2 / Mediator MR",
            div(
              class = "card",
              div(
                class = "card-header",
                div(class = "card-title", "Two-sample MR: Exposure 1 \u2192 Exposure 2 / Mediator"),
                div(class = "card-subtitle", "Run MR for the selected Exposure 1 against the Exposure 2 / Mediator")
              ),
              div(
                class = "card-body",
                fluidRow(
                  column(
                    4,
                    # Tab5 instrument settings
                    div(
                      class = "panel",
                      div(class = "panel-title", "Instrument settings"),
                      checkboxInput("tab5_do_clump", "Use LD clumping (OpenGWAS)", value = TRUE),
                      numericInput("tab5_p_thresh", "Instrument p-value (p1)", value = 5e-8, min = 0, max = 1, step = 1e-8),
                      sliderInput("tab5_clump_r2", "Clump r\u00b2", min = 0.0001, max = 0.2, value = 0.01, step = 0.0005),
                      numericInput("tab5_clump_kb", "Clump window (kb)", value = 10000, min = 250, max = 50000, step = 250)
                    ),
                    div(class = "sep"),
                    # Tab5 outcome extraction options
                    div(
                      class = "panel",
                      div(class = "panel-title", "Outcome extraction options"),
                      checkboxInput("tab5_use_proxies", "Use proxies if SNP missing", value = TRUE),
                      sliderInput("tab5_proxy_rsq", "Min proxy r\u00b2", min = 0.1, max = 1.0, value = 0.8, step = 0.05),
                      checkboxInput("tab5_allow_palindromes", "Allow palindromic SNPs", value = TRUE),
                      sliderInput("tab5_maf_thresh", "MAF threshold", min = 0.0, max = 0.5, value = 0.30, step = 0.01),
                      radioButtons(
                        "tab5_harm_action",
                        "Harmonise action",
                        choices = c(
                          "Positive strand" = "pos",
                          "Align (recommended)" = "align",
                          "Exclude palindromes" = "exclude"
                        ),
                        selected = "align"
                      )
                    ),
                    div(class = "sep"),
                    # Tab5 MR methods & downloads
                    div(
                      class = "panel",
                      div(class = "panel-title", "MR methods"),
                      checkboxGroupInput(
                        "tab5_mr_methods",
                        "Methods",
                        choices = mr_method_choices(),
                        selected = c(
                          "mr_wald_ratio",
                          "mr_ivw",
                          "mr_egger_regression",
                          "mr_weighted_median",
                          "mr_weighted_mode"
                        )
                      ),
                      actionButton("run_mr_tab5", "\u26a1 Run MR (Exposure1\u2192Mediator)", class = "btn-primary w100"),
                      div(class = "sep"),
                      div(class = "panel-title", "Downloads (current selected pair)"),
                      
                      conditionalPanel(
                        condition = "output.tab5_has_results === true",
                        
                        downloadButton(
                          "dl_tab5_mr_xlsx",
                          "⬇ MR Results (Excel)",
                          class = "btn-ghost w100"
                        ),
                        div(class="sep"),
                        
                        downloadButton(
                          "dl_tab5_sens_xlsx",
                          "⬇ Sensitivity Tables (Excel)",
                          class = "btn-ghost w100"
                        ),
                        div(class="sep"),
                        
                        downloadButton(
                          "dl_tab5_plot_scatter",
                          "⬇ Method comparison plot (PDF)",
                          class = "btn-ghost w100"
                        ),
                        div(class="sep"),
                        
                        downloadButton(
                          "dl_tab5_plot_forest",
                          "⬇ Single SNP analysis plot (PDF)",
                          class = "btn-ghost w100"
                        ),
                        div(class="sep"),
                        
                        downloadButton(
                          "dl_tab5_plot_loo",
                          "⬇ Leave-one-out analysis plot (PDF)",
                          class = "btn-ghost w100"
                        ),
                        div(class="sep"),
                        
                        downloadButton(
                          "dl_tab5_plot_funnel",
                          "⬇ Funnel plot (PDF)",
                          class = "btn-ghost w100"
                        ),
                        div(class="sep"),
                        
                        downloadButton(
                          "dl_tab5_snp_table_xlsx",
                          "⬇ SNP Effects Table (Excel)",
                          class = "btn-ghost w100"
                        )
                      ),
                      
                      conditionalPanel(
                        condition = "output.tab5_has_results === false",
                        
                        tags$p(
                          style="opacity:0.6; font-style:italic;",
                          "Downloads will appear after running the MR analysis."
                        )
                      )
                    )
                  ),
                  column(
                    8,
                    # Tab5 status and pair picker
                    div(class = "sep"),
                    div(
                      class = "panel",
                      div(class = "panel-title", "Status"),
                      verbatimTextOutput("tab5_mr_status"),
                      tags$hr(),
                      verbatimTextOutput("opengwas_countdown")
                    ),
                    div(
                      class = "panel",
                      div(class = "panel-title", "Choose pair (X → M)"),
                      
                      uiOutput("tab5_pair_ui"),
                      
                      conditionalPanel(
                        condition = "output.tab5_has_results === false",
                        tags$p(
                          style = "margin-top:10px; opacity:0.7; font-style:italic;",
                          "No MR results yet."
                        )
                      )
                    ),
                    div(class = "sep"),
                    # Tab5 results & statistics
                    div(
                      class = "panel",
                      div(class = "panel-title", "MR Results & Statistics"),
                      
                      # Message BEFORE MR run
                      conditionalPanel(
                        condition = "output.tab5_has_results === false",
                        tags$p(
                          style = "opacity:0.7; font-style:italic;",
                          "MR Results & Statistics will appear here after running the MR analysis."
                        )
                      ),
                      
                      # Actual tables AFTER MR run
                      conditionalPanel(
                        condition = "output.tab5_has_results === true",
                        
                        # ===============================
                        # MR Results
                        # ===============================
                        div(
                          style = "display:flex; align-items:center; justify-content:space-between;",
                          tags$h4("MR Results"),
                          actionButton("info_mr_results",
                                       label = NULL,
                                       icon = icon("info-circle"),
                                       class = "btn-link")
                        ),
                        div(
                          class = "note",
                          "The MR Results table summarizes causal effect estimates obtained using different Mendelian Randomization methods. 
  The null hypothesis (H₀) states that there is no causal effect of the exposure on the outcome (β = 0; OR = 1). 
  The p-value evaluates the probability of observing the estimated effect under H₀. 
  A small p-value (e.g., < 0.05) provides evidence to reject H₀, supporting a potential causal association."
                        ),
                        DTOutput("tbl_tab5_mr_results"),
                        div(class = "sep"),
                        
                        # ===============================
                        # Heterogeneity
                        # ===============================
                        div(
                          style = "display:flex; align-items:center; justify-content:space-between;",
                          tags$h4("Heterogeneity statistics"),
                          actionButton("input$info_het",
                                       label = NULL,
                                       icon = icon("info-circle"),
                                       class = "btn-link")
                        ),
                        div(
                          class = "note",
                          "Heterogeneity refers to the variability in SNP-specific causal estimates. Cochran’s Q statistic tests the null hypothesis that all genetic instruments estimate a common causal effect. A small Q-test p-value (e.g., < 0.05) indicates statistically significant heterogeneity, which may suggest horizontal pleiotropy, invalid instruments, or model misspecification."
                        ),
                        DTOutput("tbl_tab5_het"),
                        div(class = "sep"),
                        
                        # ===============================
                        # Horizontal pleiotropy
                        # ===============================
                        div(
                          style = "display:flex; align-items:center; justify-content:space-between;",
                          tags$h4("Horizontal pleiotropy"),
                          actionButton("input$info_pleio",
                                       label = NULL,
                                       icon = icon("info-circle"),
                                       class = "btn-link")
                        ),
                        div(
                          class = "note",
                          p("Horizontal pleiotropy occurs when genetic instruments (SNPs) affect the outcome through pathways other than the exposure of interest, thereby violating the exclusion restriction assumption of Mendelian Randomization. If pleiotropic effects are directionally unbalanced, they can bias causal effect estimates. The MR-Egger regression intercept provides an estimate of the average directional pleiotropic effect across instruments. A statistically significant intercept (p < 0.05) suggests the presence of directional horizontal pleiotropy.")
                        ),
                        verbatimTextOutput("input$info_dir"),
                        DTOutput("tbl_tab5_pleio"),
                        div(class = "sep"),
                        
                        # ===============================
                        # Causal direction
                        # ===============================
                        div(
                          style = "display:flex; align-items:center; justify-content:space-between;",
                          tags$h4("Causal direction test"),
                          actionButton("info_tab5_dir",
                                       label = NULL,
                                       icon = icon("info-circle"),
                                       class = "btn-link")
                        ),
                        div(
                          class = "note",
                          "The causal direction (Steiger) test evaluates whether genetic instruments explain more variance in the exposure than in the outcome. The null hypothesis states that the instruments explain equal or greater variance in the outcome, suggesting possible reverse causation. A small p-value supports the specified direction (Exposure → Outcome). For binary outcomes, variance is estimated on the liability scale."
                        ),
                        DTOutput("tbl_tab5_dir")
                      )
                    ),
                    div(class = "sep"),
                    # Tab5 plots
                    # Tab5 plots
                    div(
                      class = "panel",
                      div(class = "panel-title", "Plots"),
                      
                      conditionalPanel(
                        condition = "output.tab5_has_results === true",
                        
                        # Method comparison
                        tags$h4("Method comparison plot"),
                        plotOutput("tab5_plot_scatter", height = "360px"),
                        div(class = "note",
                            "SNP effects on the outcome are plotted against SNP effects on the exposure (all SNPs with negative effects on the exposure are shown to be positive, with the sign of the effect on the outcome flipped). The slope of the line represents the causal association, and each method has a different line. The Egger estimate is the only line which does not automatically pass through the origin."
                        ),
                        
                        # Single SNP
                        tags$h4("Single SNP analysis"),
                        plotOutput("tab5_plot_forest", height = "360px"),
                        div(class = "note",
                            "The causal effect of exposure on outcome is estimated using each SNP singly using the Wald ratio, and represented in a forest plot. The MR estimate using all SNPs with MR Egger and IVW methods are also shown. Formal estimates of heterogeneity are reported in the tables above."
                        ),
                        
                        # Leave-one-out
                        tags$h4("Leave-one-out analysis"),
                        plotOutput("tab5_plot_loo", height = "360px"),
                        div(class = "note",
                            "Leave-one-out sensitivity analysis is performed to assess whether the association is disproportionately influenced by a single SNP. Each point represents the MR analysis excluding one SNP at a time. The overall estimate using all SNPs is shown for comparison."
                        ),
                        
                        # Funnel plot
                        tags$h4("Funnel plot"),
                        plotOutput("tab5_plot_funnel", height = "360px"),
                        div(class = "note",
                            "Funnel plot to assess heterogeneity. Less precise estimates (lower precision) funnel inward as precision increases. Larger spread suggests higher heterogeneity, which may reflect horizontal pleiotropy. Asymmetry may indicate directional pleiotropy. MR Egger regression can help account for this."
                        )
                      ),
                      
                      conditionalPanel(
                        condition = "output.tab5_has_results === false",
                        tags$p(
                          style = "opacity:0.7; font-style:italic;",
                          "Plots will appear here after running the MR analysis."
                        )
                      )
                    )
                  )
                )
              )
            )
          ),
          # D) Build SNP Table
          tabPanel(
            "D) Build SNP Table",
            div(
              class = "card",
              div(
                class = "card-header",
                div(class = "card-title", "Build SNP Table for MVMR"),
                div(class = "card-subtitle", "Build the combined SNP effects table used for the MVMR model.")
              ),
              div(
                class = "card-body",
                div(
                  class = "panel",
                  div(class = "panel-title", ""),
                  div(class = "note", ""),
                  div(class = "sep"),
                  actionButton("build_tab5", "\ud83e\udde9 Build SNP Table", class = "btn-primary w100"),
                  div(class="sep"),
                  
                  conditionalPanel(
                    condition = "output.tab5_snp_ready === true",
                    downloadButton(
                      "dl_tab5_snp_mvmr_xlsx",
                      "⬇ Download SNP Table (Excel)",
                      class = "btn-ghost w100"
                    )
                  ),
                  
                  conditionalPanel(
                    condition = "output.tab5_snp_ready === false",
                    tags$p(
                      style="opacity:0.6; font-style:italic;",
                      "Download available after building SNP table."
                    )
                  ),
                  div(class = "sep"),
                  verbatimTextOutput("tab5_warn"),
                  DTOutput("tbl_tab5_snp_mvmr")
                )
              )
              
            )
          ),
          # E) Run MVMR
          tabPanel(
            "E) Run MVMR",
            div(
              class = "card",
              div(
                class = "card-header",
                div(class = "card-title", "Multivariable MR"),
                div(class = "card-subtitle", "MVMR estimates the direct causal effects of Exposure 1 and Exposure 2 / Mediator on the Outcome using summary-level GWAS data.")
              ),
              div(
                class = "card-body",
                fluidRow(
                  column(
                    4,
                    div(
                      class = "panel",
                      div(class = "panel-title", ""),
                      div(class = "note", ""),
                      div(class = "sep"),
                      actionButton("run_mvmr", "\u26a1 Run MVMR", class = "btn-primary w100")
                    ),
                    div(class = "sep"),
                    div(
                      class = "panel",
                      div(class = "panel-title", "Downloads"),
                      
                      conditionalPanel(
                        condition = "output.mvmr_has_results === true",
                        downloadButton("dl_mvmr_full_xlsx", "⬇ Download Full MVMR Results", class = "btn-ghost w100")
                      ),
                      
                      conditionalPanel(
                        condition = "output.mvmr_has_results === false",
                        tags$p(
                          style = "opacity:0.7; font-style:italic;",
                          "Download will appear after MVMR results are available."
                        )
                      )
                    )
                  ),
                  column(
                    8,
                    div(class = "panel", div(class = "panel-title", "Status"), verbatimTextOutput("mvmr_status")),
                    div(class = "sep"),
                    ## New MVMR Results panel (without tabs)
                    div(
                      class = "panel",
                      div(class = "panel-title", "MVMR Results & Statistics"),
                      conditionalPanel(
                        condition = "output.mvmr_has_results === false",
                        tags$p(
                          style = "opacity:0.7; font-style:italic;",
                          "MR Results & Statistics will appear here after running the MVMR analysis."
                        )
                      ),
                      conditionalPanel(
                        condition = "output.mvmr_has_results === true",
                        # MVMR IVW result section
                        div(
                          tags$h4("MVMR IVW result")
                        ),
                        div(
                          class = "note",
                          "MVMR IVW results represent direct causal effects of each exposure on the outcome, adjusted for other exposures. 
The null hypothesis states that the direct effect equals zero (β = 0). 
A small p-value (< 0.05) suggests statistical evidence of a direct causal effect."
                        ),
                        DTOutput("tbl_mvmr_ivw"),
                        div(class = "sep"),
                        # Instrument strength section
                        div(
                          tags$h4("Instrument strength")
                        ),
                        div(
                          class = "note",
                          "Instrument strength assessed using F-statistics. Conventional rule: F > 10 suggests low weak-instrument bias.If your conditional F-statistic < 10 in MVMR, it means your instruments are weak for that exposure (conditional on the other exposure), so estimates can be biased and SE/p-values can be unreliable."
                        ),
                        DTOutput("tbl_mvmr_strength"),
                        div(class = "sep"),
                        # Pleiotropy section
                        div(
                          tags$h4("Test for horizontal pleiotropy")
                        ),
                        div(
                          class = "note",
                          "Evaluated using heterogeneity (Q-statistic).Heterogeneity reflects the variability in SNP-specific causal estimates. Cochran’s Q statistic tests the null hypothesis that all genetic instruments estimate a common causal effect. A significant Q-test (p < 0.05) suggests excess heterogeneity, which may indicate horizontal pleiotropy or violations of MR assumptions."
                        ),
                        DTOutput("tbl_mvmr_pleio")
                      )
                    )
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

# Tab 5: Mediation Analysis
create_mediation_tab <- function() {
  tabPanel(
    "Mediation Analysis",
    div(
      class = "card",
      div(
        class = "card-header",
        div(class = "card-title", "Mediation Analysis"),
        div(class = "card-subtitle", "Mediation analysis is used to understand how or through which pathway an exposure affects an outcome. Instead of asking only whether an exposure causes an outcome, mediation analysis asks whether part of that effect operates through an intermediate variable (the mediator).")
      ),
      div(
        class = "card-body",
        fluidRow(
          column(
            12,
            # The contents of this tab are generated dynamically in the server based on the mediator source.
            uiOutput("ui_tab7_main_tabs")
          )
        )
      )
    )
  )
}

##----------------------------------------------------------------------------
##  Assemble the UI
##----------------------------------------------------------------------------

# Compose the full UI by combining the styles, scripts, header card and tab
# panels.  Encapsulating this in a function allows for lazy evaluation and
# clear separation of concerns.
build_ui <- function() {
  fluidPage(
    tags$head(
      tags$style(HTML(custom_css)),
      custom_scroll_script(),
      tags$script(HTML("
    Shiny.addCustomMessageHandler('playSound', function(message) {
      var audio = new Audio('https://actions.google.com/sounds/v1/alarms/beep_short.ogg');
      audio.play();
    });
  "))
    ),
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

# Instantiate the UI.  Calling the function here ensures that everything is
# defined in the global environment before the server is sourced or defined.
ui <- build_ui()

##----------------------------------------------------------------------------
##  Run the Application
##----------------------------------------------------------------------------

# Instantiate the Shiny app.  
shinyApp(ui, server)