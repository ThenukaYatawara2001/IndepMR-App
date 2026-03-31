# =============================================================================
# R/07_pdf_helpers.R
# PDF report helpers (no LaTeX) using grid + gridExtra
# =============================================================================

suppressPackageStartupMessages({
  library(grid)
  library(gridExtra)
})

# -------------------------------------------------------------------------
# Fallback: if %||% not loaded yet (from your 00 utils),
# define a minimal version here so PDF doesn't crash.
# -------------------------------------------------------------------------
if (!exists("%||%")) {
  `%||%` <- function(a, b) {
    if (!is.null(a) && length(a) > 0 && !is.na(a)) a else b
  }
}

# -------------------------------------------------------------------------
# Fallback: if round_df_display() not loaded yet (from your 06 helpers),
# define a minimal version here so PDF doesn't crash.
# -------------------------------------------------------------------------
if (!exists("round_df_display")) {
  round_df_display <- function(df, digits = 4) {
    if (is.null(df) || nrow(df) == 0) return(df)
    out <- df
    num_cols <- names(out)[vapply(out, is.numeric, logical(1))]
    for (cc in num_cols) out[[cc]] <- round(as.numeric(out[[cc]]), digits)
    out
  }
}

# -------------------------------------------------------------------------
# Header band used on each PDF page
# -------------------------------------------------------------------------
pdf_header <- function(title, subtitle = NULL, right_note = NULL) {
  grid.rect(
    x = 0.5, y = 0.97, width = 1, height = 0.10,
    gp = gpar(fill = "#203a78", col = NA)
  )
  grid.text(
    title, x = 0.02, y = 0.985, just = c("left", "top"),
    gp = gpar(fontsize = 16, fontface = "bold", col = "white")
  )
  if (!is.null(subtitle)) {
    grid.text(
      subtitle, x = 0.02, y = 0.955, just = c("left", "top"),
      gp = gpar(fontsize = 9.5, col = "white")
    )
  }
  if (!is.null(right_note)) {
    grid.text(
      right_note, x = 0.98, y = 0.985, just = c("right", "top"),
      gp = gpar(fontsize = 9, col = "white")
    )
  }
}

# -------------------------------------------------------------------------
# Draw a table chunk. If too many columns, split into 2 stacked tables.
# (Top = first 6 columns, Bottom = remaining columns)
# -------------------------------------------------------------------------
draw_table_chunk_split_cols <- function(df_chunk, theme_obj) {
  if (is.null(df_chunk) || nrow(df_chunk) == 0) {
    grid.text("No data.", x = 0.02, y = 0.86, just = c("left", "top"),
              gp = gpar(fontsize = 11, col = "grey20"))
    return(invisible())
  }
  
  # safer for tableGrob
  df_chunk <- as.data.frame(df_chunk, stringsAsFactors = FALSE)
  
  if (ncol(df_chunk) <= 6) {
    tg <- tableGrob(df_chunk, rows = NULL, theme = theme_obj)
    pushViewport(viewport(x = 0.02, y = 0.88, width = 0.96, height = 0.78, just = c("left", "top")))
    grid.draw(tg)
    popViewport()
    return(invisible())
  }
  
  top_df <- df_chunk[, 1:6, drop = FALSE]
  bot_df <- df_chunk[, 7:ncol(df_chunk), drop = FALSE]
  
  tg1 <- tableGrob(top_df, rows = NULL, theme = theme_obj)
  tg2 <- tableGrob(bot_df, rows = NULL, theme = theme_obj)
  
  pushViewport(viewport(x = 0.02, y = 0.88, width = 0.96, height = 0.42, just = c("left", "top")))
  grid.draw(tg1)
  popViewport()
  
  pushViewport(viewport(x = 0.02, y = 0.43, width = 0.96, height = 0.42, just = c("left", "top")))
  grid.draw(tg2)
  popViewport()
  
  invisible()
}

# -------------------------------------------------------------------------
# Make PDF pages for a table, splitting into multiple pages by rows.
# Automatically rounds numeric columns for nicer display.
# -------------------------------------------------------------------------
pdf_table_pages <- function(title, df, subtitle = NULL, footer = NULL, rows_per_page = 25) {
  grid.newpage()
  pdf_header(title, subtitle = subtitle, right_note = format(Sys.time(), "%Y-%m-%d %H:%M"))
  
  if (is.null(df) || nrow(df) == 0) {
    grid.text("No data.", x = 0.02, y = 0.86, just = c("left", "top"),
              gp = gpar(fontsize = 11, col = "grey20"))
    return(invisible())
  }
  
  df2 <- round_df_display(as.data.frame(df, stringsAsFactors = FALSE), digits = 4)
  n <- nrow(df2)
  starts <- seq(1, n, by = rows_per_page)
  
  thm <- ttheme_minimal(
    base_size = 9,
    core = list(bg_params = list(fill = c("#ffffff", "#f6f8ff"), col = NA)),
    colhead = list(
      fg_params = list(col = "white", fontface = "bold"),
      bg_params = list(fill = "#2f55b9", col = NA)
    )
  )
  
  for (i in seq_along(starts)) {
    st <- starts[i]
    en <- min(st + rows_per_page - 1, n)
    chunk <- df2[st:en, , drop = FALSE]
    
    if (i > 1) {
      grid.newpage()
      pdf_header(
        title,
        subtitle = paste0(
          subtitle %||% "",
          if (!is.null(subtitle)) " — " else "",
          "Page ", i, " (rows ", st, "-", en, " of ", n, ")"
        ),
        right_note = format(Sys.time(), "%Y-%m-%d %H:%M")
      )
    }
    
    draw_table_chunk_split_cols(chunk, thm)
    
    if (!is.null(footer)) {
      grid.text(
        footer, x = 0.02, y = 0.06, just = c("left", "bottom"),
        gp = gpar(fontsize = 9, col = "grey30")
      )
    }
  }
  
  invisible()
}