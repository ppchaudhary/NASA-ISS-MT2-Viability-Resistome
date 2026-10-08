library(ggplot2)
library(patchwork)

outdir <- "figures/figure6"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Helper: source panel script in isolated environment
# and retrieve ggplot object "p"
# ============================================================

get_panel <- function(script) {

  e <- new.env(parent = globalenv())

  cat("Loading:", script, "\n")

  sys.source(
    script,
    envir = e
  )

  if (!exists("p", envir = e)) {
    stop(
      paste(
        "No ggplot object named 'p' found in",
        script
      )
    )
  }

  get("p", envir = e)
}

# ============================================================
# Load TRUE VECTOR ggplot objects
# ============================================================

cat("\nLoading original Figure 6 ggplot objects...\n")

pA <- get_panel(
  "scripts/08_figures/15_Figure6A_ARG_taxon_composition.R"
)

pB <- get_panel(
  "scripts/08_figures/16_Figure6B_recurrent_ARG_taxon.R"
)

pC <- get_panel(
  "scripts/08_figures/17_Figure6C_ARG_taxon_flow.R"
)

# ============================================================
# Small panel-specific adjustments for compiled figure
# ============================================================

pA <- pA +
  theme(
    plot.margin = margin(8, 8, 8, 8)
  )

pB <- pB +
  theme(
    plot.margin = margin(8, 8, 8, 8)
  )

pC <- pC +
  theme(
    plot.margin = margin(8, 15, 8, 15)
  )

# ============================================================
# Layout
#
# A = large left panel
# B = upper right
# C = lower right
# ============================================================

right_side <- pB / pC +
  plot_layout(
    heights = c(1.05, 0.95)
  )

fig6 <- pA | right_side

fig6 <- fig6 +
  plot_layout(
    widths = c(1.05, 1.35)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 18,
        color = "black"
      )
    )
  )

# ============================================================
# Save TRUE VECTOR PDF
# ============================================================

ggsave(
  file.path(
    outdir,
    "Figure6_FINAL_VECTOR.pdf"
  ),
  fig6,
  width = 16,
  height = 11.5,
  device = cairo_pdf
)

# High-resolution preview
ggsave(
  file.path(
    outdir,
    "Figure6_FINAL_VECTOR.png"
  ),
  fig6,
  width = 16,
  height = 11.5,
  dpi = 500
)

cat("\n========================================\n")
cat("FIGURE 6 COMPILED\n")
cat("========================================\n")
cat(
  "PDF: figures/figure6/Figure6_FINAL_VECTOR.pdf\n"
)
cat(
  "PNG: figures/figure6/Figure6_FINAL_VECTOR.png\n"
)
cat("TRUE VECTOR PDF: yes\n")
cat("Layout: A left | B/C right\n")
cat("========================================\n")
