library(ggplot2)
library(patchwork)

outdir <- "figures/figure7"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Source panel scripts and retrieve ggplot object "p"
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

cat("\nLoading Figure 7 panels...\n")

pA <- get_panel(
  "scripts/08_figures/19_Figure7A_microbiome_resistome.R"
)

pB <- get_panel(
  "scripts/08_figures/20_Figure7B_microbiome_ARGburden.R"
)

pC <- get_panel(
  "scripts/08_figures/21_Figure7C_resistome_ARGburden.R"
)

# ============================================================
# Consistent compiled margins
# ============================================================

pA <- pA +
  theme(
    plot.margin = margin(8, 10, 8, 8)
  )

pB <- pB +
  theme(
    plot.margin = margin(8, 10, 8, 8)
  )

pC <- pC +
  theme(
    plot.margin = margin(8, 10, 8, 8)
  )

# ============================================================
# Horizontal A | B | C layout
# ============================================================

fig7 <- pA | pB | pC

fig7 <- fig7 +
  plot_layout(
    widths = c(1, 1, 1)
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
# Save true-vector PDF + high-resolution PNG
# ============================================================

ggsave(
  file.path(
    outdir,
    "Figure7_FINAL_VECTOR.pdf"
  ),
  fig7,
  width = 17,
  height = 5.8,
  device = cairo_pdf
)

ggsave(
  file.path(
    outdir,
    "Figure7_FINAL_VECTOR.png"
  ),
  fig7,
  width = 17,
  height = 5.8,
  dpi = 500
)

cat("\n========================================\n")
cat("FIGURE 7 COMPILED\n")
cat("========================================\n")
cat("A: Microbiome shift vs resistome shift\n")
cat("B: Microbiome shift vs signed ARG burden change\n")
cat("C: Resistome shift vs absolute ARG burden change\n")
cat("\nPDF: figures/figure7/Figure7_FINAL_VECTOR.pdf\n")
cat("PNG: figures/figure7/Figure7_FINAL_VECTOR.png\n")
cat("TRUE VECTOR PDF: yes\n")
cat("Layout: A | B | C\n")
cat("========================================\n")
