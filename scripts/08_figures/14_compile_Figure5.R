
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# ============================================================
# Source each panel script and retrieve its ggplot object "p"
# ============================================================

get_panel <- function(script) {

  e <- new.env(parent = globalenv())

  sys.source(
    script,
    envir = e
  )

  if (!exists("p", envir = e)) {
    stop(
      "Plot object 'p' not found in: ",
      script
    )
  }

  get("p", envir = e)
}

cat("Loading original ggplot objects...\n")

A <- get_panel(
  "scripts/08_figures/10_Figure5A_ARG_concordance.R"
)

B <- get_panel(
  "scripts/08_figures/11_Figure5B_AMRFinder_classes.R"
)

C <- get_panel(
  "scripts/08_figures/12_Figure5C_four_method_evidence.R"
)

D <- get_panel(
  "scripts/08_figures/13_Figure5D_class_concordance.R"
)

# ============================================================
# Compile directly from ggplot objects
# NO PDF import
# NO rasterization
# ============================================================

top <- A + B +
  plot_layout(
    widths = c(0.95, 1.25)
  )

bottom <- C + D +
  plot_layout(
    widths = c(1.08, 1.00)
  )

final <- (top / bottom) +
  plot_layout(
    heights = c(0.95, 1.10)
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 18
    )
  )

# ============================================================
# TRUE VECTOR PDF
# ============================================================

ggsave(
  "figures/figure5/Figure5_FINAL_VECTOR.pdf",
  final,
  width = 16,
  height = 12,
  units = "in",
  device = cairo_pdf
)

# High-resolution preview
ggsave(
  "figures/figure5/Figure5_FINAL_VECTOR.png",
  final,
  width = 16,
  height = 12,
  units = "in",
  dpi = 500,
  bg = "white"
)

cat("\n====================================\n")
cat("TRUE-VECTOR FIGURE 5 COMPLETE\n")
cat("====================================\n")
cat("No panel PDFs were rasterized.\n")
cat("No draw_image() used.\n\n")
cat("figures/figure5/Figure5_FINAL_VECTOR.pdf\n")
cat("figures/figure5/Figure5_FINAL_VECTOR.png\n")

