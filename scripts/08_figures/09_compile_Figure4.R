
# ============================================================
# Compile final Figure 4
# A: Total 16S-normalized ARG abundance
# B: ARG-group richness
# C: ARG-class persistence
# D: Recurrently detected PMA ARG groups
# ============================================================

.libPaths(c(path.expand("~/R/4.5"), .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

OUTDIR <- "figures/figure4"

# ------------------------------------------------------------
# Read individual PDF panels
# using patchwork + magick if available
# ------------------------------------------------------------

if (!requireNamespace("magick", quietly = TRUE)) {
  stop("R package 'magick' is required for PDF panel compilation.")
}

read_panel <- function(path) {
  img <- magick::image_read_pdf(
    path,
    density = 300
  )

  patchwork::wrap_elements(
    full = grid::rasterGrob(
      as.raster(img),
      interpolate = TRUE
    )
  )
}

A <- read_panel(
  file.path(OUTDIR, "Figure4A_ARG_burden.pdf")
)

B <- read_panel(
  file.path(OUTDIR, "Figure4B_ARG_richness.pdf")
)

C <- read_panel(
  file.path(OUTDIR, "Figure4C_ARG_class_persistence.pdf")
)

D <- read_panel(
  file.path(OUTDIR, "Figure4D_persistent_ARGs.pdf")
)

# ------------------------------------------------------------
# Layout
#
# Top row slightly shorter.
# Bottom row gets more vertical space for labels.
# ------------------------------------------------------------

final_fig <-
  (A | B) /
  (C | D) +
  plot_layout(
    widths = c(1, 1),
    heights = c(0.88, 1.12)
  ) &
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

ggsave(
  file.path(
    OUTDIR,
    "Figure4_FINAL.pdf"
  ),
  final_fig,
  width = 11.2,
  height = 9.0,
  units = "in",
  bg = "white"
)

ggsave(
  file.path(
    OUTDIR,
    "Figure4_FINAL.png"
  ),
  final_fig,
  width = 11.2,
  height = 9.0,
  units = "in",
  dpi = 400,
  bg = "white"
)

cat("\n===== FIGURE 4 FINAL =====\n")
cat("Layout: A | B / C | D\n")
cat("Colors: non-PMA #1479E8; PMA #FF4B4B\n")
cat("Saved:\n")
cat("  figures/figure4/Figure4_FINAL.pdf\n")
cat("  figures/figure4/Figure4_FINAL.png\n")

