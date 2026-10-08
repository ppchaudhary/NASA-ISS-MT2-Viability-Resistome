

# ============================================================
# FINAL Figure 3 assembly
#
# A = Stable differential genera heatmap
# B = Top-genus composition
# C = ANCOM-BC2 forest plot
#
# Layout:
# Large heatmap on left
# Composition upper-right
# Forest lower-right
# ============================================================

.libPaths(c(path.expand("~/R/4.5"), .libPaths()))

suppressPackageStartupMessages({
  library(png)
  library(grid)
})

indir  <- "figures/figure3"
outdir <- "figures/figure3"

# ------------------------------------------------------------
# Read existing high-resolution panels
# ------------------------------------------------------------

A <- png::readPNG(
  file.path(indir, "Fig3B_StableGeneraHeatmap.png")
)

B <- png::readPNG(
  file.path(indir, "Fig3A_TopGenusComposition.png")
)

C <- png::readPNG(
  file.path(indir, "Fig3C_ANCOMBC2_Forest.png")
)

# ------------------------------------------------------------
# Drawing function
# ------------------------------------------------------------

draw_final <- function() {
  
  grid.newpage()
  
  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = 2,
        ncol = 2,
        
        # Heatmap gets ~64% of width
        widths = unit(
          c(0.64, 0.36),
          "npc"
        ),
        
        heights = unit(
          c(0.47, 0.53),
          "npc"
        )
      )
    )
  )
  
  # ========================================================
  # A — Heatmap
  # spans both rows on LEFT
  # ========================================================
  
  pushViewport(
    viewport(
      layout.pos.row = 1:2,
      layout.pos.col = 1
    )
  )
  
  grid.raster(
    A,
    width = unit(0.98, "npc"),
    height = unit(0.98, "npc")
  )
  
  grid.text(
    "A",
    x = unit(0.01, "npc"),
    y = unit(0.99, "npc"),
    just = c("left", "top"),
    gp = gpar(
      fontsize = 20,
      fontface = "bold"
    )
  )
  
  popViewport()
  
  # ========================================================
  # B — Taxonomic composition
  # upper RIGHT
  # ========================================================
  
  pushViewport(
    viewport(
      layout.pos.row = 1,
      layout.pos.col = 2
    )
  )
  
  grid.raster(
    B,
    width = unit(0.96, "npc"),
    height = unit(0.96, "npc")
  )
  
  grid.text(
    "B",
    x = unit(0.01, "npc"),
    y = unit(0.99, "npc"),
    just = c("left", "top"),
    gp = gpar(
      fontsize = 20,
      fontface = "bold"
    )
  )
  
  popViewport()
  
  # ========================================================
  # C — ANCOM-BC2 forest plot
  # lower RIGHT
  # ========================================================
  
  pushViewport(
    viewport(
      layout.pos.row = 2,
      layout.pos.col = 2
    )
  )
  
  grid.raster(
    C,
    width = unit(0.98, "npc"),
    height = unit(0.98, "npc")
  )
  
  grid.text(
    "C",
    x = unit(0.01, "npc"),
    y = unit(0.99, "npc"),
    just = c("left", "top"),
    gp = gpar(
      fontsize = 20,
      fontface = "bold"
    )
  )
  
  popViewport()
  
  popViewport()
}

# ------------------------------------------------------------
# PDF
# ------------------------------------------------------------

pdf(
  file.path(
    outdir,
    "Figure3_FINAL_reoriented.pdf"
  ),
  width = 15,
  height = 9,
  useDingbats = FALSE
)

draw_final()

dev.off()

# ------------------------------------------------------------
# 600-dpi PNG
# ------------------------------------------------------------

png(
  file.path(
    outdir,
    "Figure3_FINAL_reoriented.png"
  ),
  width = 15,
  height = 9,
  units = "in",
  res = 600
)

draw_final()

dev.off()

cat("\n===== FIGURE 3 REORIENTED =====\n")
cat("A = stable differential genera heatmap\n")
cat("B = top-genus composition\n")
cat("C = ANCOM-BC2 forest plot\n\n")

cat(
  "PDF: figures/figure3/Figure3_FINAL_reoriented.pdf\n"
)

cat(
  "PNG: figures/figure3/Figure3_FINAL_reoriented.png\n"
)

