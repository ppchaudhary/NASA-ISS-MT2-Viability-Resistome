
# ============================================================
# FINAL Figure 2 assembly
#
# A = Observed richness
# B = Shannon diversity
# C = Bray-Curtis PCoA
#
# Statistics:
# A: paired Wilcoxon P = 8.34e-7
# B: paired Wilcoxon P = 0.194
# C: pair-restricted PERMANOVA
#    R2 = 0.0289, P = 0.0001
#
# PERMDISP P = 0.141 reported in caption, not panel.
# ============================================================

.libPaths(c(path.expand("~/R/4.5"), .libPaths()))

suppressPackageStartupMessages({
    library(png)
    library(grid)
})

indir  <- "figures/figure2"
outdir <- "figures/figure2"

# ------------------------------------------------------------
# Read existing high-resolution panels
# ------------------------------------------------------------

A <- png::readPNG(
    file.path(indir, "Fig2B_Observed_paired.png")
)

B <- png::readPNG(
    file.path(indir, "Fig2A_Shannon_paired.png")
)

C <- png::readPNG(
    file.path(indir, "Fig2B_Bray_PCoA.png")
)

# ------------------------------------------------------------
# Draw final figure
# ------------------------------------------------------------

draw_final <- function() {

    grid.newpage()

    pushViewport(
        viewport(
            layout = grid.layout(
                nrow = 1,
                ncol = 3,
                widths = unit(
                    c(0.32, 0.32, 0.36),
                    "npc"
                )
            )
        )
    )

    # ========================================================
    # A — Observed richness
    # ========================================================

    pushViewport(
        viewport(
            layout.pos.row = 1,
            layout.pos.col = 1
        )
    )

    grid.raster(
        A,
        width = unit(0.98, "npc"),
        height = unit(0.98, "npc")
    )

    grid.text(
        expression(
            italic(P) == 8.34 %*% 10^{-7}
        ),
        x = unit(0.50, "npc"),
        y = unit(0.94, "npc"),
        gp = gpar(
            fontsize = 11
        )
    )

    popViewport()

    # ========================================================
    # B — Shannon diversity
    # ========================================================

    pushViewport(
        viewport(
            layout.pos.row = 1,
            layout.pos.col = 2
        )
    )

    grid.raster(
        B,
        width = unit(0.98, "npc"),
        height = unit(0.98, "npc")
    )

    grid.text(
        "B",
        x = unit(0.015, "npc"),
        y = unit(0.985, "npc"),
        just = c("left", "top"),
        gp = gpar(
            fontsize = 20,
            fontface = "bold"
        )
    )

    grid.text(
        expression(
            italic(P) == 0.194
        ),
        x = unit(0.50, "npc"),
        y = unit(0.94, "npc"),
        gp = gpar(
            fontsize = 11
        )
    )

    popViewport()

    # ========================================================
    # C — Bray-Curtis PCoA
    # ========================================================

    pushViewport(
        viewport(
            layout.pos.row = 1,
            layout.pos.col = 3
        )
    )

    grid.raster(
        C,
        width = unit(0.98, "npc"),
        height = unit(0.98, "npc")
    )

    grid.text(
        "C",
        x = unit(0.015, "npc"),
        y = unit(0.985, "npc"),
        just = c("left", "top"),
        gp = gpar(
            fontsize = 20,
            fontface = "bold"
        )
    )

    grid.text(
        "PERMANOVA: R\u00B2 = 0.0289, P = 0.0001",
        x = unit(0.52, "npc"),
        y = unit(0.94, "npc"),
        gp = gpar(
            fontsize = 10
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
        "Figure2_FINAL.pdf"
    ),
    width = 15,
    height = 5.3,
    useDingbats = FALSE
)

draw_final()

dev.off()

# ------------------------------------------------------------
# High-resolution PNG
# ------------------------------------------------------------

png(
    file.path(
        outdir,
        "Figure2_FINAL.png"
    ),
    width = 15,
    height = 5.3,
    units = "in",
    res = 600
)

draw_final()

dev.off()

cat("\n===== FINAL FIGURE 2 COMPILED =====\n")
cat("A = Observed richness\n")
cat("B = Shannon diversity\n")
cat("C = Bray-Curtis PCoA\n\n")

cat(
    "PDF: figures/figure2/Figure2_FINAL.pdf\n"
)

cat(
    "PNG: figures/figure2/Figure2_FINAL.png\n"
)

