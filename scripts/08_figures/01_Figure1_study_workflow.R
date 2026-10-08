library(ggplot2)
library(patchwork)
library(grid)

outdir <- "figures/figure1"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# CONSISTENT MANUSCRIPT PALETTE
# ============================================================

BLUE   <- "#1479E8"   # non-PMA / microbiome
RED    <- "#FF4B4B"   # PMA / resistome
TEAL   <- "#10B8B8"   # taxonomy
PURPLE <- "#8B6BD9"   # integration
ORANGE <- "#F39C34"   # assembly / ARG-host
GREEN  <- "#45A85A"   # final biological interpretation

LIGHT_BLUE   <- "#EAF4FF"
LIGHT_RED    <- "#FFF0EF"
LIGHT_TEAL   <- "#E8F8F7"
LIGHT_PURPLE <- "#F1EDFB"
LIGHT_ORANGE <- "#FFF3E4"
LIGHT_GREEN  <- "#EDF8EF"
LIGHT_GREY   <- "#F4F5F6"

EDGE <- "#4D4D4D"

# ============================================================
# HELPERS
# ============================================================

box <- function(xmin, xmax, ymin, ymax,
                fill = "white",
                colour = EDGE,
                linewidth = 0.7,
                radius = NULL) {

  annotate(
    "rect",
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    fill = fill,
    colour = colour,
    linewidth = linewidth
  )
}

txt <- function(x, y, label,
                size = 3.6,
                face = "plain",
                colour = "black",
                lineheight = 0.95) {

  annotate(
    "text",
    x = x,
    y = y,
    label = label,
    size = size,
    fontface = face,
    colour = colour,
    lineheight = lineheight
  )
}

arrow_seg <- function(x, y, xend, yend,
                      colour = EDGE,
                      linewidth = 0.65) {

  annotate(
    "segment",
    x = x,
    y = y,
    xend = xend,
    yend = yend,
    colour = colour,
    linewidth = linewidth,
    arrow = arrow(
      length = unit(0.11, "inches"),
      type = "closed"
    )
  )
}

line_seg <- function(x, y, xend, yend,
                     colour = EDGE,
                     linewidth = 0.65) {

  annotate(
    "segment",
    x = x,
    y = y,
    xend = xend,
    yend = yend,
    colour = colour,
    linewidth = linewidth
  )
}

base_theme <- theme_void(base_size = 12) +
  theme(
    plot.margin = margin(10, 10, 10, 10)
  )

# ============================================================
# PANEL A
# MATCHED STUDY DESIGN
# ============================================================

pA <- ggplot() +

  # ISS source
  box(
    2.0, 8.0, 8.45, 9.55,
    fill = LIGHT_PURPLE,
    colour = PURPLE,
    linewidth = 1
  ) +

  txt(
    5, 9.00,
    "ISS environmental\nmetagenomic samples",
    size = 4.0,
    face = "bold"
  ) +

  arrow_seg(
    5, 8.40,
    5, 7.75,
    colour = PURPLE
  ) +

  # Locations
  box(
    3.0, 7.0, 6.65, 7.70,
    fill = "#FFFFFF",
    colour = PURPLE,
    linewidth = 0.9
  ) +

  txt(
    5, 7.18,
    "32 matched locations",
    size = 4.0,
    face = "bold"
  ) +

  # Split
  line_seg(5, 6.60, 5, 6.15) +
  line_seg(2.7, 6.15, 7.3, 6.15) +

  arrow_seg(
    2.7, 6.15,
    2.7, 5.55,
    colour = BLUE
  ) +

  arrow_seg(
    7.3, 6.15,
    7.3, 5.55,
    colour = RED
  ) +

  # non-PMA
  box(
    0.55, 4.75, 3.45, 5.50,
    fill = LIGHT_BLUE,
    colour = BLUE,
    linewidth = 1.1
  ) +

  txt(
    2.65, 5.02,
    "non-PMA",
    size = 4.7,
    face = "bold",
    colour = BLUE
  ) +

  txt(
    2.65, 4.47,
    "n = 32",
    size = 3.8,
    face = "bold"
  ) +

  txt(
    2.65, 3.92,
    "Total DNA-associated\nmetagenomic signal",
    size = 3.35
  ) +

  # PMA
  box(
    5.25, 9.45, 3.45, 5.50,
    fill = LIGHT_RED,
    colour = RED,
    linewidth = 1.1
  ) +

  txt(
    7.35, 5.02,
    "PMA",
    size = 4.7,
    face = "bold",
    colour = RED
  ) +

  txt(
    7.35, 4.47,
    "n = 32",
    size = 3.8,
    face = "bold"
  ) +

  txt(
    7.35, 3.92,
    "Viability-informed /\nintact-cell-associated signal",
    size = 3.25
  ) +

  # Converge
  arrow_seg(
    2.65, 3.40,
    4.30, 2.60,
    colour = BLUE
  ) +

  arrow_seg(
    7.35, 3.40,
    5.70, 2.60,
    colour = RED
  ) +

  box(
    2.7, 7.3, 1.25, 2.60,
    fill = LIGHT_PURPLE,
    colour = PURPLE,
    linewidth = 1
  ) +

  txt(
    5, 2.15,
    "64 metagenomes",
    size = 4.1,
    face = "bold"
  ) +

  txt(
    5, 1.65,
    "32 matched PMA–non-PMA pairs",
    size = 3.35
  ) +

  coord_cartesian(
    xlim = c(0, 10),
    ylim = c(0.8, 9.9),
    clip = "off"
  ) +

  base_theme


# ============================================================
# PANEL B
# BIOINFORMATICS WORKFLOW
# ============================================================

pB <- ggplot() +

  # Input
  box(
    2.8, 7.2, 9.00, 9.90,
    fill = LIGHT_PURPLE,
    colour = PURPLE,
    linewidth = 1
  ) +

  txt(
    5, 9.45,
    "Paired-end shotgun\nmetagenomic reads",
    size = 3.7,
    face = "bold"
  ) +

  arrow_seg(
    5, 8.95,
    5, 8.45,
    colour = PURPLE
  ) +

  # QC
  box(
    3.0, 7.0, 7.55, 8.40,
    fill = LIGHT_TEAL,
    colour = TEAL,
    linewidth = 0.9
  ) +

  txt(
    5, 7.98,
    "Quality control\nBBDuk",
    size = 3.45,
    face = "bold"
  ) +

  arrow_seg(
    5, 7.50,
    5, 7.00,
    colour = TEAL
  ) +

  # Human removal
  box(
    3.0, 7.0, 6.05, 6.95,
    fill = LIGHT_TEAL,
    colour = TEAL,
    linewidth = 0.9
  ) +

  txt(
    5, 6.50,
    "Human-read removal\nBowtie2",
    size = 3.35,
    face = "bold"
  ) +

  # Branch split
  line_seg(
    5, 6.00,
    5, 5.55
  ) +

  line_seg(
    1.75, 5.55,
    8.25, 5.55
  ) +

  arrow_seg(
    1.75, 5.55,
    1.75, 5.05,
    colour = BLUE
  ) +

  arrow_seg(
    5.00, 5.55,
    5.00, 5.05,
    colour = RED
  ) +

  arrow_seg(
    8.25, 5.55,
    8.25, 5.05,
    colour = ORANGE
  ) +

  # Taxonomy
  box(
    0.25, 3.25, 3.60, 5.00,
    fill = LIGHT_BLUE,
    colour = BLUE,
    linewidth = 1
  ) +

  txt(
    1.75, 4.62,
    "Taxonomic profiling",
    size = 3.55,
    face = "bold",
    colour = BLUE
  ) +

  txt(
    1.75, 4.05,
    "Kraken2\n+ Bracken",
    size = 3.35
  ) +

  # Read-based resistome
  box(
    3.50, 6.50, 3.60, 5.00,
    fill = LIGHT_RED,
    colour = RED,
    linewidth = 1
  ) +

  txt(
    5, 4.62,
    "Read-based resistome",
    size = 3.45,
    face = "bold",
    colour = RED
  ) +

  txt(
    5, 4.08,
    "DeepARG\nShortBRED • RGI",
    size = 3.05
  ) +

  # Assembly
  box(
    6.75, 9.75, 3.60, 5.00,
    fill = LIGHT_ORANGE,
    colour = ORANGE,
    linewidth = 1
  ) +

  txt(
    8.25, 4.62,
    "Assembly-based",
    size = 3.55,
    face = "bold",
    colour = "#D97706"
  ) +

  txt(
    8.25, 4.08,
    "MEGAHIT\n→ AMRFinderPlus",
    size = 3.05
  ) +

  # ARG-bearing contig taxonomy
  arrow_seg(
    8.25, 3.55,
    8.25, 3.05,
    colour = ORANGE
  ) +

  box(
    6.55, 9.95, 1.90, 3.00,
    fill = LIGHT_ORANGE,
    colour = ORANGE,
    linewidth = 0.9
  ) +

  txt(
    8.25, 2.45,
    "ARG-bearing contig\ntaxonomic assignment",
    size = 3.05,
    face = "bold"
  ) +

  # Integration arrows
  arrow_seg(
    1.75, 3.55,
    4.15, 1.35,
    colour = BLUE
  ) +

  arrow_seg(
    5.00, 3.55,
    5.00, 1.35,
    colour = RED
  ) +

  arrow_seg(
    8.25, 1.85,
    5.85, 1.35,
    colour = ORANGE
  ) +

  # Integrated analysis
  box(
    3.2, 6.8, 0.35, 1.35,
    fill = LIGHT_PURPLE,
    colour = PURPLE,
    linewidth = 1.1
  ) +

  txt(
    5, 0.85,
    "Integrated analysis",
    size = 3.8,
    face = "bold",
    colour = PURPLE
  ) +

  coord_cartesian(
    xlim = c(0, 10),
    ylim = c(0.1, 10.1),
    clip = "off"
  ) +

  base_theme


# ============================================================
# PANEL C
# ANALYTICAL FRAMEWORK
# ============================================================

pC <- ggplot() +

  # Main comparison
  box(
    2.0, 8.0, 8.65, 9.65,
    fill = LIGHT_PURPLE,
    colour = PURPLE,
    linewidth = 1
  ) +

  txt(
    5, 9.15,
    "Matched PMA vs non-PMA comparison",
    size = 3.8,
    face = "bold"
  ) +

  arrow_seg(
    5, 8.60,
    5, 8.05,
    colour = PURPLE
  ) +

  # Microbiome
  box(
    0.15, 3.15, 5.95, 8.00,
    fill = LIGHT_BLUE,
    colour = BLUE,
    linewidth = 1
  ) +

  txt(
    1.65, 7.60,
    "Microbiome",
    size = 4.0,
    face = "bold",
    colour = BLUE
  ) +

  txt(
    1.65, 6.90,
    "Diversity\ncomposition\nrobust taxa",
    size = 3.1
  ) +

  txt(
    1.65, 6.20,
    "Figures 2–3",
    size = 3.25,
    face = "bold"
  ) +

  # Resistome
  box(
    3.50, 6.50, 5.95, 8.00,
    fill = LIGHT_RED,
    colour = RED,
    linewidth = 1
  ) +

  txt(
    5, 7.60,
    "Resistome",
    size = 4.0,
    face = "bold",
    colour = RED
  ) +

  txt(
    5, 6.90,
    "Burden • breadth\npersistence\ncross-tool evidence",
    size = 3.05
  ) +

  txt(
    5, 6.20,
    "Figures 4–5",
    size = 3.25,
    face = "bold"
  ) +

  # ARG taxon
  box(
    6.85, 9.85, 5.95, 8.00,
    fill = LIGHT_ORANGE,
    colour = ORANGE,
    linewidth = 1
  ) +

  txt(
    8.35, 7.60,
    "ARG–taxon",
    size = 4.0,
    face = "bold",
    colour = "#D97706"
  ) +

  txt(
    8.35, 6.90,
    "ARG-bearing contigs\nrecurrent taxonomic\nassociations",
    size = 3.0
  ) +

  txt(
    8.35, 6.20,
    "Figure 6",
    size = 3.25,
    face = "bold"
  ) +

  # Converge
  arrow_seg(
    1.65, 5.90,
    4.20, 4.45,
    colour = BLUE
  ) +

  arrow_seg(
    5, 5.90,
    5, 4.45,
    colour = RED
  ) +

  arrow_seg(
    8.35, 5.90,
    5.80, 4.45,
    colour = ORANGE
  ) +

  # Coupling
  box(
    2.15, 7.85, 3.10, 4.45,
    fill = LIGHT_PURPLE,
    colour = PURPLE,
    linewidth = 1.1
  ) +

  txt(
    5, 3.92,
    "Microbiome–resistome coupling",
    size = 3.9,
    face = "bold",
    colour = PURPLE
  ) +

  txt(
    5, 3.42,
    "Figure 7",
    size = 3.25,
    face = "bold"
  ) +

  arrow_seg(
    5, 3.05,
    5, 2.50,
    colour = GREEN
  ) +

  # Final interpretation
  box(
    1.15, 8.85, 0.85, 2.45,
    fill = LIGHT_GREEN,
    colour = GREEN,
    linewidth = 1.1
  ) +

  txt(
    5, 1.92,
    "Viability-informed microbial and\nresistance signatures in the ISS environment",
    size = 3.65,
    face = "bold",
    colour = "#287A3A"
  ) +

  txt(
    5, 1.25,
    "PMA indicates intact-cell-associated DNA;\nnot direct evidence of activity or gene expression",
    size = 2.75,
    colour = EDGE
  ) +

  coord_cartesian(
    xlim = c(0, 10),
    ylim = c(0.5, 10),
    clip = "off"
  ) +

  base_theme


# ============================================================
# COMPILE
# ============================================================

fig1 <- pA | pB | pC

fig1 <- fig1 +
  plot_layout(
    widths = c(1, 1.12, 1.12)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 18,
        colour = "black"
      )
    )
  )

# ============================================================
# SAVE
# ============================================================

ggsave(
  file.path(
    outdir,
    "Figure1_FINAL_VECTOR.pdf"
  ),
  fig1,
  width = 17,
  height = 6.8,
  device = cairo_pdf
)

ggsave(
  file.path(
    outdir,
    "Figure1_FINAL_VECTOR.png"
  ),
  fig1,
  width = 17,
  height = 6.8,
  dpi = 500
)

cat("\n========================================\n")
cat("REVISED FIGURE 1 CREATED\n")
cat("========================================\n")
cat("A: Matched study design\n")
cat("B: Bioinformatics workflow\n")
cat("C: Analytical framework\n")
cat("PDF: figures/figure1/Figure1_FINAL_VECTOR.pdf\n")
cat("PNG: figures/figure1/Figure1_FINAL_VECTOR.png\n")
cat("========================================\n")
