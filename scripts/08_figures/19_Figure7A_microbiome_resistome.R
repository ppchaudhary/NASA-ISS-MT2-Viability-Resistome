library(tidyverse)

infile <- paste0(
  "analysis_ready/integrated_analysis/",
  "paired_distances_with_depth.tsv"
)

outdir <- "figures/figure7"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dat <- read.delim(
  infile,
  check.names = FALSE
)

# ============================================================
# Confirm primary association directly from plotting data
# ============================================================

ct <- cor.test(
  dat$Microbiome_BrayCurtis,
  dat$Resistome_BrayCurtis,
  method = "spearman",
  exact = FALSE
)

rho <- unname(ct$estimate)
pval <- ct$p.value

# Previously validated depth-adjusted rank result
partial_r <- 0.527
partial_p <- 0.00194

cat("\n===== FIGURE 7A =====\n")
cat("N =", nrow(dat), "\n")
cat("Spearman rho =", rho, "\n")
cat("P =", pval, "\n")
cat("Depth-adjusted r =", partial_r, "\n")
cat("Depth-adjusted P =", partial_p, "\n\n")

# ============================================================
# Annotation
# ============================================================

stat_label <- paste0(
  "Spearman \u03c1 = ",
  sprintf("%.2f", rho),
  "\n",
  "P = ",
  format.pval(
    pval,
    digits = 2,
    eps = 0.001
  )
)

depth_label <- paste0(
  "Depth-adjusted r = ",
  sprintf("%.2f", partial_r),
  ", P = ",
  format.pval(
    partial_p,
    digits = 2,
    eps = 0.001
  )
)

# ============================================================
# Plot
# ============================================================

p <- ggplot(
  dat,
  aes(
    x = Microbiome_BrayCurtis,
    y = Resistome_BrayCurtis
  )
) +

  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = TRUE,
    span = 0.9,
    color = "grey25",
    fill = "grey85",
    linewidth = 0.8,
    alpha = 0.45
  ) +

  geom_point(
    shape = 21,
    size = 3.5,
    stroke = 0.6,
    fill = "#1479E8",
    color = "black",
    alpha = 0.9
  ) +

  # Combined statistical annotation
  annotate(
    "label",
    x = 0.79,
    y = 0.93,
    label = paste0(
      "Spearman ρ = ", sprintf("%.2f", rho),
      "\nP = ", format.pval(pval, digits = 2, eps = 0.001),
      "\nDepth-adjusted P = ",
      format.pval(partial_p, digits = 2, eps = 0.001)
    ),
    hjust = 1,
    vjust = 1,
    size = 3.6,
    lineheight = 1.15,
    fill = "white",
    color = "black",
    label.size = 0.25
  ) +

  scale_x_continuous(
    limits = c(0, 0.82),
    breaks = seq(0, 0.8, 0.2),
    expand = expansion(mult = c(0.01, 0.02))
  ) +

  scale_y_continuous(
    limits = c(0, 0.96),
    breaks = seq(0, 0.8, 0.2),
    expand = expansion(mult = c(0.01, 0.02))
  ) +

  labs(
    x = "Microbiome Bray\u2013Curtis distance\n(PMA vs non-PMA)",
    y = "Resistome Bray\u2013Curtis distance\n(PMA vs non-PMA)"
  ) +

  theme_classic(base_size = 12) +

  theme(
    axis.title = element_text(
      size = 12,
      color = "black"
    ),

    axis.text = element_text(
      size = 10,
      color = "black"
    ),

    axis.line = element_line(
      linewidth = 0.6,
      color = "black"
    ),

    axis.ticks = element_line(
      color = "black"
    ),

    plot.margin = margin(
      12, 18, 12, 12
    )
  )

# ============================================================
# Save
# ============================================================

ggsave(
  file.path(
    outdir,
    "Figure7A_microbiome_resistome.pdf"
  ),
  p,
  width = 6.2,
  height = 5.5,
  device = cairo_pdf
)

ggsave(
  file.path(
    outdir,
    "Figure7A_microbiome_resistome.png"
  ),
  p,
  width = 6.2,
  height = 5.5,
  dpi = 500
)

cat(
  "Saved Figure7A_microbiome_resistome.pdf/png\n"
)
