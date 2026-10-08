library(tidyverse)

infile <- paste0(
  "analysis_ready/integrated_analysis/",
  "paired_ARG_burden_changes.tsv"
)

outdir <- "figures/figure7"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dat <- read.delim(
  infile,
  check.names = FALSE
)

# ============================================================
# Recalculate primary association from plotting data
# ============================================================

ct <- cor.test(
  dat$Microbiome_BrayCurtis,
  dat$Delta_ARG_burden,
  method = "spearman",
  exact = FALSE
)

rho <- unname(ct$estimate)
pval <- ct$p.value

# Validated depth-sensitivity rank-regression result
depth_beta <- -0.425
depth_p <- 0.0169

cat("\n===== FIGURE 7B =====\n")
cat("N =", nrow(dat), "\n")
cat("Spearman rho =", rho, "\n")
cat("P =", pval, "\n")
cat("Depth-adjusted beta =", depth_beta, "\n")
cat("Depth-adjusted P =", depth_p, "\n\n")

# ============================================================
# Labels
# ============================================================

stat_label <- paste0(
  "Spearman \u03c1 = ",
  sprintf("%.2f", rho),
  "\nP = ",
  format.pval(
    pval,
    digits = 2,
    eps = 0.001
  )
)

depth_label <- paste0(
  "Depth-adjusted \u03b2 = ",
  sprintf("%.2f", depth_beta),
  ", P = ",
  format.pval(
    depth_p,
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
    y = Delta_ARG_burden
  )
) +

  # Zero = no difference in normalized ARG burden
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.55,
    color = "grey55"
  ) +

  # LOESS = visualization only
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
    x = 0.77,
    y = 2.85,
    label = paste0(
      "Spearman ρ = ", sprintf("%.2f", rho),
      "\nP = ", format.pval(pval, digits = 2, eps = 0.001),
      "\nDepth-adjusted P = ",
      format.pval(depth_p, digits = 2, eps = 0.001)
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

  labs(
    x = "Microbiome Bray\u2013Curtis distance\n(PMA vs non-PMA)",
    y = "\u0394 normalized ARG burden\n(PMA \u2212 non-PMA)"
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
    "Figure7B_microbiome_ARGburden.pdf"
  ),
  p,
  width = 6.2,
  height = 5.5,
  device = cairo_pdf
)

ggsave(
  file.path(
    outdir,
    "Figure7B_microbiome_ARGburden.png"
  ),
  p,
  width = 6.2,
  height = 5.5,
  dpi = 500
)

cat(
  "Saved Figure7B_microbiome_ARGburden.pdf/png\n"
)
