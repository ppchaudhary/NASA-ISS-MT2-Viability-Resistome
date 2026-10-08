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
# Primary association
# ============================================================

ct <- cor.test(
  dat$Resistome_BrayCurtis,
  dat$Abs_Delta_ARG_burden,
  method = "spearman",
  exact = FALSE
)

rho  <- unname(ct$estimate)
pval <- ct$p.value

# Previously validated depth-adjusted rank-regression result
depth_beta <- 0.584
depth_p    <- 0.00119

cat("\n===== FIGURE 7C =====\n")
cat("N =", nrow(dat), "\n")
cat("Spearman rho =", rho, "\n")
cat("P =", pval, "\n")
cat("Depth-adjusted beta =", depth_beta, "\n")
cat("Depth-adjusted P =", depth_p, "\n\n")

# ============================================================
# Plot
# ============================================================

p <- ggplot(
  dat,
  aes(
    x = Resistome_BrayCurtis,
    y = Abs_Delta_ARG_burden
  )
) +

  # LOESS is visualization only;
  # formal inference is Spearman correlation.
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

  # Same combined statistics box as 7A and 7B
  annotate(
    "label",
    x = 0.90,
    y = max(dat$Abs_Delta_ARG_burden, na.rm = TRUE) * 0.97,
    label = paste0(
      "Spearman \u03c1 = ", sprintf("%.2f", rho),
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
    limits = c(0, 0.96),
    breaks = seq(0, 0.8, 0.2),
    expand = expansion(mult = c(0.01, 0.02))
  ) +

  scale_y_continuous(
    limits = c(
      0,
      max(dat$Abs_Delta_ARG_burden, na.rm = TRUE) * 1.05
    ),
    expand = expansion(mult = c(0.01, 0.02))
  ) +

  labs(
    x = "Resistome Bray\u2013Curtis distance\n(PMA vs non-PMA)",
    y = "Absolute change in normalized ARG burden\n|PMA \u2212 non-PMA|"
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
    "Figure7C_resistome_ARGburden.pdf"
  ),
  p,
  width = 6.2,
  height = 5.5,
  device = cairo_pdf
)

ggsave(
  file.path(
    outdir,
    "Figure7C_resistome_ARGburden.png"
  ),
  p,
  width = 6.2,
  height = 5.5,
  dpi = 500
)

cat("Saved Figure7C_resistome_ARGburden.pdf/png\n")
