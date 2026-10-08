
# ============================================================
# Figure 4A
# Paired total DeepARG 16S-normalized ARG abundance
# NASA MT2: non-PMA vs PMA
# ============================================================

.libPaths(c(path.expand("~/R/4.5"), .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

INFILE <- "analysis_ready/deeparg/figure4/DeepARG_sample_resistome_metrics.tsv"
OUTDIR <- "figures/figure4"

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

d <- read.delim(
  INFILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n===== INPUT COLUMNS =====\n")
print(names(d))

# Find total 16S-normalized abundance column robustly.
candidate_cols <- grep(
  "16S|16s",
  names(d),
  value = TRUE
)

cat("\n16S candidate columns:\n")
print(candidate_cols)

# Prefer the previously generated total abundance column.
if ("DeepARG_total_16S" %in% names(d)) {
  abundance_col <- "DeepARG_total_16S"
} else {
  abundance_col <- candidate_cols[1]
}

cat("\nUsing abundance column:", abundance_col, "\n")

plotdat <- d %>%
  select(
    Pair_ID,
    STATUS,
    all_of(abundance_col)
  )

names(plotdat)[3] <- "Abundance"

plotdat$STATUS <- factor(
  plotdat$STATUS,
  levels = c("non-PMA", "PMA")
)

# ------------------------------------------------------------
# Paired Wilcoxon test
# ------------------------------------------------------------

wide <- plotdat %>%
  tidyr::pivot_wider(
    names_from = STATUS,
    values_from = Abundance
  )

wt <- wilcox.test(
  wide$PMA,
  wide$`non-PMA`,
  paired = TRUE,
  exact = FALSE
)

med_non <- median(wide$`non-PMA`, na.rm = TRUE)
med_pma <- median(wide$PMA, na.rm = TRUE)

cat("\n===== FIGURE 4A STATISTICS =====\n")
cat("Pairs:", nrow(wide), "\n")
cat("Median non-PMA:", med_non, "\n")
cat("Median PMA:", med_pma, "\n")
cat("Paired Wilcoxon P:", wt$p.value, "\n")

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

cols <- c(
  "non-PMA" = "#1479E8",
  "PMA" = "#FF4B4B"
)

ymax <- max(plotdat$Abundance, na.rm = TRUE)

p <- ggplot(
  plotdat,
  aes(
    x = STATUS,
    y = Abundance,
    group = Pair_ID
  )
) +

  geom_line(
    color = "grey80",
    linewidth = 0.45,
    alpha = 0.85
  ) +

  geom_point(
    aes(color = STATUS),
    size = 2.5,
    alpha = 0.9
  ) +

  stat_summary(
    aes(group = STATUS),
    fun = median,
    geom = "crossbar",
    width = 0.42,
    linewidth = 0.55,
    color = "black"
  ) +

  annotate(
    "segment",
    x = 1,
    xend = 2,
    y = ymax * 1.08,
    yend = ymax * 1.08,
    linewidth = 0.45
  ) +

  annotate(
    "text",
    x = 1.5,
    y = ymax * 1.14,
    label = expression(italic(P) == 0.029),
    size = 3.6
  ) +

  scale_color_manual(
    values = cols,
    guide = "none"
  ) +

  scale_y_continuous(
    expand = expansion(mult = c(0.03, 0.20))
  ) +

  labs(
    x = NULL,
    y = "Total ARG abundance\n(16S-normalized)",
    title = "A"
  ) +

  theme_classic(base_size = 11) +

  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0
    ),
    plot.margin = margin(8, 10, 8, 8)
  )

ggsave(
  file.path(OUTDIR, "Figure4A_ARG_burden.pdf"),
  p,
  width = 4.3,
  height = 4.2,
  units = "in"
)

ggsave(
  file.path(OUTDIR, "Figure4A_ARG_burden.png"),
  p,
  width = 4.3,
  height = 4.2,
  units = "in",
  dpi = 400
)

cat("\nSaved:\n")
cat(file.path(OUTDIR, "Figure4A_ARG_burden.pdf"), "\n")
cat(file.path(OUTDIR, "Figure4A_ARG_burden.png"), "\n")

