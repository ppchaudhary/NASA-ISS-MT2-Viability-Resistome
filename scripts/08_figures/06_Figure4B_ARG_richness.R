
# ============================================================
# Figure 4B
# Paired DeepARG ARG-group richness
# NASA MT2: non-PMA vs PMA
# ============================================================

.libPaths(c(path.expand("~/R/4.5"), .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

INFILE <- "analysis_ready/deeparg/figure4/DeepARG_sample_resistome_metrics.tsv"
OUTDIR <- "figures/figure4"

d <- read.delim(
  INFILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n===== INPUT COLUMNS =====\n")
print(names(d))

# Find richness column
if ("ARG_subtype_richness" %in% names(d)) {
  richness_col <- "ARG_subtype_richness"
} else {
  richness_candidates <- grep(
    "richness",
    names(d),
    value = TRUE,
    ignore.case = TRUE
  )
  richness_col <- richness_candidates[1]
}

cat("\nUsing richness column:", richness_col, "\n")

plotdat <- d %>%
  select(
    Pair_ID,
    STATUS,
    all_of(richness_col)
  )

names(plotdat)[3] <- "Richness"

plotdat$STATUS <- factor(
  plotdat$STATUS,
  levels = c("non-PMA", "PMA")
)

wide <- plotdat %>%
  pivot_wider(
    names_from = STATUS,
    values_from = Richness
  )

wt <- wilcox.test(
  wide$PMA,
  wide$`non-PMA`,
  paired = TRUE,
  exact = FALSE
)

cat("\n===== FIGURE 4B STATISTICS =====\n")
cat("Pairs:", nrow(wide), "\n")
cat(
  "Median non-PMA:",
  median(wide$`non-PMA`),
  "\n"
)
cat(
  "Median PMA:",
  median(wide$PMA),
  "\n"
)
cat(
  "PMA lower:",
  sum(wide$PMA < wide$`non-PMA`),
  "of", nrow(wide), "\n"
)
cat(
  "Paired Wilcoxon P:",
  wt$p.value,
  "\n"
)

cols <- c(
  "non-PMA" = "#1479E8",
  "PMA" = "#FF4B4B"
)

ymax <- max(plotdat$Richness, na.rm = TRUE)

p <- ggplot(
  plotdat,
  aes(
    x = STATUS,
    y = Richness,
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
    label = expression(
      italic(P) < 10^-6
    ),
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
    y = "ARG-group richness",
    title = "B"
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
  file.path(
    OUTDIR,
    "Figure4B_ARG_richness.pdf"
  ),
  p,
  width = 4.3,
  height = 4.2,
  units = "in"
)

ggsave(
  file.path(
    OUTDIR,
    "Figure4B_ARG_richness.png"
  ),
  p,
  width = 4.3,
  height = 4.2,
  units = "in",
  dpi = 400
)

cat("\nSaved Figure 4B.\n")

