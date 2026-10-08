
# ============================================================
# Figure 4C
# Resistance-class prevalence after PMA
# DeepARG, NASA MT2
# ============================================================

.libPaths(c(path.expand("~/R/4.5"), .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

INFILE <- paste0(
  "analysis_ready/deeparg/figure4/",
  "DeepARG_category_paired_results.tsv"
)

OUTDIR <- "figures/figure4"

d <- read.delim(
  INFILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Select biologically informative classes:
# all FDR-significant classes +
# major highly prevalent PMA classes
# ------------------------------------------------------------

persistent_add <- c(
  "tetracycline",
  "bacitracin",
  "aminoglycoside",
  "peptide"
)

plotdat <- d %>%
  filter(
    FDR_BH < 0.05 |
      ARG_category %in% persistent_add
  ) %>%
  mutate(
    Significant = FDR_BH < 0.05,
    Display_category = case_when(
      ARG_category == "antibacterial_free_fatty_acids" ~
        "antibacterial free\nfatty acids",
      ARG_category == "fusidic-acid" ~
        "fusidic acid",
      TRUE ~ gsub("_", " ", ARG_category)
    ),

    Label = ifelse(
      Significant,
      paste0(Display_category, "*"),
      Display_category
    )
  )

plotdat <- tibble::as_tibble(plotdat)

cat("\n===== FIGURE 4C CLASSES =====\n")

print(
  plotdat %>%
    select(
      ARG_category,
      `Prevalence_non-PMA`,
      Prevalence_PMA,
      FDR_BH,
      Significant
    ) %>%
    arrange(desc(Prevalence_PMA)),
  n = Inf
)

# ------------------------------------------------------------
# Order by PMA prevalence
# ------------------------------------------------------------

plotdat <- plotdat %>%
  arrange(Prevalence_PMA)

plotdat$Label <- factor(
  plotdat$Label,
  levels = plotdat$Label
)

long <- plotdat %>%
  select(
    Label,
    `Prevalence_non-PMA`,
    Prevalence_PMA
  ) %>%
  pivot_longer(
    cols = c(
      `Prevalence_non-PMA`,
      Prevalence_PMA
    ),
    names_to = "STATUS",
    values_to = "Prevalence"
  ) %>%
  mutate(
    STATUS = recode(
      STATUS,
      "Prevalence_non-PMA" = "non-PMA",
      "Prevalence_PMA" = "PMA"
    )
  )

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

cols <- c(
  "non-PMA" = "#1479E8",
  "PMA" = "#FF4B4B"
)

p <- ggplot() +

  geom_segment(
    data = plotdat,
    aes(
      x = `Prevalence_non-PMA` * 100,
      xend = Prevalence_PMA * 100,
      y = Label,
      yend = Label
    ),
    color = "grey72",
    linewidth = 0.8
  ) +

  geom_point(
    data = long,
    aes(
      x = Prevalence * 100,
      y = Label,
      color = STATUS
    ),
    size = 3.2
  ) +

  scale_color_manual(
    values = cols,
    breaks = c("non-PMA", "PMA"),
    name = NULL
  ) +

  scale_x_continuous(
    limits = c(0, 103),
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.01, 0.02))
  ) +

  labs(
    x = "Samples with detected ARG class",
    y = NULL,
    title = "C"
  ) +

  theme_classic(base_size = 11) +

  theme(
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(
      size = 9.4,
      color = "black"
    ),
    axis.title = element_text(color = "black"),

    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(
      face = "bold",
      size = 10
    ),
    legend.key.width = grid::unit(0.65, "cm"),

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
    "Figure4C_ARG_class_persistence.pdf"
  ),
  p,
  width = 5.2,
  height = 5.0,
  units = "in"
)

ggsave(
  file.path(
    OUTDIR,
    "Figure4C_ARG_class_persistence.png"
  ),
  p,
  width = 5.2,
  height = 5.0,
  units = "in",
  dpi = 400
)

cat("\n* indicates paired abundance FDR < 0.05\n")
cat("Saved Figure 4C.\n")

