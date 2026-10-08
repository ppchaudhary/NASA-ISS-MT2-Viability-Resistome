
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

dir.create("figures/figure5", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Read AMRFinder class prevalence
# ============================================================

x <- read_tsv(
  "analysis_ready/amr_concordance/AMRFinder_class_prevalence_wide.tsv",
  show_col_types = FALSE
)

# Keep classes detected in at least 10% of samples in either group.
# This avoids filling the main panel with extremely rare classes.
pdat <- x %>%
  filter(
    Prevalence_PMA >= 0.10 |
    `Prevalence_non-PMA` >= 0.10
  ) %>%
  mutate(
    Class_label = str_to_sentence(Class),
    Class_label = str_replace_all(
      Class_label,
      "/",
      " / "
    ),
    Class_label = str_wrap(
      Class_label,
      width = 28
    ),
    maxprev = pmax(
      Prevalence_PMA,
      `Prevalence_non-PMA`
    )
  ) %>%
  arrange(maxprev) %>%
  mutate(
    Class_label = factor(
      Class_label,
      levels = Class_label
    )
  )

cat("\n===== FIGURE 5B CLASSES =====\n")
print(
  pdat %>%
    select(
      Class,
      Prevalence_PMA,
      `Prevalence_non-PMA`
    )
)

# ============================================================
# Plot
# ============================================================

p <- ggplot(pdat) +

  geom_segment(
    aes(
      x = 100 * `Prevalence_non-PMA`,
      xend = 100 * Prevalence_PMA,
      y = Class_label,
      yend = Class_label
    ),
    color = "grey72",
    linewidth = 0.8
  ) +

  geom_point(
    aes(
      x = 100 * `Prevalence_non-PMA`,
      y = Class_label,
      color = "non-PMA"
    ),
    size = 3.2
  ) +

  geom_point(
    aes(
      x = 100 * Prevalence_PMA,
      y = Class_label,
      color = "PMA"
    ),
    size = 3.2
  ) +

  scale_color_manual(
    values = c(
      "non-PMA" = "#1479E8",
      "PMA" = "#FF4B4B"
    ),
    breaks = c("non-PMA", "PMA"),
    name = NULL
  ) +

  scale_x_continuous(
    breaks = seq(0, 100, 25),
    expand = expansion(mult = c(0, 0))
  ) +

  coord_cartesian(
    xlim = c(-3, 102),
    clip = "off"
  ) +

  labs(
    x = "Samples with assembled AMR class (%)",
    y = NULL
  ) +

  theme_classic(base_size = 11) +

  theme(
    axis.text = element_text(
      color = "black"
    ),

    axis.text.y = element_text(
      size = 9
    ),

    axis.title.x = element_text(
      size = 10.5
    ),

    legend.position = "top",

    legend.text = element_text(
      face = "bold",
      size = 10
    ),

    legend.key.width = unit(
      0.8,
      "cm"
    ),

    plot.margin = margin(
      8, 10, 8, 8
    )
  )

# ============================================================
# Save
# ============================================================

ggsave(
  "figures/figure5/Figure5B_AMRFinder_classes.pdf",
  p,
  width = 6.1,
  height = 5.0
)

ggsave(
  "figures/figure5/Figure5B_AMRFinder_classes.png",
  p,
  width = 6.1,
  height = 5.0,
  dpi = 400
)

cat("\nSaved Figure 5B.\n")

