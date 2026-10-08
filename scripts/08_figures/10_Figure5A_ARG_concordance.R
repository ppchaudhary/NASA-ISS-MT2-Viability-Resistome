
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggtext)
})

dir.create("figures/figure5", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Read harmonized concordance table
# ============================================================

x <- read_tsv(
  "analysis_ready/amr_concordance/CARD_candidate_concordance.tsv",
  show_col_types = FALSE
)

# Order by cross-method support
arg_order <- x %>%
  arrange(desc(Methods_PMA50),
          desc(Mean_PMA_prevalence)) %>%
  pull(ARG)

# ============================================================
# Long format
# ============================================================

pdat <- x %>%
  select(
    ARG,
    DeepARG_PMA,
    ShortBRED_PMA,
    RGI_PMA
  ) %>%
  pivot_longer(
    -ARG,
    names_to = "Method",
    values_to = "Prevalence"
  ) %>%
  mutate(
    Method = recode(
      Method,
      DeepARG_PMA = "DeepARG",
      ShortBRED_PMA = "ShortBRED",
      RGI_PMA = "RGI"
    ),
    Method = factor(
      Method,
      levels = c("DeepARG", "ShortBRED", "RGI")
    ),
    ARG = factor(
      ARG,
      levels = rev(arg_order)
    ),
    Percent = 100 * Prevalence,
    Label = sprintf("%.0f%%", Percent)
  )

# ============================================================
# Plot
# ============================================================

p <- ggplot(
  pdat,
  aes(
    x = Method,
    y = ARG,
    fill = Prevalence
  )
) +
  geom_tile(
    color = "white",
    linewidth = 1.2
  ) +

  geom_text(
    aes(
      label = Label,
      color = Prevalence >= 0.55
    ),
    size = 3.6,
    fontface = "bold"
  ) +

  scale_color_manual(
    values = c(
      `TRUE` = "white",
      `FALSE` = "black"
    ),
    guide = "none"
  ) +

  scale_fill_gradient2(
    low = "#1479E8",
    mid = "white",
    high = "#FF4B4B",
    midpoint = 0.50,
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = c("0", "25", "50", "75", "100"),
    name = "PMA prevalence (%)"
  ) +

  scale_y_discrete(
    labels = function(z) {
      paste0("<i>", z, "</i>")
    }
  ) +

  labs(
    x = NULL,
    y = NULL
  ) +

  theme_classic(base_size = 11) +

  theme(
    axis.text.x = element_text(color = "black", face = "bold",
      size = 10
    ),

    axis.text.y = ggtext::element_markdown(
      color = "black",
      size = 10
    ),

    axis.ticks = element_blank(),

    legend.title = element_text(
      face = "bold",
      size = 9.5
    ),

    legend.text = element_text(
      size = 9
    ),

    legend.position = "right",

    plot.margin = margin(
      8, 10, 8, 8
    )
  )

# ============================================================
# Save
# ============================================================

ggsave(
  "figures/figure5/Figure5A_ARG_concordance.pdf",
  p,
  width = 5.6,
  height = 4.7
)

ggsave(
  "figures/figure5/Figure5A_ARG_concordance.png",
  p,
  width = 5.6,
  height = 4.7,
  dpi = 400
)

cat("\n===== FIGURE 5A =====\n")
print(
  pdat %>%
    arrange(ARG, Method)
)

cat("\nSaved Figure 5A.\n")

