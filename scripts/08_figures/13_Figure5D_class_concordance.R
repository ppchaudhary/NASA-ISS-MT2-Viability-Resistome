
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
})

dir.create(
  "figures/figure5",
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Read PMA class prevalence
# ============================================================

x <- read_tsv(
  "analysis_ready/amr_concordance/class_prevalence_four_methods.tsv",
  show_col_types = FALSE
) %>%
  filter(STATUS == "PMA")

# ============================================================
# Order classes by convergence across read-based methods
# ============================================================

summary_dat <- x %>%
  filter(Method %in% c("DeepARG", "ShortBRED", "RGI")) %>%
  group_by(Class) %>%
  summarise(
    Methods_50 = sum(Prevalence >= 0.50),
    Mean_prevalence = mean(Prevalence),
    .groups = "drop"
  ) %>%
  arrange(
    desc(Methods_50),
    desc(Mean_prevalence)
  )

class_order <- summary_dat$Class

# ============================================================
# Plot data
# ============================================================

pdat <- x %>%
  mutate(
    Method = factor(
      Method,
      levels = c(
        "DeepARG",
        "ShortBRED",
        "RGI",
        "AMRFinderPlus"
      )
    ),
    Class = factor(
      Class,
      levels = rev(class_order)
    ),
    Label = paste0(
      round(Prevalence * 100),
      "%"
    )
  )

# ============================================================
# Heatmap
# ============================================================

p <- ggplot(
  pdat,
  aes(
    x = Method,
    y = Class,
    fill = Prevalence
  )
) +

  geom_tile(
    color = "white",
    linewidth = 1
  ) +

  geom_text(
    aes(label = Label),
    size = 3.2,
    color = "black"
  ) +

  # Separate assembly-based AMRFinderPlus
  geom_vline(
    xintercept = 3.5,
    linetype = "dashed",
    linewidth = 0.6,
    color = "grey45"
  ) +

  scale_x_discrete(
    limits = c(
      "DeepARG",
      "ShortBRED",
      "RGI",
      "AMRFinderPlus"
    ),
    drop = FALSE
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

  labs(
    x = NULL,
    y = NULL
  ) +

  theme_classic(
    base_size = 11
  ) +

  theme(
    axis.text.x = element_text(
      color = "black",
      face = "bold",
      size = 9.5
    ),

    axis.text.y = element_text(
      color = "black",
      face = "plain",
      size = 9.5
    ),

    axis.ticks = element_blank(),

    legend.position = "bottom",

    legend.title = element_text(
      size = 9
    ),

    legend.text = element_text(
      size = 8
    ),

    legend.key.width = grid::unit(
      1.4,
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
  "figures/figure5/Figure5D_class_concordance.pdf",
  p,
  width = 6.5,
  height = 5.2
)

ggsave(
  "figures/figure5/Figure5D_class_concordance.png",
  p,
  width = 6.5,
  height = 5.2,
  dpi = 400
)

cat("\n===== FIGURE 5D CLASS ORDER =====\n")
print(summary_dat)

cat("\nSaved Figure 5D.\n")

