
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggtext)
})

dir.create(
  "figures/figure5",
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Read evidence table
# ============================================================

x <- read_tsv(
  "analysis_ready/amr_concordance/four_method_candidate_evidence.tsv",
  show_col_types = FALSE
)

# Preserve evidence-based ordering
arg_order <- x %>%
  arrange(
    desc(Four_method_support),
    desc(Read_methods_recurrent),
    desc(DeepARG_PMA_prevalence)
  ) %>%
  pull(ARG)

# ============================================================
# Build evidence matrix
# ============================================================

read_dat <- x %>%
  select(
    ARG,
    DeepARG_recurrent,
    ShortBRED_recurrent,
    RGI_recurrent
  ) %>%
  pivot_longer(
    -ARG,
    names_to = "Method",
    values_to = "Supported"
  ) %>%
  mutate(
    Method = recode(
      Method,
      DeepARG_recurrent = "DeepARG",
      ShortBRED_recurrent = "ShortBRED",
      RGI_recurrent = "RGI"
    ),
    Evidence_type = "Recurrent read-based"
  )

amr_dat <- x %>%
  transmute(
    ARG,
    Method = "AMRFinderPlus",
    Supported = AMRFinder_PMA_detected,
    Evidence_type = "Assembly-based detection"
  )

pdat <- bind_rows(
  read_dat,
  amr_dat
) %>%
  mutate(
    Method = factor(
      Method,
      levels = c(
        "DeepARG",
        "ShortBRED",
        "RGI",
        "AMRFinderPlus"
      ),
      ordered = TRUE
    ),
    ARG = factor(
      ARG,
      levels = rev(arg_order)
    ),
    Evidence_type = factor(
      Evidence_type,
      levels = c(
        "Recurrent read-based",
        "Assembly-based detection"
      )
    )
  )

# ============================================================
# Plot
# ============================================================

# Dedicated legend entries
legend_dat <- data.frame(
  x = NA_real_,
  y = NA_real_,
  Legend = factor(
    c(
      "Recurrent PMA support (>=50%)",
      "Below recurrence threshold",
      "Assembly-based PMA detection"
    ),
    levels = c(
      "Recurrent PMA support (>=50%)",
      "Below recurrence threshold",
      "Assembly-based PMA detection"
    )
  )
)


p <- ggplot(
  pdat,
  aes(
    x = Method,
    y = ARG
  )
) +

  # faint background circles = evaluated but criterion not met
  geom_point(
    data = pdat %>% filter(!Supported),
    shape = 21,
    size = 4.5,
    stroke = 0.7,
    fill = "white",
    color = "grey70",
    aes(shape = "Below recurrence threshold")
  ) +

  # read-based recurrent support
  geom_point(
    data = pdat %>%
      filter(
        Supported,
        Evidence_type == "Recurrent read-based"
      ),
    shape = 21,
    size = 4.7,
    stroke = 0.8,
    fill = "#FF4B4B",
    color = "black",
    aes(shape = "Recurrent PMA support (>=50%)")
  ) +

  # AMRFinder assembly-based support: diamond
  geom_point(
    data = pdat %>%
      filter(
        Supported,
        Evidence_type == "Assembly-based detection"
      ),
    shape = 23,
    size = 5.2,
    stroke = 0.8,
    fill = "#FF4B4B",
    color = "black",
    aes(shape = "Assembly-based PMA detection")
  ) +

  # Separate AMRFinder visually
  geom_vline(
    xintercept = 3.5,
    linetype = "dashed",
    linewidth = 0.55,
    color = "grey60"
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

  geom_point(
    data = legend_dat,
    aes(x = x, y = y, shape = Legend),
    inherit.aes = FALSE,
    size = 4,
    show.legend = TRUE
  ) +

  scale_shape_manual(
    name = NULL,
    values = c(
      "Recurrent PMA support (>=50%)" = 16,
      "Below recurrence threshold" = 1,
      "Assembly-based PMA detection" = 18
    ),
    drop = FALSE
  ) +

  guides(
    shape = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        size = c(4, 4, 4.5),
        color = c("#FF4B4B", "grey65", "#FF4B4B")
      )
    )
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

  theme_classic(
    base_size = 11
  ) +

  theme(
    axis.text.x = element_text(
      color = "black",
      face = "bold",
      size = 9.5
    ),

    axis.text.y = ggtext::element_markdown(
      color = "black",
      size = 10
    ),

    axis.ticks = element_blank(),

    panel.grid.major.y = element_line(
      color = "grey92",
      linewidth = 0.4
    ),

    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 8.5),
    legend.box.margin = margin(t = 3),

    plot.margin = margin(
      8, 10, 8, 8
    )
  )

# ============================================================
# Save
# ============================================================

ggsave(
  "figures/figure5/Figure5C_four_method_evidence.pdf",
  p,
  width = 6.5,
  height = 4.7
)

ggsave(
  "figures/figure5/Figure5C_four_method_evidence.png",
  p,
  width = 6.5,
  height = 4.7,
  dpi = 400
)

cat("\n===== FIGURE 5C EVIDENCE MATRIX =====\n")

print(
  pdat %>%
    arrange(ARG, Method)
)

cat("\nLegend interpretation:\n")
cat("Filled circle = >=50% PMA prevalence in read-based method\n")
cat("Open circle   = criterion not met\n")
cat("Filled diamond = PMA detection by assembly-based AMRFinderPlus\n")

cat("\nSaved Figure 5C.\n")

