
# ============================================================
# Figure 4D
# Recurrently detected ARG groups after PMA
# DeepARG, NASA MT2
# Inclusion: >=75% PMA prevalence
# ============================================================

.libPaths(c(path.expand("~/R/4.5"), .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggtext)
})

INFILE <- paste0(
  "analysis_ready/deeparg/figure4/",
  "DeepARG_ARG_group_persistence.tsv"
)

OUTDIR <- "figures/figure4"

d <- read.delim(
  INFILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

plotdat <- d %>%
  filter(Prevalence_PMA >= 0.75) %>%
  mutate(
    ARG_label = case_when(
      ARG_group == "RPOB2" ~ "rpoB2",
      ARG_group == "BACA"  ~ "bacA",
      ARG_group == "NORA"  ~ "norA",
      ARG_group == "MARR"  ~ "marR",
      ARG_group == "MSRA"  ~ "msrA",
      ARG_group == "FMTC"  ~ "fmtC",
      ARG_group == "UGD"   ~ "ugd",
      ARG_group == "TETM"  ~ "tetM",
      ARG_group == "DFRC"  ~ "dfrC",
      ARG_group == "ACRB"  ~ "acrB",

      ARG_group == "MULTIDRUG_ABC_TRANSPORTER" ~
        "Multidrug ABC\ntransporter",

      ARG_group ==
        "MAJOR_FACILITATOR_SUPERFAMILY_TRANSPORTER" ~
        "Major facilitator\ntransporter",

      ARG_group ==
        "EMRB-QACA_FAMILY_MAJOR_FACILITATOR_TRANSPORTER" ~
        "EmrB/QacA-family\ntransporter",

      TRUE ~ ARG_group
    )
  ) %>%
  arrange(Prevalence_PMA) %>%
  as_tibble()

cat("\n===== FIGURE 4D ARG GROUPS =====\n")

print(
  plotdat %>%
    select(
      ARG_group,
      `Prevalence_non-PMA`,
      Prevalence_PMA,
      Median_nonPMA,
      Median_PMA,
      FDR_BH,
      Direction
    ),
  n = Inf
)

# ------------------------------------------------------------
# Order
# ------------------------------------------------------------

plotdat$ARG_label <- factor(
  plotdat$ARG_label,
  levels = plotdat$ARG_label
)

# Publication-style labels:
# italicize gene symbols; retain generic transporter labels
# as regular text.
gene_labels <- c(
  "rpoB2", "bacA", "norA", "marR", "msrA",
  "fmtC", "ugd", "tetM", "dfrC", "acrB"
)

label_map <- setNames(
  as.character(levels(plotdat$ARG_label)),
  levels(plotdat$ARG_label)
)

for (g in gene_labels) {
  if (g %in% names(label_map)) {
    label_map[g] <- paste0("italic('", g, "')")
  }
}

long <- plotdat %>%
  select(
    ARG_label,
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
      y = ARG_label,
      yend = ARG_label
    ),
    color = "grey72",
    linewidth = 0.8
  ) +

  geom_point(
    data = long,
    aes(
      x = Prevalence * 100,
      y = ARG_label,
      color = STATUS
    ),
    size = 3.2
  ) +

  scale_color_manual(
    values = cols,
    breaks = c("non-PMA", "PMA"),
    name = NULL
  ) +

  scale_y_discrete(
    labels = function(x) {
      sapply(
        x,
        function(z) {
          if (z %in% gene_labels) {
            paste0("<i>", z, "</i>")
          } else {
            gsub("\\n", "<br>", z)
          }
        }
      )
    }
  ) +

  scale_x_continuous(
    limits = c(0, 103),
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.01, 0.02))
  ) +

  labs(
    x = "Samples with detected ARG group",
    y = NULL,
    title = "D"
  ) +

  theme_classic(base_size = 11) +

  theme(
    axis.text = element_text(color = "black"),
    axis.text.y = ggtext::element_markdown(
      color = "black",
      size = 9.4,
      lineheight = 0.95
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
    "Figure4D_persistent_ARGs.pdf"
  ),
  p,
  width = 5.2,
  height = 5.0,
  units = "in"
)

ggsave(
  file.path(
    OUTDIR,
    "Figure4D_persistent_ARGs.png"
  ),
  p,
  width = 5.2,
  height = 5.0,
  units = "in",
  dpi = 400
)

cat("\nThreshold: >=75% PMA prevalence\n")
cat("Saved Figure 4D.\n")

