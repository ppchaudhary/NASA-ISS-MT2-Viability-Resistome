library(tidyverse)
library(ggtext)

infile <- "analysis_ready/amrfinder/amrfinder_ARG_contig_taxonomy.tsv"
outdir <- "figures/figure6"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

x <- read.delim(infile, check.names = FALSE)

# ============================================================
# Build unique sample-level ARG-taxon observations
# ============================================================

d <- x %>%
  filter(
    !is.na(Taxon_name),
    Taxon_name != "",
    Taxon_name != "unclassified"
  ) %>%
  distinct(
    Run,
    STATUS,
    `Element symbol`,
    Taxon_name
  ) %>%
  rename(
    ARG = `Element symbol`,
    Taxon = Taxon_name
  )

# ============================================================
# Calculate sample prevalence separately for PMA/non-PMA
# 32 samples in each group
# ============================================================

prev <- d %>%
  count(
    STATUS,
    ARG,
    Taxon,
    name = "Samples"
  ) %>%
  mutate(
    Prevalence = Samples / 32
  )

# ============================================================
# Select associations recurrent in >=2 PMA samples
#
# Selection is based ONLY on PMA recurrence.
# Then compare those same associations in non-PMA.
# ============================================================

selected <- prev %>%
  filter(
    STATUS == "PMA",
    Samples >= 2
  ) %>%
  select(ARG, Taxon)

cat("\nSelected recurrent PMA associations:", nrow(selected), "\n")

# Complete both statuses, including true zeros
plot_df <- selected %>%
  crossing(
    STATUS = c("non-PMA", "PMA")
  ) %>%
  left_join(
    prev,
    by = c("ARG", "Taxon", "STATUS")
  ) %>%
  mutate(
    Samples = replace_na(Samples, 0L),
    Prevalence = replace_na(Prevalence, 0),
    Percent = 100 * Prevalence
  )

# ============================================================
# Order associations by PMA recurrence
# ============================================================

order_df <- plot_df %>%
  filter(STATUS == "PMA") %>%
  arrange(
    Samples,
    ARG,
    Taxon
  ) %>%
  mutate(
    Pair = paste(ARG, Taxon, sep = "|||")
  )

pair_levels <- order_df$Pair

plot_df <- plot_df %>%
  mutate(
    Pair = paste(ARG, Taxon, sep = "|||"),
    Pair = factor(Pair, levels = pair_levels),

    # italic ARG only; taxonomic assignment stays regular
    Pair_label = paste0(
      "<i>", ARG, "</i>",
      " \u2013 ",
      Taxon
    ),

    STATUS = factor(
      STATUS,
      levels = c("non-PMA", "PMA")
    ),

    Cell_label = paste0(
      Samples,
      "/32\n",
      round(Percent),
      "%"
    )
  )

# Labels must follow factor ordering
label_lookup <- plot_df %>%
  distinct(Pair, Pair_label) %>%
  arrange(Pair)

pair_labels <- setNames(
  label_lookup$Pair_label,
  as.character(label_lookup$Pair)
)

# ============================================================
# Plot
# ============================================================

p <- ggplot(
  plot_df,
  aes(
    x = STATUS,
    y = Pair,
    fill = Prevalence
  )
) +

  geom_tile(
    color = "white",
    linewidth = 1.2
  ) +

  geom_text(
    aes(label = Cell_label),
    size = 3.7,
    lineheight = 0.9,
    fontface = "bold",
    color = "black"
  ) +

  scale_y_discrete(
    labels = pair_labels
  ) +

  scale_fill_gradient2(
    low = "white",
    mid = "#DDEEFF",
    high = "#1479E8",
    midpoint = 0.5,
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::percent_format(accuracy = 1),
    name = "Sample\nprevalence"
  ) +

  labs(
    x = NULL,
    y = NULL
  ) +

  theme_classic(base_size = 12) +

  theme(
    axis.text.x = element_text(
      color = "black",
      face = "bold",
      size = 11
    ),

    axis.text.y = ggtext::element_markdown(
      color = "black",
      size = 10
    ),

    axis.ticks = element_blank(),

    legend.title = element_text(
      face = "bold",
      size = 10
    ),

    legend.text = element_text(
      size = 9
    ),

    legend.position = "right",

    plot.margin = margin(
      10, 15, 10, 10
    )
  )

# ============================================================
# Save
# ============================================================

ggsave(
  file.path(
    outdir,
    "Figure6B_ARG_taxon_PMA_vs_nonPMA.pdf"
  ),
  p,
  width = 7.5,
  height = 6.0,
  device = cairo_pdf
)

ggsave(
  file.path(
    outdir,
    "Figure6B_ARG_taxon_PMA_vs_nonPMA.png"
  ),
  p,
  width = 7.5,
  height = 6.0,
  dpi = 500
)

# ============================================================
# Save underlying comparative table
# ============================================================

comparison <- plot_df %>%
  select(
    ARG,
    Taxon,
    STATUS,
    Samples,
    Prevalence
  ) %>%
  arrange(
    ARG,
    Taxon,
    STATUS
  )

write.table(
  comparison,
  "analysis_ready/arg_host/ARG_taxon_PMA_vs_nonPMA_selected.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# Print comparison
# ============================================================

cat("\n===== FIGURE 6B: PMA vs non-PMA ARG-TAXON PREVALENCE =====\n\n")

wide <- comparison %>%
  mutate(
    Value = paste0(
      Samples,
      "/32 (",
      round(100 * Prevalence, 1),
      "%)"
    )
  ) %>%
  select(
    ARG,
    Taxon,
    STATUS,
    Value
  ) %>%
  pivot_wider(
    names_from = STATUS,
    values_from = Value
  ) %>%
  arrange(ARG, Taxon)

print(wide, n = Inf)

cat(
  "\nSelection criterion:",
  "association detected in >=2 independent PMA samples.\n"
)

cat(
  "Cell values = samples detected / 32 and prevalence percentage.\n"
)

cat(
  "Important: prevalence differences are descriptive and",
  "may be influenced by substantially different assembly depth.\n"
)

cat("\nSaved revised Figure 6B.\n")
