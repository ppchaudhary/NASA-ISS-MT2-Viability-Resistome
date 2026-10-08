library(tidyverse)
library(ggplot2)

infile <- "analysis_ready/arg_host/ARG_taxon_recurrent_PMA.tsv"
outdir <- "figures/figure6"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load recurrent PMA ARG-taxon associations
# Same criterion as Figure 6B: >=2 independent PMA samples
# ============================================================

x <- read.delim(infile, check.names = FALSE) %>%
  rename(
    ARG = `Element symbol`,
    Taxon = Taxon_name
  ) %>%
  filter(Samples >= 2)

# ============================================================
# Order nodes to minimize crossings and emphasize structure
# ============================================================

# Taxa ordered by total recurrent sample support
tax_order <- x %>%
  group_by(Taxon) %>%
  summarise(
    Support = sum(Samples),
    Max_support = max(Samples),
    .groups = "drop"
  ) %>%
  arrange(desc(Support), desc(Max_support), Taxon)

# For each ARG, determine its dominant taxonomic assignment
arg_order <- x %>%
  group_by(ARG) %>%
  mutate(
    ARG_total = sum(Samples)
  ) %>%
  arrange(ARG, desc(Samples), Taxon) %>%
  slice(1) %>%
  ungroup() %>%
  left_join(
    tax_order %>%
      mutate(Taxon_order = row_number()) %>%
      select(Taxon, Taxon_order),
    by = "Taxon"
  ) %>%
  arrange(Taxon_order, desc(ARG_total), desc(Samples), ARG)

# ============================================================
# Coordinates
# ============================================================

# ARG nodes on left
arg_nodes <- arg_order %>%
  transmute(
    ARG,
    x = 0,
    y = rev(seq_len(n()))
  )

# Taxon nodes on right
tax_nodes <- tax_order %>%
  mutate(
    x = 1,
    y = rev(seq_len(n()))
  ) %>%
  select(Taxon, Support, x, y)

# Rescale the 4 taxon positions so they span the ARG panel height
if (nrow(tax_nodes) > 1) {
  tax_nodes$y <- scales::rescale(
    tax_nodes$y,
    to = c(1.3, nrow(arg_nodes) - 0.3)
  )
} else {
  tax_nodes$y <- mean(arg_nodes$y)
}

# ============================================================
# Build edges
# ============================================================

edges <- x %>%
  left_join(
    arg_nodes %>% select(ARG, x_arg = x, y_arg = y),
    by = "ARG"
  ) %>%
  left_join(
    tax_nodes %>% select(Taxon, x_tax = x, y_tax = y),
    by = "Taxon"
  )

# ============================================================
# Plot
# ============================================================

p <- ggplot() +

  # ----------------------------------------------------------
  # Flow connections
  # Width = number of independent PMA samples
  # ----------------------------------------------------------

  geom_curve(
    data = edges,
    aes(
      x = x_arg + 0.025,
      y = y_arg,
      xend = x_tax - 0.025,
      yend = y_tax,
      linewidth = Samples
    ),
    curvature = 0.18,
    color = "#FF4B4B",
    alpha = 0.48,
    lineend = "round"
  ) +

  # ----------------------------------------------------------
  # Left ARG anchor points
  # ----------------------------------------------------------

  geom_point(
    data = arg_nodes,
    aes(x = x, y = y),
    shape = 21,
    size = 4.3,
    stroke = 0.55,
    fill = "#FF4B4B",
    color = "black"
  ) +

  # ----------------------------------------------------------
  # Right taxonomic-assignment anchor points
  # ----------------------------------------------------------

  geom_point(
    data = tax_nodes,
    aes(x = x, y = y),
    shape = 21,
    size = 4.3,
    stroke = 0.55,
    fill = "#1479E8",
    color = "black"
  ) +

  # ----------------------------------------------------------
  # ARG labels: italic
  # ----------------------------------------------------------

  geom_text(
    data = arg_nodes,
    aes(
      x = x - 0.045,
      y = y,
      label = ARG
    ),
    hjust = 1,
    size = 4.1,
    fontface = "italic",
    color = "black"
  ) +

  # ----------------------------------------------------------
  # Taxonomic assignment labels: regular
  # Mixed taxonomic ranks, therefore do not call these genera.
  # ----------------------------------------------------------

  geom_text(
    data = tax_nodes,
    aes(
      x = x + 0.045,
      y = y,
      label = Taxon
    ),
    hjust = 0,
    size = 4.0,
    color = "black"
  ) +

  # ----------------------------------------------------------
  # Column headings
  # ----------------------------------------------------------

  annotate(
    "text",
    x = 0,
    y = nrow(arg_nodes) + 0.75,
    label = "ARG",
    fontface = "bold",
    size = 4.2
  ) +

  annotate(
    "text",
    x = 1,
    y = nrow(arg_nodes) + 0.75,
    label = "Taxonomic assignment",
    fontface = "bold",
    size = 4.2
  ) +

  # ----------------------------------------------------------
  # Edge-width legend
  # ----------------------------------------------------------

  scale_linewidth_continuous(
    name = "PMA samples",
    breaks = sort(unique(edges$Samples)),
    range = c(1.2, 5.2)
  ) +

  coord_cartesian(
    xlim = c(-0.30, 1.42),
    ylim = c(0.35, nrow(arg_nodes) + 1.05),
    clip = "off"
  ) +

  guides(
    linewidth = guide_legend(
      title.position = "top",
      override.aes = list(
        color = "grey45",
        alpha = 0.8
      )
    )
  ) +

  theme_void(base_size = 12) +

  theme(
    legend.position = "bottom",

    legend.direction = "horizontal",

    legend.title = element_text(
      face = "bold",
      size = 10
    ),

    legend.text = element_text(
      size = 9
    ),

    plot.margin = margin(
      15, 35, 10, 35
    )
  )

# ============================================================
# Save
# ============================================================

ggsave(
  file.path(outdir, "Figure6C_ARG_taxon_flow.pdf"),
  p,
  width = 8.2,
  height = 5.6,
  device = cairo_pdf
)

ggsave(
  file.path(outdir, "Figure6C_ARG_taxon_flow.png"),
  p,
  width = 8.2,
  height = 5.6,
  dpi = 500
)

# ============================================================
# Report exactly what was plotted
# ============================================================

cat("\n===== FIGURE 6C: RECURRENT PMA ARG-TAXON FLOW =====\n\n")

cat("Criterion: ARG-taxon association observed in >=2 independent PMA samples\n\n")

cat("Associations:\n")
print(
  edges %>%
    select(ARG, Taxon, Samples, Contigs) %>%
    arrange(desc(Samples), ARG, Taxon)
)

cat("\nNetwork dimensions:\n")
cat("Associations:", nrow(edges), "\n")
cat("ARGs:", n_distinct(edges$ARG), "\n")
cat("Taxonomic assignments:", n_distinct(edges$Taxon), "\n")

cat("\nEdge width = number of independent PMA samples.\n")
cat("ARG labels are italic; taxonomic assignments are regular.\n")
cat("\nSaved Figure 6C flow diagram.\n")
