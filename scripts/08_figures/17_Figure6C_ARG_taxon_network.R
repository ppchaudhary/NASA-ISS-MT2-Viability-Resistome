library(tidyverse)
library(ggplot2)

infile <- "analysis_ready/arg_host/ARG_taxon_recurrent_PMA.tsv"
outdir <- "figures/figure6"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

x <- read.delim(infile, check.names = FALSE) %>%
  rename(
    ARG = `Element symbol`,
    Taxon = Taxon_name
  ) %>%
  filter(Samples >= 2)

# ------------------------------------------------------------
# Node support
# Sum recurrence across recurrent edges for visual sizing only.
# ------------------------------------------------------------

arg_nodes <- x %>%
  group_by(ARG) %>%
  summarise(
    Support = sum(Samples),
    Max_support = max(Samples),
    .groups = "drop"
  ) %>%
  arrange(desc(Support), desc(Max_support), ARG)

tax_nodes <- x %>%
  group_by(Taxon) %>%
  summarise(
    Support = sum(Samples),
    Max_support = max(Samples),
    .groups = "drop"
  ) %>%
  arrange(desc(Support), desc(Max_support), Taxon)

# ------------------------------------------------------------
# Coordinates
# ARGs left; taxa right
# ------------------------------------------------------------

arg_nodes <- arg_nodes %>%
  mutate(
    x = 0,
    y = seq(
      from = n(),
      to = 1,
      length.out = n()
    ),
    Type = "ARG"
  )

tax_nodes <- tax_nodes %>%
  mutate(
    x = 1,
    y = seq(
      from = n(),
      to = 1,
      length.out = n()
    ),
    Type = "Taxonomic assignment"
  )

# ------------------------------------------------------------
# Edge coordinates
# ------------------------------------------------------------

edges <- x %>%
  left_join(
    arg_nodes %>% select(ARG, x_arg = x, y_arg = y),
    by = "ARG"
  ) %>%
  left_join(
    tax_nodes %>% select(Taxon, x_tax = x, y_tax = y),
    by = "Taxon"
  )

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p <- ggplot() +

  # recurrent ARG-taxon edges
  geom_curve(
    data = edges,
    aes(
      x = x_arg,
      y = y_arg,
      xend = x_tax,
      yend = y_tax,
      linewidth = Samples
    ),
    curvature = 0.08,
    color = "grey45",
    alpha = 0.75,
    lineend = "round"
  ) +

  # ARG nodes
  geom_point(
    data = arg_nodes,
    aes(
      x = x,
      y = y,
      size = Support
    ),
    shape = 21,
    fill = "#FF4B4B",
    color = "black",
    stroke = 0.5
  ) +

  # taxon nodes
  geom_point(
    data = tax_nodes,
    aes(
      x = x,
      y = y,
      size = Support
    ),
    shape = 21,
    fill = "#1479E8",
    color = "black",
    stroke = 0.5
  ) +

  # ARG labels — italic
  geom_text(
    data = arg_nodes,
    aes(
      x = x - 0.045,
      y = y,
      label = ARG
    ),
    hjust = 1,
    size = 4.2,
    fontface = "italic",
    color = "black"
  ) +

  # Taxon labels — regular because assignments are mixed ranks
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

  scale_linewidth_continuous(
    name = "PMA samples",
    breaks = 2:5,
    range = c(0.8, 3.2)
  ) +

  scale_size_continuous(
    name = "Total recurrent\nsample support",
    range = c(4.5, 8.5)
  ) +

  coord_cartesian(
    xlim = c(-0.35, 1.42),
    clip = "off"
  ) +

  labs(
    x = NULL,
    y = NULL
  ) +

  theme_void(base_size = 12) +

  theme(
    legend.position = "right",

    legend.title = element_text(
      face = "bold",
      size = 10
    ),

    legend.text = element_text(
      size = 9
    ),

    plot.margin = margin(
      15, 30, 15, 30
    )
  )

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

ggsave(
  file.path(outdir, "Figure6C_ARG_taxon_network.pdf"),
  p,
  width = 8.2,
  height = 5.6,
  device = cairo_pdf
)

ggsave(
  file.path(outdir, "Figure6C_ARG_taxon_network.png"),
  p,
  width = 8.2,
  height = 5.6,
  dpi = 500
)

# ------------------------------------------------------------
# Print exact network contents
# ------------------------------------------------------------

cat("\n===== FIGURE 6C RECURRENT PMA NETWORK =====\n\n")

cat("Edges:\n")
print(
  edges %>%
    select(ARG, Taxon, Samples, Contigs) %>%
    arrange(desc(Samples), ARG, Taxon)
)

cat("\nARG nodes:\n")
print(
  arg_nodes %>%
    select(ARG, Support)
)

cat("\nTaxonomic-assignment nodes:\n")
print(
  tax_nodes %>%
    select(Taxon, Support)
)

cat("\nEdges:", nrow(edges), "\n")
cat("ARG nodes:", nrow(arg_nodes), "\n")
cat("Taxonomic-assignment nodes:", nrow(tax_nodes), "\n")

cat("\nSaved Figure 6C.\n")
