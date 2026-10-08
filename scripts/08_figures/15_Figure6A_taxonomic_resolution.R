library(tidyverse)

infile <- "analysis_ready/amrfinder/amrfinder_ARG_contig_taxonomy.tsv"
outdir <- "figures/figure6"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

x <- read.delim(infile, check.names = FALSE)

# ------------------------------------------------------------
# Rank labels
# ------------------------------------------------------------
rank_labels <- c(
  D  = "Domain",
  D1 = "Domain/subrank",
  P  = "Phylum",
  C  = "Class",
  O  = "Order",
  F  = "Family",
  G  = "Genus",
  G1 = "Genus/subrank",
  S  = "Species",
  S1 = "Species/subrank",
  U  = "Unclassified"
)

# Count unique ARG-bearing contigs within each sample/status.
# A contig carrying multiple AMRFinder hits should not be counted
# repeatedly for taxonomic-resolution composition.
d <- x %>%
  distinct(Run, STATUS, `Contig id`, Taxon_rank) %>%
  mutate(
    Rank = recode(Taxon_rank, !!!rank_labels),
    Rank = ifelse(is.na(Rank), Taxon_rank, Rank)
  ) %>%
  count(STATUS, Rank, name = "ARG_contigs") %>%
  group_by(STATUS) %>%
  mutate(
    Total = sum(ARG_contigs),
    Proportion = ARG_contigs / Total
  ) %>%
  ungroup()

# Order from broad to specific
rank_order <- c(
  "Domain",
  "Domain/subrank",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Genus/subrank",
  "Species",
  "Species/subrank",
  "Unclassified"
)

d$Rank <- factor(d$Rank, levels = rank_order)

# Keep desired group order
d$STATUS <- factor(d$STATUS, levels = c("non-PMA", "PMA"))

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
p <- ggplot(d, aes(x = STATUS, y = Proportion, fill = Rank)) +
  geom_col(
    width = 0.66,
    color = "white",
    linewidth = 0.35
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Proportion of classified ARG-bearing contigs",
    fill = "Taxonomic rank"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      color = "black",
      face = "bold",
      size = 11
    ),
    axis.text.y = element_text(
      color = "black",
      size = 10
    ),
    axis.title.y = element_text(
      color = "black",
      size = 11
    ),
    legend.title = element_text(
      face = "bold",
      size = 10
    ),
    legend.text = element_text(
      size = 9
    ),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  file.path(outdir, "Figure6A_taxonomic_resolution.pdf"),
  p,
  width = 6.5,
  height = 4.8,
  device = cairo_pdf
)

ggsave(
  file.path(outdir, "Figure6A_taxonomic_resolution.png"),
  p,
  width = 6.5,
  height = 4.8,
  dpi = 500
)

# ------------------------------------------------------------
# Print summary
# ------------------------------------------------------------
cat("\n===== FIGURE 6A TAXONOMIC RESOLUTION =====\n")

print(
  d %>%
    mutate(Percent = round(100 * Proportion, 1)) %>%
    arrange(STATUS, Rank)
)

cat("\nSaved Figure 6A.\n")
