library(tidyverse)
library(ggplot2)

infile <- "analysis_ready/amrfinder/amrfinder_ARG_contig_taxonomy.tsv"
outdir <- "figures/figure6"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

x <- read.delim(infile, check.names = FALSE)

# ============================================================
# One taxonomic assignment per unique sample-contig
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
    `Contig id`,
    Taxon_name
  )

# ============================================================
# Keep biologically informative assignments
# Everything else -> Other assignments
#
# Mixed classifier-reported ranks are intentionally preserved.
# ============================================================

keep_taxa <- c(
  "Staphylococcus",
  "Staphylococcus epidermidis",
  "Actinomycetes",
  "Staphylococcaceae",
  "Corynebacterium",
  "Streptococcus",
  "Bacteria"
)

d <- d %>%
  mutate(
    Taxon_plot = if_else(
      Taxon_name %in% keep_taxa,
      Taxon_name,
      "Other assignments"
    )
  )

# ============================================================
# Composition based on unique ARG-bearing sample-contigs
# ============================================================

comp <- d %>%
  count(
    STATUS,
    Taxon_plot,
    name = "ARG_contigs"
  ) %>%
  group_by(STATUS) %>%
  mutate(
    Total_ARG_contigs = sum(ARG_contigs),
    Proportion = ARG_contigs / Total_ARG_contigs
  ) %>%
  ungroup()

# Order legend biologically
taxon_levels <- c(
  "Staphylococcus",
  "Staphylococcus epidermidis",
  "Actinomycetes",
  "Staphylococcaceae",
  "Corynebacterium",
  "Streptococcus",
  "Bacteria",
  "Other assignments"
)

comp <- comp %>%
  mutate(
    Taxon_plot = factor(
      Taxon_plot,
      levels = taxon_levels
    ),

    # rev() makes non-PMA appear above PMA
    STATUS = factor(
      STATUS,
      levels = c("PMA", "non-PMA")
    )
  )

# ============================================================
# Manuscript palette
# ============================================================

tax_colors <- c(
  "Staphylococcus" = "#F8766D",
  "Staphylococcus epidermidis" = "#E78600",
  "Actinomycetes" = "#B79F00",
  "Staphylococcaceae" = "#53B400",
  "Corynebacterium" = "#00BA38",
  "Streptococcus" = "#00C094",
  "Bacteria" = "#00BFC4",
  "Other assignments" = "grey75"
)

# ============================================================
# Plot
# ============================================================

p <- ggplot(
  comp,
  aes(
    x = Proportion,
    y = STATUS,
    fill = Taxon_plot
  )
) +

  geom_col(
    width = 0.48,
    color = "white",
    linewidth = 0.25
  ) +

  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +

  scale_fill_manual(
    values = tax_colors,
    drop = FALSE,
    name = "Taxonomic assignment"
  ) +

  labs(
    x = "Proportion of classified ARG-bearing contigs",
    y = NULL
  ) +

  theme_classic(base_size = 12) +

  theme(
    axis.text.x = element_text(
      color = "black",
      size = 10
    ),

    axis.text.y = element_text(
      color = "black",
      face = "bold",
      size = 11
    ),

    axis.title.x = element_text(
      size = 11
    ),

    legend.position = "right",

    legend.title = element_text(
      face = "bold",
      size = 10
    ),

    legend.text = element_text(
      size = 9
    ),

    legend.key.height = grid::unit(0.45, "cm"),

    plot.margin = margin(
      10, 10, 10, 10
    )
  )

# ============================================================
# Save
# ============================================================

ggsave(
  file.path(
    outdir,
    "Figure6A_ARG_taxon_composition.pdf"
  ),
  p,
  width = 7.2,
  height = 5.2,
  device = cairo_pdf
)

ggsave(
  file.path(
    outdir,
    "Figure6A_ARG_taxon_composition.png"
  ),
  p,
  width = 7.2,
  height = 5.2,
  dpi = 500
)

# ============================================================
# Save exact underlying values
# ============================================================

write.table(
  comp %>%
    arrange(STATUS, desc(Proportion)),
  "analysis_ready/arg_host/Figure6A_ARG_taxon_composition.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n===== FIGURE 6A TAXONOMIC COMPOSITION =====\n\n")

print(
  comp %>%
    select(
      STATUS,
      Taxon_plot,
      ARG_contigs,
      Total_ARG_contigs,
      Proportion
    ) %>%
    arrange(STATUS, desc(Proportion)),
  n = Inf
)

cat("\nSaved revised Figure 6A.\n")
