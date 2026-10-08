library(tidyverse)
library(vegan)

outdir <- "analysis_ready/integrated_analysis"
figdir <- "figures/integrated_analysis"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# INPUTS
# ============================================================

meta_file <- "analysis_ready/metadata/master_sample_table.tsv"

tax_file <- paste0(
  "analysis_ready/taxonomy/",
  "bracken_species_relative_abundance.tsv"
)

arg_file <- paste0(
  "analysis_ready/deeparg/",
  "deeparg_subtype_16S_normalized_matrix.tsv"
)

meta <- read.delim(
  meta_file,
  check.names = FALSE
)

cat("\n===== METADATA =====\n")
cat("Samples:", nrow(meta), "\n")
cat("Pairs:", n_distinct(meta$Pair_ID), "\n")
print(table(meta$STATUS))

# ============================================================
# HELPER
# Convert feature x sample matrix to sample x feature matrix
# ============================================================

read_feature_matrix <- function(file, meta) {

  x <- read.delim(
    file,
    check.names = FALSE
  )

  cat("\nReading:", file, "\n")
  cat("Raw dimensions:", dim(x), "\n")

  # Identify columns matching Run IDs
  sample_cols <- intersect(
    colnames(x),
    meta$Run
  )

  cat(
    "Matched sample columns:",
    length(sample_cols),
    "\n"
  )

  if (length(sample_cols) != 64) {
    stop(
      paste(
        "Expected 64 Run columns but found",
        length(sample_cols)
      )
    )
  }

  mat <- x %>%
    select(all_of(sample_cols)) %>%
    as.matrix()

  storage.mode(mat) <- "numeric"

  # transpose: samples x features
  mat <- t(mat)

  # retain Run IDs
  rownames(mat) <- sample_cols

  mat
}

# ============================================================
# MICROBIOME MATRIX
# Bracken species relative abundance
# ============================================================

tax_mat <- read_feature_matrix(
  tax_file,
  meta
)

cat(
  "\nTaxonomy sample x feature:",
  dim(tax_mat),
  "\n"
)

# Remove zero-sum features
tax_mat <- tax_mat[
  ,
  colSums(tax_mat, na.rm = TRUE) > 0,
  drop = FALSE
]

# ============================================================
# RESISTOME MATRIX
# DeepARG subtype 16S-normalized abundance
# ============================================================

arg_mat <- read_feature_matrix(
  arg_file,
  meta
)

cat(
  "DeepARG sample x feature:",
  dim(arg_mat),
  "\n"
)

arg_mat <- arg_mat[
  ,
  colSums(arg_mat, na.rm = TRUE) > 0,
  drop = FALSE
]

# ============================================================
# Ensure identical sample order
# ============================================================

common_runs <- meta$Run[
  meta$Run %in% rownames(tax_mat) &
  meta$Run %in% rownames(arg_mat)
]

if (length(common_runs) != 64) {
  stop(
    paste(
      "Expected 64 common samples but found",
      length(common_runs)
    )
  )
}

tax_mat <- tax_mat[common_runs, , drop = FALSE]
arg_mat <- arg_mat[common_runs, , drop = FALSE]

meta2 <- meta %>%
  filter(Run %in% common_runs) %>%
  arrange(match(Run, common_runs))

stopifnot(
  identical(meta2$Run, rownames(tax_mat)),
  identical(meta2$Run, rownames(arg_mat))
)

# ============================================================
# Calculate WITHIN-PAIR Bray-Curtis distance
#
# One value per Pair_ID:
# distance between PMA and its matched non-PMA sample.
# ============================================================

pair_ids <- unique(meta2$Pair_ID)

results <- vector(
  "list",
  length(pair_ids)
)

for (i in seq_along(pair_ids)) {

  pid <- pair_ids[i]

  m <- meta2 %>%
    filter(Pair_ID == pid)

  if (
    nrow(m) != 2 ||
    !all(c("PMA", "non-PMA") %in% m$STATUS)
  ) {
    stop(
      paste(
        "Incomplete pair:",
        pid
      )
    )
  }

  run_non <- m$Run[
    m$STATUS == "non-PMA"
  ]

  run_pma <- m$Run[
    m$STATUS == "PMA"
  ]

  # --------------------------------------------
  # Microbiome Bray-Curtis
  # --------------------------------------------

  tax_pair <- tax_mat[
    c(run_non, run_pma),
    ,
    drop = FALSE
  ]

  micro_bc <- as.numeric(
    vegdist(
      tax_pair,
      method = "bray"
    )
  )

  # --------------------------------------------
  # Resistome Bray-Curtis
  # --------------------------------------------

  arg_pair <- arg_mat[
    c(run_non, run_pma),
    ,
    drop = FALSE
  ]

  arg_bc <- as.numeric(
    vegdist(
      arg_pair,
      method = "bray"
    )
  )

  results[[i]] <- tibble(
    Pair_ID = pid,
    Run_nonPMA = run_non,
    Run_PMA = run_pma,
    Microbiome_BrayCurtis = micro_bc,
    Resistome_BrayCurtis = arg_bc
  )
}

pair_dist <- bind_rows(results)

# ============================================================
# Add useful pair metadata
# ============================================================

pair_meta <- meta2 %>%
  group_by(Pair_ID) %>%
  summarise(
    Flight = first(Flight),
    LocationCode = first(LocationCode),
    Location = first(Location),
    .groups = "drop"
  )

pair_dist <- pair_dist %>%
  left_join(
    pair_meta,
    by = "Pair_ID"
  )

# ============================================================
# Spearman association
# ============================================================

cor_test <- cor.test(
  pair_dist$Microbiome_BrayCurtis,
  pair_dist$Resistome_BrayCurtis,
  method = "spearman",
  exact = FALSE
)

rho <- unname(
  cor_test$estimate
)

pval <- cor_test$p.value

stats <- tibble(
  N_pairs = nrow(pair_dist),
  Spearman_rho = rho,
  P_value = pval
)

# ============================================================
# Save tables
# ============================================================

write.table(
  pair_dist,
  file.path(
    outdir,
    "paired_microbiome_resistome_distances.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  stats,
  file.path(
    outdir,
    "microbiome_resistome_coupling_stats.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# Exploratory plot
# NOT yet a manuscript figure
# ============================================================

label_txt <- paste0(
  "Spearman \u03c1 = ",
  sprintf("%.2f", rho),
  "\nP = ",
  format.pval(
    pval,
    digits = 2,
    eps = 0.001
  )
)

p <- ggplot(
  pair_dist,
  aes(
    x = Microbiome_BrayCurtis,
    y = Resistome_BrayCurtis
  )
) +

  geom_point(
    size = 3,
    alpha = 0.8,
    color = "#1479E8"
  ) +

  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    linewidth = 0.7
  ) +

  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = label_txt,
    hjust = 1.1,
    vjust = 1.2,
    size = 4.2
  ) +

  labs(
    x = "Within-pair microbiome Bray\u2013Curtis distance",
    y = "Within-pair resistome Bray\u2013Curtis distance"
  ) +

  theme_classic(base_size = 12) +

  theme(
    axis.title = element_text(
      color = "black"
    ),

    axis.text = element_text(
      color = "black"
    )
  )

ggsave(
  file.path(
    figdir,
    "microbiome_resistome_coupling_exploratory.pdf"
  ),
  p,
  width = 6.2,
  height = 5.2,
  device = cairo_pdf
)

ggsave(
  file.path(
    figdir,
    "microbiome_resistome_coupling_exploratory.png"
  ),
  p,
  width = 6.2,
  height = 5.2,
  dpi = 400
)

# ============================================================
# REPORT
# ============================================================

cat("\n============================================\n")
cat("MICROBIOME-RESISTOME COUPLING\n")
cat("============================================\n\n")

cat(
  "Complete pairs:",
  nrow(pair_dist),
  "\n\n"
)

cat("Microbiome Bray-Curtis:\n")
print(
  summary(
    pair_dist$Microbiome_BrayCurtis
  )
)

cat("\nResistome Bray-Curtis:\n")
print(
  summary(
    pair_dist$Resistome_BrayCurtis
  )
)

cat("\nAssociation:\n")
cat(
  "Spearman rho =",
  rho,
  "\n"
)

cat(
  "P =",
  pval,
  "\n"
)

cat("\nSaved pair-level distances and statistics.\n")
cat("Exploratory plot saved; not yet designated as a manuscript figure.\n")
cat("============================================\n")
