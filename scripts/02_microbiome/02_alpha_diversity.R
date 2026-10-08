#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ============================================================
# Paths
# ============================================================

ROOT <- getwd()

count_file <- file.path(
  ROOT,
  "analysis_ready/taxonomy/bracken_species_counts.tsv"
)

master_file <- file.path(
  ROOT,
  "analysis_ready/metadata/master_sample_table.tsv"
)

outdir <- file.path(
  ROOT,
  "analysis_ready/taxonomy/alpha_diversity"
)

figdir <- file.path(
  ROOT,
  "figures/figure2"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Figure colors
# ============================================================

COL_NONPMA <- "#1479E8"
COL_PMA    <- "#FF4B4B"
COL_PAIR   <- "grey80"

# ============================================================
# Load Bracken counts
# ============================================================

raw <- read.delim(
  count_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("========================================\n")
cat("ALPHA DIVERSITY ANALYSIS\n")
cat("========================================\n")

cat("Raw Bracken dimensions:",
    nrow(raw), "x", ncol(raw), "\n")

required_tax_cols <- c("taxonomy_id", "name")

if (!all(required_tax_cols %in% colnames(raw))) {
  stop("Expected taxonomy_id and name columns not found.")
}

sample_cols <- setdiff(
  colnames(raw),
  required_tax_cols
)

# species x samples
count_species_sample <- as.matrix(
  raw[, sample_cols, drop = FALSE]
)

storage.mode(count_species_sample) <- "numeric"

# samples x species
counts <- t(count_species_sample)

cat("Samples:", nrow(counts), "\n")
cat("Species:", ncol(counts), "\n")

# ============================================================
# Metadata
# ============================================================

meta <- read.delim(
  master_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(rownames(counts) %in% meta$Run)) {
  stop("Some Bracken sample IDs are missing from metadata.")
}

meta2 <- meta[
  match(rownames(counts), meta$Run),
]

if (!all(meta2$Run == rownames(counts))) {
  stop("Sample order mismatch.")
}

cat("Metadata matched:", nrow(meta2), "\n")
cat("Unique pairs:", length(unique(meta2$Pair_ID)), "\n")

# ============================================================
# Alpha diversity
# ============================================================

observed <- specnumber(counts)

shannon <- diversity(
  counts,
  index = "shannon"
)

simpson <- diversity(
  counts,
  index = "simpson"
)

alpha <- data.frame(
  Run = rownames(counts),
  Observed = observed,
  Shannon = shannon,
  Simpson = simpson,
  stringsAsFactors = FALSE
)

alpha <- alpha %>%
  left_join(
    meta2 %>%
      select(
        Run,
        SampleName,
        Pair_ID,
        STATUS,
        Flight,
        LocationCode,
        Location,
        HostRemoved_read_pairs
      ),
    by = "Run"
  )

alpha$STATUS <- factor(
  alpha$STATUS,
  levels = c("non-PMA", "PMA")
)

write.table(
  alpha,
  file.path(outdir, "alpha_diversity_metrics.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ============================================================
# Group summaries
# ============================================================

summary_table <- alpha %>%
  group_by(STATUS) %>%
  summarise(
    n = n(),

    Observed_median = median(Observed),
    Observed_IQR = IQR(Observed),

    Shannon_median = median(Shannon),
    Shannon_IQR = IQR(Shannon),

    Simpson_median = median(Simpson),
    Simpson_IQR = IQR(Simpson),

    .groups = "drop"
  )

cat("\n========================================\n")
cat("GROUP SUMMARIES\n")
cat("========================================\n")

print(summary_table)

write.table(
  summary_table,
  file.path(outdir, "alpha_diversity_group_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ============================================================
# Paired dataset
# ============================================================

paired <- alpha %>%
  select(
    Pair_ID,
    STATUS,
    Observed,
    Shannon,
    Simpson
  ) %>%
  pivot_wider(
    names_from = STATUS,
    values_from = c(
      Observed,
      Shannon,
      Simpson
    )
  )

if (nrow(paired) != 32) {
  stop("Expected exactly 32 matched pairs.")
}

cat("\nMatched pairs:", nrow(paired), "\n")

# ============================================================
# Paired tests
# ============================================================

run_paired_test <- function(metric) {

  nonpma <- paired[[paste0(metric, "_non-PMA")]]
  pma    <- paired[[paste0(metric, "_PMA")]]

  wt <- suppressWarnings(
    wilcox.test(
      pma,
      nonpma,
      paired = TRUE,
      exact = FALSE,
      conf.int = FALSE
    )
  )

  delta <- pma - nonpma

  data.frame(
    Metric = metric,
    N_pairs = length(delta),

    Median_nonPMA = median(nonpma),
    Median_PMA = median(pma),

    Median_paired_difference_PMA_minus_nonPMA =
      median(delta),

    Pairs_PMA_higher = sum(delta > 0),
    Pairs_PMA_lower = sum(delta < 0),
    Pairs_equal = sum(delta == 0),

    Wilcoxon_V = unname(wt$statistic),
    P_value = wt$p.value,

    stringsAsFactors = FALSE
  )
}

test_results <- bind_rows(
  run_paired_test("Observed"),
  run_paired_test("Shannon"),
  run_paired_test("Simpson")
)

test_results$FDR_BH <- p.adjust(
  test_results$P_value,
  method = "BH"
)

cat("\n========================================\n")
cat("PAIRED WILCOXON RESULTS\n")
cat("PMA vs matched non-PMA\n")
cat("========================================\n")

print(test_results)

write.table(
  test_results,
  file.path(outdir, "alpha_diversity_paired_tests.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ============================================================
# Depth correlations
# ============================================================

depth_results <- bind_rows(
  data.frame(
    Metric = "Observed",
    Spearman_rho = cor(
      alpha$Observed,
      alpha$HostRemoved_read_pairs,
      method = "spearman"
    )
  ),

  data.frame(
    Metric = "Shannon",
    Spearman_rho = cor(
      alpha$Shannon,
      alpha$HostRemoved_read_pairs,
      method = "spearman"
    )
  ),

  data.frame(
    Metric = "Simpson",
    Spearman_rho = cor(
      alpha$Simpson,
      alpha$HostRemoved_read_pairs,
      method = "spearman"
    )
  )
)

cat("\n========================================\n")
cat("DEPTH CORRELATIONS\n")
cat("Spearman rho with host-removed read pairs\n")
cat("========================================\n")

print(depth_results)

write.table(
  depth_results,
  file.path(outdir, "alpha_diversity_depth_correlations.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ============================================================
# Long format for plotting
# ============================================================

plotdata <- alpha %>%
  select(
    Run,
    Pair_ID,
    STATUS,
    Observed,
    Shannon,
    Simpson
  ) %>%
  pivot_longer(
    cols = c(
      Observed,
      Shannon,
      Simpson
    ),
    names_to = "Metric",
    values_to = "Value"
  )

# ============================================================
# Plotting function
# ============================================================

make_plot <- function(metric_name, y_label) {

  d <- plotdata %>%
    filter(Metric == metric_name)

  ggplot(
    d,
    aes(
      x = STATUS,
      y = Value
    )
  ) +

    geom_line(
      aes(group = Pair_ID),
      color = COL_PAIR,
      linewidth = 0.45,
      alpha = 0.75
    ) +

    geom_boxplot(
      aes(fill = STATUS, group = STATUS),
      width = 0.46,
      alpha = 0.18,
      outlier.shape = NA,
      linewidth = 0.65
    ) +

    geom_point(
      aes(color = STATUS),
      size = 2.2,
      alpha = 0.90,
      position = position_jitter(
        width = 0.045,
        height = 0
      )
    ) +

    scale_color_manual(
      values = c(
        "non-PMA" = COL_NONPMA,
        "PMA" = COL_PMA
      )
    ) +

    scale_fill_manual(
      values = c(
        "non-PMA" = COL_NONPMA,
        "PMA" = COL_PMA
      )
    ) +

    labs(
      x = NULL,
      y = y_label
    ) +

    theme_classic(base_size = 12) +

    theme(
      legend.position = "none",
      axis.text.x = element_text(
        color = "black",
        size = 11
      ),
      axis.text.y = element_text(
        color = "black"
      ),
      axis.title = element_text(
        color = "black",
        face = "bold"
      )
    )
}

p_shannon <- make_plot(
  "Shannon",
  "Shannon diversity"
)

p_observed <- make_plot(
  "Observed",
  "Observed species richness"
)

p_simpson <- make_plot(
  "Simpson",
  "Simpson diversity"
)

# ============================================================
# Save plots
# ============================================================

ggsave(
  file.path(figdir, "Fig2A_Shannon_paired.pdf"),
  p_shannon,
  width = 4.1,
  height = 4.0
)

ggsave(
  file.path(figdir, "Fig2A_Shannon_paired.png"),
  p_shannon,
  width = 4.1,
  height = 4.0,
  dpi = 600
)

ggsave(
  file.path(figdir, "Fig2B_Observed_paired.pdf"),
  p_observed,
  width = 4.1,
  height = 4.0
)

ggsave(
  file.path(figdir, "Fig2B_Observed_paired.png"),
  p_observed,
  width = 4.1,
  height = 4.0,
  dpi = 600
)

ggsave(
  file.path(figdir, "FigS_alpha_Simpson_paired.pdf"),
  p_simpson,
  width = 4.1,
  height = 4.0
)

ggsave(
  file.path(figdir, "FigS_alpha_Simpson_paired.png"),
  p_simpson,
  width = 4.1,
  height = 4.0,
  dpi = 600
)

cat("\n========================================\n")
cat("FILES SAVED\n")
cat("========================================\n")

cat(
  file.path(outdir, "alpha_diversity_metrics.tsv"),
  "\n"
)

cat(
  file.path(outdir, "alpha_diversity_group_summary.tsv"),
  "\n"
)

cat(
  file.path(outdir, "alpha_diversity_paired_tests.tsv"),
  "\n"
)

cat(
  file.path(outdir, "alpha_diversity_depth_correlations.tsv"),
  "\n"
)

cat(
  file.path(figdir, "Fig2A_Shannon_paired.pdf"),
  "\n"
)

cat(
  file.path(figdir, "Fig2B_Observed_paired.pdf"),
  "\n"
)

cat("\nDONE\n")
