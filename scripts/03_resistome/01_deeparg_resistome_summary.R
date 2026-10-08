
# ============================================================
# DeepARG resistome summary
# NASA MT2 PMA vs non-PMA
#
# Primary quantitative abundance:
#   DeepARG 16S-normalized abundance
#
# Design:
#   32 matched PMA/non-PMA pairs
# ============================================================

.libPaths(c(path.expand("~/R/4.5"), .libPaths()))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

OUTDIR <- "analysis_ready/deeparg/figure4"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Read data
# ------------------------------------------------------------

sub <- read.delim(
  "analysis_ready/deeparg/deeparg_subtype_long.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

catg <- read.delim(
  "analysis_ready/deeparg/deeparg_category_long.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

master <- read.delim(
  "analysis_ready/metadata/master_sample_table.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n===== INPUT =====\n")
cat("Subtype rows:", nrow(sub), "\n")
cat("Category rows:", nrow(catg), "\n")
cat("Runs:", length(unique(sub$Run)), "\n")
cat("Pairs:", length(unique(sub$Pair_ID)), "\n")

# ------------------------------------------------------------
# Sample-level subtype metrics
# ------------------------------------------------------------

sample_sub <- sub %>%
  group_by(
    Run,
    STATUS,
    SampleName,
    Pair_ID
  ) %>%
  summarise(
    DeepARG_total_16S =
      sum(`16s-NormalizedReadCount`, na.rm = TRUE),

    DeepARG_total_reads =
      sum(ReadCount, na.rm = TRUE),

    ARG_subtype_richness =
      n_distinct(
        `ARG-group`[
          ReadCount > 0 &
          !is.na(ReadCount)
        ]
      ),

    .groups = "drop"
  )

# ------------------------------------------------------------
# Category richness
# ------------------------------------------------------------

sample_cat <- catg %>%
  group_by(Run) %>%
  summarise(
    ARG_category_richness =
      n_distinct(
        `ARG-category`[
          ReadCount > 0 &
          !is.na(ReadCount)
        ]
      ),
    .groups = "drop"
  )

sample <- sample_sub %>%
  left_join(sample_cat, by = "Run") %>%
  left_join(
    master %>%
      select(
        Run,
        HostRemoved_read_pairs,
        Assembly_total_bp
      ),
    by = "Run"
  )

sample$STATUS <- factor(
  sample$STATUS,
  levels = c("non-PMA", "PMA")
)

write.table(
  sample,
  file.path(
    OUTDIR,
    "DeepARG_sample_resistome_metrics.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Group summaries
# ------------------------------------------------------------

summary_tab <- sample %>%
  group_by(STATUS) %>%
  summarise(
    n = n(),

    median_total16S =
      median(DeepARG_total_16S, na.rm = TRUE),

    IQR_total16S =
      IQR(DeepARG_total_16S, na.rm = TRUE),

    median_subtype_richness =
      median(ARG_subtype_richness, na.rm = TRUE),

    IQR_subtype_richness =
      IQR(ARG_subtype_richness, na.rm = TRUE),

    median_category_richness =
      median(ARG_category_richness, na.rm = TRUE),

    IQR_category_richness =
      IQR(ARG_category_richness, na.rm = TRUE),

    .groups = "drop"
  )

write.table(
  summary_tab,
  file.path(
    OUTDIR,
    "DeepARG_group_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n===== GROUP SUMMARY =====\n")
print(summary_tab)

# ------------------------------------------------------------
# Paired helper
# ------------------------------------------------------------

paired_test <- function(data, variable) {

  w <- data %>%
    select(
      Pair_ID,
      STATUS,
      all_of(variable)
    ) %>%
    pivot_wider(
      names_from = STATUS,
      values_from = all_of(variable)
    ) %>%
    filter(
      !is.na(PMA),
      !is.na(`non-PMA`)
    ) %>%
    mutate(
      delta = PMA - `non-PMA`
    )

  tst <- wilcox.test(
    w$PMA,
    w$`non-PMA`,
    paired = TRUE,
    exact = FALSE
  )

  data.frame(
    Metric = variable,
    N_pairs = nrow(w),

    Median_nonPMA =
      median(w$`non-PMA`, na.rm = TRUE),

    Median_PMA =
      median(w$PMA, na.rm = TRUE),

    Median_delta =
      median(w$delta, na.rm = TRUE),

    PMA_higher =
      sum(w$delta > 0, na.rm = TRUE),

    PMA_lower =
      sum(w$delta < 0, na.rm = TRUE),

    Equal =
      sum(w$delta == 0, na.rm = TRUE),

    Wilcoxon_V =
      unname(tst$statistic),

    P_value =
      tst$p.value
  )
}

tests <- bind_rows(
  paired_test(
    sample,
    "DeepARG_total_16S"
  ),
  paired_test(
    sample,
    "ARG_subtype_richness"
  ),
  paired_test(
    sample,
    "ARG_category_richness"
  )
)

tests$FDR_BH <- p.adjust(
  tests$P_value,
  method = "BH"
)

write.table(
  tests,
  file.path(
    OUTDIR,
    "DeepARG_paired_tests.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n===== PAIRED TESTS =====\n")
print(tests)

# ------------------------------------------------------------
# Sequencing-depth correlations
# ------------------------------------------------------------

depth_metrics <- c(
  "DeepARG_total_16S",
  "ARG_subtype_richness",
  "ARG_category_richness"
)

depth_cor <- lapply(
  depth_metrics,
  function(v) {

    ok <- complete.cases(
      sample[[v]],
      sample$HostRemoved_read_pairs
    )

    z <- suppressWarnings(
      cor.test(
        sample[[v]][ok],
        sample$HostRemoved_read_pairs[ok],
        method = "spearman",
        exact = FALSE
      )
    )

    data.frame(
      Metric = v,
      Spearman_rho =
        unname(z$estimate),
      P_value =
        z$p.value
    )
  }
) %>%
  bind_rows()

depth_cor$FDR_BH <- p.adjust(
  depth_cor$P_value,
  method = "BH"
)

write.table(
  depth_cor,
  file.path(
    OUTDIR,
    "DeepARG_depth_correlations.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n===== DEPTH CORRELATIONS =====\n")
print(depth_cor)

# ------------------------------------------------------------
# Pair-level abundance ratios
# ------------------------------------------------------------

pair_abundance <- sample %>%
  select(
    Pair_ID,
    STATUS,
    DeepARG_total_16S
  ) %>%
  pivot_wider(
    names_from = STATUS,
    values_from = DeepARG_total_16S
  ) %>%
  mutate(
    PMA_to_nonPMA_ratio =
      ifelse(
        `non-PMA` > 0,
        PMA / `non-PMA`,
        NA_real_
      ),

    log2_ratio =
      log2(
        (PMA + 1e-8) /
        (`non-PMA` + 1e-8)
      )
  )

write.table(
  pair_abundance,
  file.path(
    OUTDIR,
    "DeepARG_pair_abundance_ratios.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n===== PAIR RATIO SUMMARY =====\n")

print(
  summary(
    pair_abundance$PMA_to_nonPMA_ratio
  )
)

cat("\n===== COMPLETE =====\n")
cat("Outputs:", OUTDIR, "\n")

