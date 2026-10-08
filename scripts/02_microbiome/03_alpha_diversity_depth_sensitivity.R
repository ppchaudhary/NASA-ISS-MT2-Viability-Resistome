#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vegan)
  library(dplyr)
  library(tidyr)
})

set.seed(20261005)

ROOT <- getwd()

COUNT <- file.path(
  ROOT, "analysis_ready/taxonomy/bracken_species_counts.tsv"
)

MASTER <- file.path(
  ROOT, "analysis_ready/metadata/master_sample_table.tsv"
)

OUTDIR <- file.path(
  ROOT, "analysis_ready/taxonomy/alpha_diversity"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

N_ITER <- 100

# ------------------------------------------------------------
# Load counts
# ------------------------------------------------------------

raw <- read.delim(
  COUNT,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sample_cols <- setdiff(
  colnames(raw),
  c("taxonomy_id", "name")
)

counts <- t(
  as.matrix(raw[, sample_cols, drop = FALSE])
)

storage.mode(counts) <- "numeric"

meta <- read.delim(
  MASTER,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

meta <- meta[
  match(rownames(counts), meta$Run),
]

stopifnot(all(meta$Run == rownames(counts)))
stopifnot(length(unique(meta$Pair_ID)) == 32)

# Bracken counts are estimated counts and should be integer-valued
# for rrarefy. Round defensively.
counts <- round(counts)

# ------------------------------------------------------------
# Pair-specific rarefaction depths
# ------------------------------------------------------------

libsize <- rowSums(counts)

depth_table <- data.frame(
  Run = rownames(counts),
  Pair_ID = meta$Pair_ID,
  STATUS = meta$STATUS,
  Bracken_reads = libsize,
  stringsAsFactors = FALSE
)

pair_depth <- depth_table %>%
  group_by(Pair_ID) %>%
  summarise(
    Rarefaction_depth = min(Bracken_reads),
    .groups = "drop"
  )

depth_table <- depth_table %>%
  left_join(pair_depth, by = "Pair_ID")

write.table(
  depth_table,
  file.path(
    OUTDIR,
    "pair_specific_rarefaction_depths.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("========================================\n")
cat("PAIR-SPECIFIC DEPTH SENSITIVITY\n")
cat("========================================\n")

cat("Samples:", nrow(counts), "\n")
cat("Pairs:", length(unique(meta$Pair_ID)), "\n")
cat("Iterations:", N_ITER, "\n")

cat("\nPair-specific rarefaction depths:\n")
print(summary(pair_depth$Rarefaction_depth))

# ------------------------------------------------------------
# Repeated pair-specific rarefaction
# ------------------------------------------------------------

results <- vector("list", N_ITER)

for (iter in seq_len(N_ITER)) {

  rare <- counts

  for (pair in unique(meta$Pair_ID)) {

    idx <- which(meta$Pair_ID == pair)

    if (length(idx) != 2) {
      stop(paste("Pair does not contain 2 samples:", pair))
    }

    target <- min(libsize[idx])

    for (i in idx) {

      if (libsize[i] > target) {

        rare[i, ] <- as.numeric(
          rrarefy(
            matrix(
              counts[i, ],
              nrow = 1
            ),
            sample = target
          )
        )

      } else {

        rare[i, ] <- counts[i, ]
      }
    }
  }

  observed <- specnumber(rare)

  shannon <- diversity(
    rare,
    index = "shannon"
  )

  simpson <- diversity(
    rare,
    index = "simpson"
  )

  tmp <- data.frame(
    Pair_ID = meta$Pair_ID,
    STATUS = meta$STATUS,
    Observed = observed,
    Shannon = shannon,
    Simpson = simpson,
    stringsAsFactors = FALSE
  )

  wide <- tmp %>%
    pivot_wider(
      names_from = STATUS,
      values_from = c(
        Observed,
        Shannon,
        Simpson
      )
    )

  get_result <- function(metric) {

    nonpma <- wide[[paste0(metric, "_non-PMA")]]
    pma <- wide[[paste0(metric, "_PMA")]]

    wt <- suppressWarnings(
      wilcox.test(
        pma,
        nonpma,
        paired = TRUE,
        exact = FALSE
      )
    )

    delta <- pma - nonpma

    data.frame(
      Iteration = iter,
      Metric = metric,
      Median_nonPMA = median(nonpma),
      Median_PMA = median(pma),
      Median_difference_PMA_minus_nonPMA =
        median(delta),
      P_value = wt$p.value,
      PMA_higher_pairs = sum(delta > 0),
      PMA_lower_pairs = sum(delta < 0),
      Equal_pairs = sum(delta == 0)
    )
  }

  results[[iter]] <- bind_rows(
    get_result("Observed"),
    get_result("Shannon"),
    get_result("Simpson")
  )
}

results <- bind_rows(results)

write.table(
  results,
  file.path(
    OUTDIR,
    "repeated_pair_rarefaction_results.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------
# Summarize repeated analyses
# ------------------------------------------------------------

summary_results <- results %>%
  group_by(Metric) %>%
  summarise(
    Iterations = n(),

    Median_nonPMA_across_iterations =
      median(Median_nonPMA),

    Median_PMA_across_iterations =
      median(Median_PMA),

    Median_paired_difference =
      median(Median_difference_PMA_minus_nonPMA),

    Median_P_value =
      median(P_value),

    Min_P_value =
      min(P_value),

    Max_P_value =
      max(P_value),

    Fraction_P_lt_0_05 =
      mean(P_value < 0.05),

    Median_PMA_higher_pairs =
      median(PMA_higher_pairs),

    Median_PMA_lower_pairs =
      median(PMA_lower_pairs),

    .groups = "drop"
  )

cat("\n========================================\n")
cat("REPEATED RAREFACTION RESULTS\n")
cat("========================================\n")

print(
  summary_results,
  width = Inf
)

write.table(
  summary_results,
  file.path(
    OUTDIR,
    "repeated_pair_rarefaction_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("\nSaved:\n")
cat(
  file.path(
    OUTDIR,
    "pair_specific_rarefaction_depths.tsv"
  ),
  "\n"
)

cat(
  file.path(
    OUTDIR,
    "repeated_pair_rarefaction_results.tsv"
  ),
  "\n"
)

cat(
  file.path(
    OUTDIR,
    "repeated_pair_rarefaction_summary.tsv"
  ),
  "\n"
)

cat("\nDONE\n")
