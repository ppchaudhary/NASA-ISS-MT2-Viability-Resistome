suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

outdir <- "analysis_ready/integrated_analysis"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

micro <- read.delim(
  "analysis_ready/taxonomy/alpha_diversity/alpha_diversity_metrics.tsv",
  check.names = FALSE
)

arg <- read.delim(
  "analysis_ready/deeparg/figure4/DeepARG_sample_resistome_metrics.tsv",
  check.names = FALSE
)

metadata <- micro %>%
  distinct(Pair_ID, Flight, LocationCode)

stopifnot(
  nrow(micro) == 64,
  nrow(arg) == 64,
  nrow(metadata) == 32,
  all(table(metadata$LocationCode) == 4),
  all(table(metadata$Flight) == 8)
)

prepare_pairs <- function(df, metric, label) {

  stopifnot(
    all(c("Pair_ID", "STATUS", metric) %in% names(df))
  )

  paired <- df %>%
    select(Pair_ID, STATUS, all_of(metric)) %>%
    pivot_wider(
      names_from = STATUS,
      values_from = all_of(metric)
    )

  stopifnot(
    nrow(paired) == 32,
    all(c("PMA", "non-PMA") %in% names(paired)),
    !anyNA(paired$PMA),
    !anyNA(paired$`non-PMA`)
  )

  paired %>%
    left_join(metadata, by = "Pair_ID") %>%
    mutate(
      Metric = label,
      PMA_minus_nonPMA = PMA - `non-PMA`
    )
}

micro_pairs <- prepare_pairs(
  micro, "Observed", "Microbial_species_richness"
)

arg_pairs <- prepare_pairs(
  arg, "ARG_subtype_richness", "ARG_group_richness"
)

all_pairs <- bind_rows(micro_pairs, arg_pairs)

location_results <- all_pairs %>%
  group_by(Metric, LocationCode) %>%
  summarise(
    N_flights = n(),
    Median_PMA_minus_nonPMA = median(PMA_minus_nonPMA),
    Mean_PMA_minus_nonPMA = mean(PMA_minus_nonPMA),
    Pairs_lower_after_PMA = sum(PMA_minus_nonPMA < 0),
    Pairs_higher_after_PMA = sum(PMA_minus_nonPMA > 0),
    .groups = "drop"
  )

stopifnot(
  nrow(location_results) == 16,
  all(location_results$N_flights == 4)
)

summary_results <- location_results %>%
  group_by(Metric) %>%
  summarise(
    N_locations = n(),
    Locations_negative = sum(Median_PMA_minus_nonPMA < 0),
    Locations_positive = sum(Median_PMA_minus_nonPMA > 0),
    Locations_zero = sum(Median_PMA_minus_nonPMA == 0),
    Median_location_effect = median(Median_PMA_minus_nonPMA),
    .groups = "drop"
  )

summary_results$Exact_sign_test_P <- vapply(
  seq_len(nrow(summary_results)),
  function(i) {
    row <- summary_results[i, ]
    n <- row$Locations_negative + row$Locations_positive
    if (n == 0) return(NA_real_)
    binom.test(
      row$Locations_negative,
      n,
      p = 0.5,
      alternative = "two.sided"
    )$p.value
  },
  numeric(1)
)

write.table(
  all_pairs,
  file.path(outdir, "richness_location_paired_differences.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  location_results,
  file.path(outdir, "richness_location_sensitivity.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  summary_results,
  file.path(outdir, "richness_location_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\n===== LOCATION-LEVEL RICHNESS RESULTS =====\n")
print(location_results, n = Inf)

cat("\n===== SUMMARY =====\n")
print(summary_results, width = Inf)

cat("\n===== STEP 59 COMPLETE =====\n")
