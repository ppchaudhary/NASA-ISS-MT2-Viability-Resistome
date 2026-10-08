
# ============================================================
# DeepARG resistance-category paired analysis
# NASA MT2 PMA vs non-PMA
#
# Quantitative metric:
#   16S-normalized DeepARG abundance
#
# Statistical design:
#   matched PMA/non-PMA samples
#   paired Wilcoxon tests
#   BH multiple-testing correction
# ============================================================

.libPaths(c(path.expand("~/R/4.5"), .libPaths()))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

INFILE <- "analysis_ready/deeparg/deeparg_category_long.tsv"
OUTDIR <- "analysis_ready/deeparg/figure4"

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

x <- read.delim(
  INFILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Complete sample × category table
# ------------------------------------------------------------

samples <- x %>%
  distinct(
    Run,
    STATUS,
    SampleName,
    Pair_ID
  )

categories <- sort(
  unique(x$`ARG-category`)
)

complete <- expand_grid(
  Run = samples$Run,
  `ARG-category` = categories
) %>%
  left_join(
    samples,
    by = "Run"
  ) %>%
  left_join(
    x %>%
      select(
        Run,
        `ARG-category`,
        ReadCount,
        `16s-NormalizedReadCount`
      ),
    by = c(
      "Run",
      "ARG-category"
    )
  ) %>%
  mutate(
    ReadCount =
      replace_na(ReadCount, 0),

    `16s-NormalizedReadCount` =
      replace_na(
        `16s-NormalizedReadCount`,
        0
      )
  )

# ------------------------------------------------------------
# Prevalence
# ------------------------------------------------------------

prevalence <- complete %>%
  group_by(
    `ARG-category`,
    STATUS
  ) %>%
  summarise(
    N_detected =
      sum(ReadCount > 0),

    Prevalence =
      mean(ReadCount > 0),

    Median_16S =
      median(
        `16s-NormalizedReadCount`
      ),

    Mean_16S =
      mean(
        `16s-NormalizedReadCount`
      ),

    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = STATUS,
    values_from = c(
      N_detected,
      Prevalence,
      Median_16S,
      Mean_16S
    )
  )

# ------------------------------------------------------------
# Paired category tests
# ------------------------------------------------------------

paired <- complete %>%
  select(
    Pair_ID,
    STATUS,
    `ARG-category`,
    `16s-NormalizedReadCount`
  ) %>%
  pivot_wider(
    names_from = STATUS,
    values_from =
      `16s-NormalizedReadCount`,
    values_fill = 0
  )

results <- lapply(
  categories,
  function(cat) {

    z <- paired %>%
      filter(
        `ARG-category` == cat
      )

    delta <-
      z$PMA - z$`non-PMA`

    # Paired Wilcoxon cannot provide meaningful
    # inference if every paired value is identical.
    if (
      length(unique(delta)) <= 1 &&
      unique(delta)[1] == 0
    ) {

      p <- 1
      V <- NA_real_

    } else {

      tst <- suppressWarnings(
        wilcox.test(
          z$PMA,
          z$`non-PMA`,
          paired = TRUE,
          exact = FALSE
        )
      )

      p <- tst$p.value
      V <- unname(tst$statistic)
    }

    data.frame(
      ARG_category = cat,

      N_pairs =
        nrow(z),

      Median_nonPMA =
        median(z$`non-PMA`),

      Median_PMA =
        median(z$PMA),

      Median_delta =
        median(delta),

      Median_log2_ratio =
        median(
          log2(
            (z$PMA + 1e-6) /
            (z$`non-PMA` + 1e-6)
          )
        ),

      PMA_higher =
        sum(delta > 0),

      PMA_lower =
        sum(delta < 0),

      Equal =
        sum(delta == 0),

      Wilcoxon_V = V,

      P_value = p
    )
  }
) %>%
  bind_rows()

results$FDR_BH <- p.adjust(
  results$P_value,
  method = "BH"
)

results <- results %>%
  left_join(
    prevalence,
    by = c(
      "ARG_category" =
        "ARG-category"
    )
  ) %>%
  mutate(

    Prevalence_change =
      Prevalence_PMA -
      `Prevalence_non-PMA`,

    Direction =
      case_when(
        FDR_BH < 0.05 &
          Median_delta > 0
          ~ "Higher in PMA",

        FDR_BH < 0.05 &
          Median_delta < 0
          ~ "Lower in PMA",

        TRUE ~ "Not significant"
      )
  ) %>%
  arrange(
    FDR_BH,
    P_value
  )

write.table(
  results,
  file.path(
    OUTDIR,
    "DeepARG_category_paired_results.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  complete,
  file.path(
    OUTDIR,
    "DeepARG_category_complete_long.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Console summary
# ------------------------------------------------------------

results <- tibble::as_tibble(results)

cat("\n===== CATEGORY ANALYSIS =====\n")
cat(
  "Categories tested:",
  nrow(results),
  "\n"
)

cat(
  "FDR < 0.05:",
  sum(
    results$FDR_BH < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Significantly lower in PMA:",
  sum(
    results$Direction ==
      "Lower in PMA",
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Significantly higher in PMA:",
  sum(
    results$Direction ==
      "Higher in PMA",
    na.rm = TRUE
  ),
  "\n"
)

cat("\n===== TOP RESULTS =====\n")

print(
  results %>%
    select(
      ARG_category,
      Median_nonPMA,
      Median_PMA,
      Median_log2_ratio,
      PMA_higher,
      PMA_lower,
      FDR_BH,
      `Prevalence_non-PMA`,
      Prevalence_PMA,
      Direction
    ) %>%
    head(25),
  n = 25
)

cat("\n===== CORE PMA CATEGORIES >=75% PREVALENCE =====\n")

print(
  results %>%
    filter(
      Prevalence_PMA >= 0.75
    ) %>%
    arrange(
      desc(Prevalence_PMA),
      FDR_BH
    ) %>%
    select(
      ARG_category,
      `Prevalence_non-PMA`,
      Prevalence_PMA,
      Median_nonPMA,
      Median_PMA,
      FDR_BH,
      Direction
    ),
  n = Inf
)

