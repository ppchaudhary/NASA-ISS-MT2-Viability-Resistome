
# ============================================================
# DeepARG ARG-group persistence analysis
# NASA MT2 PMA vs non-PMA
#
# Goal:
# Identify specific ARG groups that remain recurrently detected
# after PMA treatment and quantify paired abundance changes.
#
# PMA persistence is detection-based and should NOT be
# interpreted as proof of ARG expression or activity.
# ============================================================

.libPaths(c(path.expand("~/R/4.5"), .libPaths()))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

INFILE <- "analysis_ready/deeparg/deeparg_subtype_long.tsv"
OUTDIR <- "analysis_ready/deeparg/figure4"

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

x <- read.delim(
  INFILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

samples <- x %>%
  distinct(Run, STATUS, SampleName, Pair_ID)

args <- sort(unique(x$`ARG-group`))

cat("\n===== INPUT =====\n")
cat("ARG groups:", length(args), "\n")
cat("Samples:", nrow(samples), "\n")
cat("Pairs:", length(unique(samples$Pair_ID)), "\n")

# ------------------------------------------------------------
# Complete ARG × sample matrix
# Missing ARG/sample combinations become zero.
# ------------------------------------------------------------

complete <- expand_grid(
  Run = samples$Run,
  `ARG-group` = args
) %>%
  left_join(samples, by = "Run") %>%
  left_join(
    x %>%
      select(
        Run,
        `ARG-group`,
        ReadCount,
        `16s-NormalizedReadCount`
      ),
    by = c("Run", "ARG-group")
  ) %>%
  mutate(
    ReadCount = replace_na(ReadCount, 0),
    `16s-NormalizedReadCount` =
      replace_na(`16s-NormalizedReadCount`, 0)
  )

# ------------------------------------------------------------
# Prevalence and abundance
# ------------------------------------------------------------

prev <- complete %>%
  group_by(`ARG-group`, STATUS) %>%
  summarise(
    N_detected = sum(ReadCount > 0),
    Prevalence = mean(ReadCount > 0),
    Median_16S = median(`16s-NormalizedReadCount`),
    Mean_16S = mean(`16s-NormalizedReadCount`),
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
# Paired abundance analysis
# ------------------------------------------------------------

wide <- complete %>%
  select(
    Pair_ID,
    STATUS,
    `ARG-group`,
    `16s-NormalizedReadCount`
  ) %>%
  pivot_wider(
    names_from = STATUS,
    values_from = `16s-NormalizedReadCount`,
    values_fill = 0
  )

tests <- lapply(
  args,
  function(a) {

    z <- wide %>%
      filter(`ARG-group` == a)

    delta <- z$PMA - z$`non-PMA`

    if (all(delta == 0)) {
      V <- NA_real_
      p <- 1
    } else {
      wt <- suppressWarnings(
        wilcox.test(
          z$PMA,
          z$`non-PMA`,
          paired = TRUE,
          exact = FALSE
        )
      )

      V <- unname(wt$statistic)
      p <- wt$p.value
    }

    data.frame(
      ARG_group = a,
      N_pairs = nrow(z),
      Median_nonPMA = median(z$`non-PMA`),
      Median_PMA = median(z$PMA),
      Median_delta = median(delta),

      Median_log2_ratio = median(
        log2(
          (z$PMA + 1e-6) /
          (z$`non-PMA` + 1e-6)
        )
      ),

      PMA_higher = sum(delta > 0),
      PMA_lower = sum(delta < 0),
      Equal = sum(delta == 0),

      Wilcoxon_V = V,
      P_value = p
    )
  }
) %>%
  bind_rows()

tests$FDR_BH <- p.adjust(
  tests$P_value,
  method = "BH"
)

# ------------------------------------------------------------
# Combine
# ------------------------------------------------------------

res <- tests %>%
  left_join(
    prev,
    by = c("ARG_group" = "ARG-group")
  ) %>%
  mutate(
    Prevalence_change =
      Prevalence_PMA - `Prevalence_non-PMA`,

    Persistent_PMA_50 =
      Prevalence_PMA >= 0.50,

    Persistent_PMA_75 =
      Prevalence_PMA >= 0.75,

    Direction = case_when(
      FDR_BH < 0.05 & Median_delta < 0 ~ "Lower in PMA",
      FDR_BH < 0.05 & Median_delta > 0 ~ "Higher in PMA",
      TRUE ~ "Not significant"
    )
  ) %>%
  arrange(
    desc(Prevalence_PMA),
    FDR_BH
  ) %>%
  as_tibble()

write.table(
  res,
  file.path(
    OUTDIR,
    "DeepARG_ARG_group_persistence.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Persistent subset
# ------------------------------------------------------------

persistent <- res %>%
  filter(Persistent_PMA_50) %>%
  arrange(
    desc(Prevalence_PMA),
    desc(Mean_16S_PMA)
  )

write.table(
  persistent,
  file.path(
    OUTDIR,
    "DeepARG_ARG_group_persistent_PMA50.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

cat("\n===== ARG PERSISTENCE SUMMARY =====\n")

cat(
  "ARG groups tested:",
  nrow(res), "\n"
)

cat(
  "Detected in >=50% PMA samples:",
  sum(res$Persistent_PMA_50, na.rm = TRUE),
  "\n"
)

cat(
  "Detected in >=75% PMA samples:",
  sum(res$Persistent_PMA_75, na.rm = TRUE),
  "\n"
)

cat(
  "FDR < 0.05 lower in PMA:",
  sum(res$Direction == "Lower in PMA"),
  "\n"
)

cat(
  "FDR < 0.05 higher in PMA:",
  sum(res$Direction == "Higher in PMA"),
  "\n"
)

cat("\n===== TOP PMA-PERSISTENT ARG GROUPS =====\n")

print(
  persistent %>%
    select(
      ARG_group,
      `Prevalence_non-PMA`,
      Prevalence_PMA,
      Median_nonPMA,
      Median_PMA,
      Median_log2_ratio,
      FDR_BH,
      Direction
    ) %>%
    head(25),
  n = 25
)

cat("\n===== ARG GROUPS >=75% PMA PREVALENCE =====\n")

print(
  res %>%
    filter(Persistent_PMA_75) %>%
    select(
      ARG_group,
      `Prevalence_non-PMA`,
      Prevalence_PMA,
      Median_nonPMA,
      Median_PMA,
      FDR_BH,
      Direction
    ),
  n = Inf
)

