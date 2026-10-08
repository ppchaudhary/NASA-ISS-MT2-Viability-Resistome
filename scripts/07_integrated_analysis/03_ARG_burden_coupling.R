library(tidyverse)

indir  <- "analysis_ready/integrated_analysis"
outdir <- "analysis_ready/integrated_analysis"
figdir <- "figures/integrated_analysis"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# INPUTS
# ============================================================

dist <- read.delim(
  file.path(
    indir,
    "paired_distances_with_depth.tsv"
  ),
  check.names = FALSE
)

metrics_file <- paste0(
  "analysis_ready/deeparg/figure4/",
  "DeepARG_sample_resistome_metrics.tsv"
)

arg <- read.delim(
  metrics_file,
  check.names = FALSE
)

cat("\n===== DEEPARG METRIC COLUMNS =====\n")
print(names(arg))

# ============================================================
# Identify total 16S-normalized ARG burden column
# ============================================================

candidates <- c(
  "DeepARG_total_16S",
  "Total_16S_normalized",
  "total_16S_normalized",
  "Total16S"
)

burden_col <- candidates[
  candidates %in% names(arg)
][1]

if (is.na(burden_col)) {

  possible <- grep(
    "total.*16|16.*total",
    names(arg),
    ignore.case = TRUE,
    value = TRUE
  )

  if (length(possible) != 1) {
    stop(
      paste(
        "Could not uniquely identify total",
        "DeepARG 16S-normalized burden column.",
        "Inspect columns printed above."
      )
    )
  }

  burden_col <- possible[1]
}

cat(
  "\nUsing ARG burden column:",
  burden_col,
  "\n"
)

# ============================================================
# Build paired ARG-burden table
# ============================================================

burden <- arg %>%
  select(
    Pair_ID,
    STATUS,
    all_of(burden_col)
  ) %>%
  rename(
    ARG_burden = all_of(burden_col)
  ) %>%
  mutate(
    ARG_burden = as.numeric(ARG_burden)
  ) %>%
  pivot_wider(
    names_from = STATUS,
    values_from = ARG_burden
  )

if (
  !all(
    c("PMA", "non-PMA") %in% names(burden)
  )
) {
  stop("Could not construct paired PMA/non-PMA ARG burden.")
}

burden <- burden %>%
  mutate(
    # Signed change
    Delta_ARG_burden =
      PMA - `non-PMA`,

    # Magnitude of change
    Abs_Delta_ARG_burden =
      abs(Delta_ARG_burden),

    # Ratio; useful descriptive measure
    ARG_burden_ratio =
      if_else(
        `non-PMA` > 0,
        PMA / `non-PMA`,
        NA_real_
      ),

    Log2_ARG_burden_ratio =
      if_else(
        PMA > 0 & `non-PMA` > 0,
        log2(PMA / `non-PMA`),
        NA_real_
      )
  )

dat <- dist %>%
  left_join(
    burden,
    by = "Pair_ID"
  )

if (nrow(dat) != 32) {
  stop(
    paste(
      "Expected 32 pairs; found",
      nrow(dat)
    )
  )
}

# ============================================================
# Correlation helper
# ============================================================

run_spearman <- function(x, y, label) {

  ok <- complete.cases(x, y)

  test <- cor.test(
    x[ok],
    y[ok],
    method = "spearman",
    exact = FALSE
  )

  tibble(
    Analysis = label,
    N = sum(ok),
    Spearman_rho =
      unname(test$estimate),
    P_value =
      test$p.value
  )
}

# ============================================================
# PRIMARY / SECONDARY ASSOCIATIONS
# ============================================================

results <- bind_rows(

  # ----------------------------------------------------------
  # Signed burden change
  # Does the direction/magnitude of PMA burden change track
  # community/resistome compositional displacement?
  # ----------------------------------------------------------

  run_spearman(
    dat$Microbiome_BrayCurtis,
    dat$Delta_ARG_burden,
    "Microbiome distance vs signed ARG burden change"
  ),

  run_spearman(
    dat$Resistome_BrayCurtis,
    dat$Delta_ARG_burden,
    "Resistome distance vs signed ARG burden change"
  ),

  # ----------------------------------------------------------
  # Absolute burden change
  # Better matched conceptually to an unsigned distance metric.
  # ----------------------------------------------------------

  run_spearman(
    dat$Microbiome_BrayCurtis,
    dat$Abs_Delta_ARG_burden,
    "Microbiome distance vs absolute ARG burden change"
  ),

  run_spearman(
    dat$Resistome_BrayCurtis,
    dat$Abs_Delta_ARG_burden,
    "Resistome distance vs absolute ARG burden change"
  ),

  # ----------------------------------------------------------
  # Fold-change sensitivity
  # ----------------------------------------------------------

  run_spearman(
    dat$Microbiome_BrayCurtis,
    dat$Log2_ARG_burden_ratio,
    "Microbiome distance vs log2 ARG burden ratio"
  ),

  run_spearman(
    dat$Resistome_BrayCurtis,
    dat$Log2_ARG_burden_ratio,
    "Resistome distance vs log2 ARG burden ratio"
  )
)

# BH adjustment across these 6 exploratory tests
results <- results %>%
  mutate(
    FDR = p.adjust(
      P_value,
      method = "BH"
    )
  )

# ============================================================
# Depth sensitivity for burden-change relationships
#
# Rank regression:
# burden-change rank ~ biological-distance rank + depth rank
# ============================================================

rank_sensitivity <- function(
    response,
    biological,
    depth,
    label
) {

  dd <- tibble(
    response = response,
    biological = biological,
    depth = depth
  ) %>%
    drop_na() %>%
    mutate(
      rank_response =
        rank(response),
      rank_biological =
        rank(biological),
      rank_depth =
        rank(depth)
    )

  fit <- lm(
    rank_response ~
      rank_biological +
      rank_depth,
    data = dd
  )

  sm <- summary(fit)$coefficients

  tibble(
    Analysis = label,
    N = nrow(dd),
    Biological_beta =
      sm[
        "rank_biological",
        "Estimate"
      ],
    Biological_P =
      sm[
        "rank_biological",
        "Pr(>|t|)"
      ],
    Depth_beta =
      sm[
        "rank_depth",
        "Estimate"
      ],
    Depth_P =
      sm[
        "rank_depth",
        "Pr(>|t|)"
      ],
    R_squared =
      summary(fit)$r.squared
  )
}

depth_results <- bind_rows(

  rank_sensitivity(
    dat$Delta_ARG_burden,
    dat$Microbiome_BrayCurtis,
    dat$Abs_log10_depth_imbalance,
    "Signed burden change ~ microbiome distance + depth"
  ),

  rank_sensitivity(
    dat$Abs_Delta_ARG_burden,
    dat$Microbiome_BrayCurtis,
    dat$Abs_log10_depth_imbalance,
    "Absolute burden change ~ microbiome distance + depth"
  ),

  rank_sensitivity(
    dat$Delta_ARG_burden,
    dat$Resistome_BrayCurtis,
    dat$Abs_log10_depth_imbalance,
    "Signed burden change ~ resistome distance + depth"
  ),

  rank_sensitivity(
    dat$Abs_Delta_ARG_burden,
    dat$Resistome_BrayCurtis,
    dat$Abs_log10_depth_imbalance,
    "Absolute burden change ~ resistome distance + depth"
  )
)

# ============================================================
# SAVE
# ============================================================

write.table(
  dat,
  file.path(
    outdir,
    "paired_ARG_burden_changes.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  results,
  file.path(
    outdir,
    "ARG_burden_coupling_results.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  depth_results,
  file.path(
    outdir,
    "ARG_burden_coupling_depth_sensitivity.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# EXPLORATORY PLOT
# Microbiome distance vs signed burden change
# ============================================================

r1 <- results %>%
  filter(
    Analysis ==
      "Microbiome distance vs signed ARG burden change"
  )

lab1 <- paste0(
  "Spearman \u03c1 = ",
  sprintf(
    "%.2f",
    r1$Spearman_rho
  ),
  "\nP = ",
  format.pval(
    r1$P_value,
    digits = 2,
    eps = 0.001
  )
)

p1 <- ggplot(
  dat,
  aes(
    x = Microbiome_BrayCurtis,
    y = Delta_ARG_burden
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = 2,
    color = "grey60"
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
    label = lab1,
    hjust = 1.1,
    vjust = 1.2,
    size = 4
  ) +

  labs(
    x = "Within-pair microbiome Bray\u2013Curtis distance",
    y = "\u0394 DeepARG 16S-normalized burden\n(PMA \u2212 non-PMA)"
  ) +

  theme_classic(base_size = 12)

ggsave(
  file.path(
    figdir,
    "ARG_burden_vs_microbiome_shift_exploratory.pdf"
  ),
  p1,
  width = 6.2,
  height = 5.2,
  device = cairo_pdf
)

ggsave(
  file.path(
    figdir,
    "ARG_burden_vs_microbiome_shift_exploratory.png"
  ),
  p1,
  width = 6.2,
  height = 5.2,
  dpi = 400
)

# ============================================================
# REPORT
# ============================================================

cat("\n============================================\n")
cat("ARG BURDEN COUPLING ANALYSIS\n")
cat("============================================\n\n")

cat("Paired burden-change summary:\n\n")

print(
  summary(
    dat$Delta_ARG_burden
  )
)

cat("\nPMA/non-PMA burden ratio:\n\n")

print(
  summary(
    dat$ARG_burden_ratio
  )
)

cat("\nAssociation tests:\n\n")
print(
  results,
  n = Inf
)

cat("\nDepth-sensitivity rank regressions:\n\n")
print(
  depth_results,
  n = Inf
)

cat("\n============================================\n")
