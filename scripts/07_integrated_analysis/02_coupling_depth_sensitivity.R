library(tidyverse)

indir  <- "analysis_ready/integrated_analysis"
outdir <- "analysis_ready/integrated_analysis"
figdir <- "figures/integrated_analysis"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# INPUT
# ============================================================

dist <- read.delim(
  file.path(
    indir,
    "paired_microbiome_resistome_distances.tsv"
  ),
  check.names = FALSE
)

master <- read.delim(
  "analysis_ready/metadata/master_sample_table.tsv",
  check.names = FALSE
)

cat("\n===== DEPTH COLUMNS AVAILABLE =====\n")
print(
  grep(
    "Host|host|read|Read|depth|Depth",
    names(master),
    value = TRUE
  )
)

# ============================================================
# Identify host-removed sequencing-depth column
# ============================================================

candidate_names <- c(
  "HostRemoved_read_pairs",
  "HostRemoved_pairs",
  "Host_removed_pairs",
  "HostRemovedPairs",
  "host_removed_pairs"
)

depth_col <- candidate_names[
  candidate_names %in% names(master)
][1]

if (is.na(depth_col)) {

  possible <- grep(
    "host.*remov.*pair|hostremoved",
    names(master),
    ignore.case = TRUE,
    value = TRUE
  )

  if (length(possible) == 0) {
    stop(
      paste(
        "Could not identify host-removed depth column.",
        "Inspect the printed column names above."
      )
    )
  }

  depth_col <- possible[1]
}

cat(
  "\nUsing depth column:",
  depth_col,
  "\n"
)

# ============================================================
# Construct pair-level depth metrics
# ============================================================

depth <- master %>%
  select(
    Run,
    Pair_ID,
    STATUS,
    all_of(depth_col)
  ) %>%
  rename(
    HostRemoved_pairs = all_of(depth_col)
  ) %>%
  mutate(
    HostRemoved_pairs = as.numeric(
      HostRemoved_pairs
    )
  )

depth_pair <- depth %>%
  select(
    Pair_ID,
    STATUS,
    HostRemoved_pairs
  ) %>%
  pivot_wider(
    names_from = STATUS,
    values_from = HostRemoved_pairs
  )

if (
  !all(
    c("PMA", "non-PMA") %in%
      names(depth_pair)
  )
) {
  stop("Could not construct PMA/non-PMA depth pairs.")
}

depth_pair <- depth_pair %>%
  mutate(
    Depth_ratio_PMA_nonPMA =
      PMA / `non-PMA`,

    Log10_depth_ratio =
      log10(PMA / `non-PMA`),

    # Primary sensitivity covariate:
    # magnitude of within-pair depth inequality
    Abs_log10_depth_imbalance =
      abs(Log10_depth_ratio)
  )

dat <- dist %>%
  left_join(
    depth_pair,
    by = "Pair_ID"
  )

if (anyNA(dat$Abs_log10_depth_imbalance)) {
  stop("Missing pair-level depth information.")
}

# ============================================================
# Descriptive depth imbalance
# ============================================================

cat("\n===== WITHIN-PAIR DEPTH RATIO =====\n")

print(
  summary(
    dat$Depth_ratio_PMA_nonPMA
  )
)

cat("\nAbsolute log10 depth imbalance:\n")

print(
  summary(
    dat$Abs_log10_depth_imbalance
  )
)

# ============================================================
# 1. Does depth imbalance predict microbiome distance?
# ============================================================

micro_depth <- cor.test(
  dat$Microbiome_BrayCurtis,
  dat$Abs_log10_depth_imbalance,
  method = "spearman",
  exact = FALSE
)

# ============================================================
# 2. Does depth imbalance predict resistome distance?
# ============================================================

res_depth <- cor.test(
  dat$Resistome_BrayCurtis,
  dat$Abs_log10_depth_imbalance,
  method = "spearman",
  exact = FALSE
)

# ============================================================
# 3. Original microbiome-resistome association
# ============================================================

original <- cor.test(
  dat$Microbiome_BrayCurtis,
  dat$Resistome_BrayCurtis,
  method = "spearman",
  exact = FALSE
)

# ============================================================
# 4. Partial Spearman correlation controlling for depth
#
# Rank-transform all variables, regress each biological
# distance on ranked depth imbalance, then correlate residuals.
# ============================================================

dat <- dat %>%
  mutate(
    rank_micro =
      rank(
        Microbiome_BrayCurtis,
        ties.method = "average"
      ),

    rank_res =
      rank(
        Resistome_BrayCurtis,
        ties.method = "average"
      ),

    rank_depth =
      rank(
        Abs_log10_depth_imbalance,
        ties.method = "average"
      )
  )

fit_micro_depth <- lm(
  rank_micro ~ rank_depth,
  data = dat
)

fit_res_depth <- lm(
  rank_res ~ rank_depth,
  data = dat
)

dat <- dat %>%
  mutate(
    micro_depth_residual =
      residuals(fit_micro_depth),

    res_depth_residual =
      residuals(fit_res_depth)
  )

partial_test <- cor.test(
  dat$micro_depth_residual,
  dat$res_depth_residual,
  method = "pearson"
)

# ============================================================
# 5. Rank regression sensitivity model
#
# Resistome rank ~ microbiome rank + depth-imbalance rank
# ============================================================

rank_model <- lm(
  rank_res ~ rank_micro + rank_depth,
  data = dat
)

coef_table <- summary(rank_model)$coefficients

# standardized rank variables for interpretable beta
dat <- dat %>%
  mutate(
    z_rank_micro = as.numeric(
      scale(rank_micro)
    ),
    z_rank_res = as.numeric(
      scale(rank_res)
    ),
    z_rank_depth = as.numeric(
      scale(rank_depth)
    )
  )

std_model <- lm(
  z_rank_res ~
    z_rank_micro +
    z_rank_depth,
  data = dat
)

std_coef <- summary(std_model)$coefficients

# ============================================================
# Save summary
# ============================================================

results <- tibble(
  Analysis = c(
    "Original microbiome-resistome Spearman",
    "Microbiome distance vs depth imbalance",
    "Resistome distance vs depth imbalance",
    "Partial rank correlation controlling depth"
  ),

  Estimate = c(
    unname(original$estimate),
    unname(micro_depth$estimate),
    unname(res_depth$estimate),
    unname(partial_test$estimate)
  ),

  P_value = c(
    original$p.value,
    micro_depth$p.value,
    res_depth$p.value,
    partial_test$p.value
  )
)

write.table(
  results,
  file.path(
    outdir,
    "microbiome_resistome_depth_sensitivity.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  dat,
  file.path(
    outdir,
    "paired_distances_with_depth.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

capture.output(
  summary(rank_model),
  file = file.path(
    outdir,
    "coupling_rank_regression.txt"
  )
)

capture.output(
  summary(std_model),
  file = file.path(
    outdir,
    "coupling_standardized_rank_regression.txt"
  )
)

# ============================================================
# Diagnostic figure
# ============================================================

label_txt <- paste0(
  "Partial rank r = ",
  sprintf(
    "%.2f",
    unname(partial_test$estimate)
  ),
  "\nP = ",
  format.pval(
    partial_test$p.value,
    digits = 2,
    eps = 0.001
  )
)

p <- ggplot(
  dat,
  aes(
    x = micro_depth_residual,
    y = res_depth_residual
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
    size = 4.1
  ) +

  labs(
    x = "Microbiome distance\n(depth-adjusted rank residual)",
    y = "Resistome distance\n(depth-adjusted rank residual)"
  ) +

  theme_classic(base_size = 12) +

  theme(
    axis.text = element_text(
      color = "black"
    ),
    axis.title = element_text(
      color = "black"
    )
  )

ggsave(
  file.path(
    figdir,
    "microbiome_resistome_coupling_depth_adjusted.pdf"
  ),
  p,
  width = 6.2,
  height = 5.2,
  device = cairo_pdf
)

ggsave(
  file.path(
    figdir,
    "microbiome_resistome_coupling_depth_adjusted.png"
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
cat("COUPLING DEPTH-SENSITIVITY ANALYSIS\n")
cat("============================================\n\n")

print(results)

cat("\nStandardized rank regression:\n\n")

print(
  coef(
    summary(std_model)
  )
)

cat(
  "\nModel R-squared =",
  summary(std_model)$r.squared,
  "\n"
)

cat(
  "Adjusted R-squared =",
  summary(std_model)$adj.r.squared,
  "\n"
)

cat("\nInterpretation check:\n")

cat(
  "Microbiome standardized beta =",
  std_coef[
    "z_rank_micro",
    "Estimate"
  ],
  "; P =",
  std_coef[
    "z_rank_micro",
    "Pr(>|t|)"
  ],
  "\n"
)

cat(
  "Depth standardized beta =",
  std_coef[
    "z_rank_depth",
    "Estimate"
  ],
  "; P =",
  std_coef[
    "z_rank_depth",
    "Pr(>|t|)"
  ],
  "\n"
)

cat("\n============================================\n")
