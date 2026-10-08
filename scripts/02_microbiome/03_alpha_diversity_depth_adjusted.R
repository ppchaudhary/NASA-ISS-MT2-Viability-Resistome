#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(glmmTMB)
  library(dplyr)
  library(emmeans)
})

ROOT <- getwd()

ALPHA <- file.path(
  ROOT,
  "analysis_ready/taxonomy/alpha_diversity/alpha_diversity_metrics.tsv"
)

BRACKEN_QC <- file.path(
  ROOT,
  "analysis_ready/taxonomy/bracken_sample_QC.tsv"
)

OUTDIR <- file.path(
  ROOT,
  "analysis_ready/taxonomy/alpha_diversity"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load data
# ============================================================

alpha <- read.delim(
  ALPHA,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

qc <- read.delim(
  BRACKEN_QC,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("========================================\n")
cat("DEPTH-ADJUSTED RICHNESS ANALYSIS\n")
cat("========================================\n")

cat("Alpha rows:", nrow(alpha), "\n")
cat("QC rows:", nrow(qc), "\n")

cat("\nQC columns:\n")
print(colnames(qc))

# ============================================================
# Identify Bracken library-size column
# ============================================================

candidate_cols <- c(
  "Bracken_total_estimated_reads",
  "Bracken_total_reads",
  "Total_Bracken_reads"
)

depth_col <- candidate_cols[
  candidate_cols %in% colnames(qc)
][1]

if (is.na(depth_col)) {
  stop(
    paste(
      "Could not identify Bracken depth column.",
      "Available columns:",
      paste(colnames(qc), collapse = ", ")
    )
  )
}

cat("\nUsing Bracken depth column:", depth_col, "\n")

depth_df <- qc %>%
  select(
    Run,
    Bracken_depth = all_of(depth_col)
  )

dat <- alpha %>%
  left_join(
    depth_df,
    by = "Run"
  )

if (any(is.na(dat$Bracken_depth))) {
  stop("Missing Bracken depth after merge.")
}

if (any(dat$Bracken_depth <= 0)) {
  stop("Bracken depth contains zero/negative values.")
}

dat$STATUS <- factor(
  dat$STATUS,
  levels = c("non-PMA", "PMA")
)

dat$Pair_ID <- factor(dat$Pair_ID)

dat$log10_Bracken_depth <- log10(
  dat$Bracken_depth
)

cat("Samples used:", nrow(dat), "\n")
cat("Pairs:", nlevels(dat$Pair_ID), "\n")

cat("\nStatus counts:\n")
print(table(dat$STATUS))

# ============================================================
# Model 0:
# PMA effect with matched-pair random intercept
# ============================================================

m0 <- glmmTMB(
  Observed ~ STATUS + (1 | Pair_ID),
  family = nbinom2,
  data = dat
)

# ============================================================
# Model 1:
# PMA + sequencing depth + matched-pair random intercept
# ============================================================

m1 <- glmmTMB(
  Observed ~ STATUS + log10_Bracken_depth + (1 | Pair_ID),
  family = nbinom2,
  data = dat
)

# ============================================================
# Model diagnostics
# ============================================================

pearson_overdispersion <- function(model) {

  rdf <- df.residual(model)

  rp <- residuals(
    model,
    type = "pearson"
  )

  Pearson_chisq <- sum(rp^2)

  ratio <- Pearson_chisq / rdf

  p <- pchisq(
    Pearson_chisq,
    df = rdf,
    lower.tail = FALSE
  )

  data.frame(
    Pearson_chisq = Pearson_chisq,
    Residual_df = rdf,
    Dispersion_ratio = ratio,
    Dispersion_P = p
  )
}

diag0 <- pearson_overdispersion(m0)
diag1 <- pearson_overdispersion(m1)

cat("\n========================================\n")
cat("MODEL DIAGNOSTICS\n")
cat("========================================\n")

cat("\nUnadjusted paired NB model:\n")
print(diag0)

cat("\nDepth-adjusted paired NB model:\n")
print(diag1)

# ============================================================
# Coefficient extraction
# ============================================================

coef_table <- function(model, model_name) {

  x <- summary(model)$coefficients$cond

  out <- data.frame(
    Model = model_name,
    Term = rownames(x),
    Estimate_log = x[, "Estimate"],
    SE = x[, "Std. Error"],
    Z = x[, "z value"],
    P_value = x[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )

  out$Rate_ratio <- exp(out$Estimate_log)

  out$CI_low <- exp(
    out$Estimate_log - 1.96 * out$SE
  )

  out$CI_high <- exp(
    out$Estimate_log + 1.96 * out$SE
  )

  out
}

coef_results <- bind_rows(
  coef_table(
    m0,
    "Paired_NB_unadjusted"
  ),
  coef_table(
    m1,
    "Paired_NB_depth_adjusted"
  )
)

cat("\n========================================\n")
cat("MODEL COEFFICIENTS\n")
cat("Rate ratio < 1 for STATUSPMA = lower richness in PMA\n")
cat("========================================\n")

print(
  coef_results,
  row.names = FALSE
)

# ============================================================
# Estimated marginal means from adjusted model
# ============================================================

emm <- emmeans(
  m1,
  ~ STATUS,
  type = "response"
)

emm_table <- as.data.frame(emm)

contrast_result <- contrast(
  emm,
  method = "revpairwise"
)

contrast_table <- as.data.frame(
  summary(
    contrast_result,
    infer = c(TRUE, TRUE),
    type = "response"
  )
)

cat("\n========================================\n")
cat("DEPTH-ADJUSTED ESTIMATED MARGINAL MEANS\n")
cat("========================================\n")

print(emm_table)

cat("\n========================================\n")
cat("DEPTH-ADJUSTED PMA CONTRAST\n")
cat("========================================\n")

print(contrast_table)

# ============================================================
# Likelihood-ratio test:
# does PMA improve model after accounting for depth?
# ============================================================

m_depth_only <- glmmTMB(
  Observed ~ log10_Bracken_depth + (1 | Pair_ID),
  family = nbinom2,
  data = dat
)

lrt <- anova(
  m_depth_only,
  m1
)

cat("\n========================================\n")
cat("LIKELIHOOD-RATIO TEST\n")
cat("Depth-only model vs Depth + PMA model\n")
cat("========================================\n")

print(lrt)

# ============================================================
# Model fit
# ============================================================

fit_table <- data.frame(
  Model = c(
    "Paired_NB_unadjusted",
    "Paired_NB_depth_adjusted",
    "Paired_NB_depth_only"
  ),
  AIC = c(
    AIC(m0),
    AIC(m1),
    AIC(m_depth_only)
  ),
  BIC = c(
    BIC(m0),
    BIC(m1),
    BIC(m_depth_only)
  )
)

cat("\n========================================\n")
cat("MODEL FIT\n")
cat("========================================\n")

print(fit_table)

# ============================================================
# Save results
# ============================================================

write.table(
  dat,
  file.path(
    OUTDIR,
    "richness_depth_adjusted_input.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  coef_results,
  file.path(
    OUTDIR,
    "richness_NB_model_coefficients.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  bind_rows(
    cbind(Model = "Unadjusted", diag0),
    cbind(Model = "Depth_adjusted", diag1)
  ),
  file.path(
    OUTDIR,
    "richness_NB_model_diagnostics.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  emm_table,
  file.path(
    OUTDIR,
    "richness_depth_adjusted_emmeans.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  contrast_table,
  file.path(
    OUTDIR,
    "richness_depth_adjusted_contrast.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

capture.output(
  lrt,
  file = file.path(
    OUTDIR,
    "richness_depth_adjusted_LRT.txt"
  )
)

write.table(
  fit_table,
  file.path(
    OUTDIR,
    "richness_NB_model_fit.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

saveRDS(
  m1,
  file.path(
    OUTDIR,
    "richness_depth_adjusted_model.rds"
  )
)

cat("\n========================================\n")
cat("DONE\n")
cat("========================================\n")

