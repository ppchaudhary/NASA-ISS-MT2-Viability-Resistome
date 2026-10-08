
# ============================================================
# DeepARG ARG richness depth-sensitivity analysis
# NASA MT2 PMA vs non-PMA
# ============================================================

.libPaths(c(path.expand("~/R/4.5"), .libPaths()))

suppressPackageStartupMessages({
  library(dplyr)
  library(glmmTMB)
})

INFILE <- "analysis_ready/deeparg/figure4/DeepARG_sample_resistome_metrics.tsv"
OUTDIR <- "analysis_ready/deeparg/figure4"

dat <- read.delim(
  INFILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

dat$STATUS <- factor(
  dat$STATUS,
  levels = c("non-PMA", "PMA")
)

dat$Pair_ID <- factor(dat$Pair_ID)

dat$log10_depth <- log10(
  dat$HostRemoved_read_pairs + 1
)

fit_metric <- function(response) {

  f0 <- as.formula(
    paste0(
      response,
      " ~ STATUS + (1|Pair_ID)"
    )
  )

  f1 <- as.formula(
    paste0(
      response,
      " ~ STATUS + log10_depth + (1|Pair_ID)"
    )
  )

  fdepth <- as.formula(
    paste0(
      response,
      " ~ log10_depth + (1|Pair_ID)"
    )
  )

  m0 <- glmmTMB(
    f0,
    family = nbinom2,
    data = dat
  )

  m1 <- glmmTMB(
    f1,
    family = nbinom2,
    data = dat
  )

  mdepth <- glmmTMB(
    fdepth,
    family = nbinom2,
    data = dat
  )

  co <- summary(m1)$coefficients$cond

  status_row <- co["STATUSPMA", ]
  depth_row  <- co["log10_depth", ]

  lrt <- anova(
    mdepth,
    m1,
    test = "LRT"
  )

  result <- data.frame(
    Metric = response,

    PMA_RR =
      exp(status_row["Estimate"]),

    PMA_CI_low =
      exp(
        status_row["Estimate"] -
        1.96 * status_row["Std. Error"]
      ),

    PMA_CI_high =
      exp(
        status_row["Estimate"] +
        1.96 * status_row["Std. Error"]
      ),

    PMA_P =
      status_row["Pr(>|z|)"],

    Depth_RR_per_10x =
      exp(depth_row["Estimate"]),

    Depth_CI_low =
      exp(
        depth_row["Estimate"] -
        1.96 * depth_row["Std. Error"]
      ),

    Depth_CI_high =
      exp(
        depth_row["Estimate"] +
        1.96 * depth_row["Std. Error"]
      ),

    Depth_P =
      depth_row["Pr(>|z|)"],

    LRT_PMA_given_depth =
      lrt$`Pr(>Chisq)`[2],

    AIC_unadjusted =
      AIC(m0),

    AIC_depth_adjusted =
      AIC(m1),

    AIC_depth_only =
      AIC(mdepth),

    pdHess =
      m1$sdr$pdHess
  )

  list(
    result = result,
    m0 = m0,
    m1 = m1,
    mdepth = mdepth
  )
}

subtype <- fit_metric(
  "ARG_subtype_richness"
)

category <- fit_metric(
  "ARG_category_richness"
)

results <- bind_rows(
  subtype$result,
  category$result
)

results$FDR_PMA <- p.adjust(
  results$PMA_P,
  method = "BH"
)

write.table(
  results,
  file.path(
    OUTDIR,
    "DeepARG_richness_depth_adjusted.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

capture.output(
  summary(subtype$m1),
  file = file.path(
    OUTDIR,
    "DeepARG_subtype_richness_depth_model.txt"
  )
)

capture.output(
  summary(category$m1),
  file = file.path(
    OUTDIR,
    "DeepARG_category_richness_depth_model.txt"
  )
)

cat("\n===== DEPTH-ADJUSTED ARG RICHNESS =====\n")
print(results)

cat("\nInterpretation:\n")
cat("RR < 1 for STATUSPMA = lower richness in PMA after including sequencing depth.\n")
cat("Depth RR = multiplicative richness change per 10-fold increase in host-removed reads.\n")
cat("pdHess should be TRUE for a well-behaved fitted model.\n")

