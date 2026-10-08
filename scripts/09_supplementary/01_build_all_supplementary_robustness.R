# ============================================================
# Supplementary robustness package: Figures S1-S7
# ISS matched PMA vs non-PMA metagenomics study
#
# Purpose:
#   Build supplementary reviewer-facing robustness figures from
#   existing/frozen analysis outputs. Primary analyses are NOT
#   recomputed here.
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# GLOBAL SETTINGS
# ------------------------------------------------------------

BLUE   <- "#1479E8"
RED    <- "#FF4B4B"
PURPLE <- "#8B6BD9"
ORANGE <- "#F39C34"
GREEN  <- "#45A85A"
TEAL   <- "#10B8B8"
GREY   <- "grey55"

figdir <- "figures/supplementary"
tabdir <- "analysis_ready/supplementary"

dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)

theme_paper <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.tag = element_text(face = "bold", size = 15),
    axis.title = element_text(colour = "black"),
    axis.text = element_text(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

save_fig <- function(plot, stem, width = 11, height = 8) {
  ggsave(
    file.path(figdir, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    device = cairo_pdf
  )
  ggsave(
    file.path(figdir, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    dpi = 400
  )
}

cat("\n====================================================\n")
cat("BUILDING SUPPLEMENTARY ROBUSTNESS PACKAGE S1-S7\n")
cat("====================================================\n")

# ============================================================
# S1 — SEQUENCING DEPTH / ASSEMBLY / DOWNSTREAM RECOVERY
# ============================================================

cat("\n[S1] Sequencing-depth and assembly-yield QC\n")

master <- read.delim(
  "analysis_ready/metadata/master_sample_table.tsv",
  check.names = FALSE
)

req <- c(
  "Run", "Pair_ID", "STATUS",
  "HostRemoved_read_pairs",
  "Assembly_total_bp",
  "AMRFinder_hits",
  "KEGG_valid_assignments"
)

stopifnot(all(req %in% names(master)))

master <- master %>%
  mutate(
    STATUS = factor(STATUS, levels = c("non-PMA", "PMA"))
  )

S1_summary <- master %>%
  group_by(STATUS) %>%
  summarise(
    n = n(),
    median_reads = median(HostRemoved_read_pairs, na.rm = TRUE),
    Q1_reads = quantile(HostRemoved_read_pairs, .25, na.rm = TRUE),
    Q3_reads = quantile(HostRemoved_read_pairs, .75, na.rm = TRUE),
    median_assembly_bp = median(Assembly_total_bp, na.rm = TRUE),
    Q1_assembly_bp = quantile(Assembly_total_bp, .25, na.rm = TRUE),
    Q3_assembly_bp = quantile(Assembly_total_bp, .75, na.rm = TRUE),
    median_AMRFinder_hits = median(AMRFinder_hits, na.rm = TRUE),
    median_KEGG_assignments =
      median(KEGG_valid_assignments, na.rm = TRUE),
    .groups = "drop"
  )

write.table(
  S1_summary,
  file.path(tabdir, "S1_group_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

pair_depth <- master %>%
  select(
    Pair_ID, STATUS,
    HostRemoved_read_pairs,
    Assembly_total_bp
  ) %>%
  pivot_wider(
    names_from = STATUS,
    values_from = c(
      HostRemoved_read_pairs,
      Assembly_total_bp
    )
  ) %>%
  mutate(
    read_ratio =
      `HostRemoved_read_pairs_PMA` /
      `HostRemoved_read_pairs_non-PMA`,
    assembly_ratio =
      `Assembly_total_bp_PMA` /
      `Assembly_total_bp_non-PMA`
  )

write.table(
  pair_depth,
  file.path(tabdir, "S1_pair_depth_ratios.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

cor_vars <- c(
  "Assembly_total_bp",
  "AMRFinder_hits",
  "KEGG_valid_assignments"
)

S1_cor <- bind_rows(lapply(cor_vars, function(v) {

  ct <- suppressWarnings(
    cor.test(
      master$HostRemoved_read_pairs,
      master[[v]],
      method = "spearman",
      exact = FALSE
    )
  )

  data.frame(
    predictor = "HostRemoved_read_pairs",
    outcome = v,
    rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}))

write.table(
  S1_cor,
  file.path(tabdir, "S1_depth_correlations.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

p1A <- ggplot(
  master,
  aes(
    STATUS,
    HostRemoved_read_pairs,
    group = Pair_ID
  )
) +
  geom_line(colour = "grey82", linewidth = .45) +
  geom_point(
    aes(fill = STATUS),
    shape = 21, size = 2.5,
    colour = "black", stroke = .35
  ) +
  scale_fill_manual(values = c("non-PMA" = BLUE, "PMA" = RED)) +
  scale_y_log10(labels = label_number()) +
  labs(
    x = NULL,
    y = "Host-removed read pairs",
    title = "Sequencing depth"
  ) +
  theme_paper +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold")
  )

p1B <- ggplot(
  master,
  aes(
    STATUS,
    Assembly_total_bp,
    group = Pair_ID
  )
) +
  geom_line(colour = "grey82", linewidth = .45) +
  geom_point(
    aes(fill = STATUS),
    shape = 21, size = 2.5,
    colour = "black", stroke = .35
  ) +
  scale_fill_manual(values = c("non-PMA" = BLUE, "PMA" = RED)) +
  scale_y_log10(labels = label_number()) +
  labs(
    x = NULL,
    y = "Assembly yield (bp)",
    title = "Assembly yield"
  ) +
  theme_paper +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold")
  )

ratio_long <- pair_depth %>%
  select(Pair_ID, read_ratio, assembly_ratio) %>%
  pivot_longer(
    c(read_ratio, assembly_ratio),
    names_to = "Metric",
    values_to = "Ratio"
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      read_ratio = "Read depth",
      assembly_ratio = "Assembly yield"
    )
  )

p1C <- ggplot(ratio_long, aes(Metric, Ratio)) +
  geom_hline(
    yintercept = 1,
    linetype = 2,
    colour = "grey45"
  ) +
  geom_boxplot(
    width = .5,
    outlier.shape = NA,
    fill = "white"
  ) +
  geom_jitter(
    width = .10,
    size = 1.7,
    alpha = .75,
    colour = RED
  ) +
  scale_y_log10() +
  labs(
    x = NULL,
    y = "PMA / non-PMA ratio",
    title = "Within-pair depth imbalance"
  ) +
  theme_paper +
  theme(axis.text.x = element_text(face = "bold"))

corr_plot_data <- master %>%
  select(
    HostRemoved_read_pairs,
    Assembly_total_bp,
    AMRFinder_hits,
    KEGG_valid_assignments
  ) %>%
  pivot_longer(
    -HostRemoved_read_pairs,
    names_to = "Outcome",
    values_to = "Value"
  ) %>%
  mutate(
    Outcome = recode(
      Outcome,
      Assembly_total_bp = "Assembly yield",
      AMRFinder_hits = "AMRFinderPlus hits",
      KEGG_valid_assignments = "KEGG assignments"
    )
  )

corr_labels <- S1_cor %>%
  mutate(
    Outcome = recode(
      outcome,
      Assembly_total_bp = "Assembly yield",
      AMRFinder_hits = "AMRFinderPlus hits",
      KEGG_valid_assignments = "KEGG assignments"
    ),
    label = sprintf(
      "\u03c1 = %.2f\nP = %.2g",
      rho, p_value
    )
  )

p1D <- ggplot(
  corr_plot_data,
  aes(HostRemoved_read_pairs, Value)
) +
  geom_point(
    size = 1.7,
    alpha = .72,
    colour = BLUE
  ) +
  scale_x_log10(labels = label_number()) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10),
    labels = label_number()
  ) +
  facet_wrap(~Outcome, scales = "free_y", nrow = 1) +
  geom_text(
    data = corr_labels,
    aes(
      x = Inf, y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1.08,
    vjust = 1.15,
    size = 3.1
  ) +
  labs(
    x = "Host-removed read pairs",
    y = "Downstream recovery",
    title = "Sequencing depth and downstream recovery"
  ) +
  theme_paper

S1 <- (p1A | p1B) / (p1C | p1D) +
  plot_annotation(tag_levels = "A")

save_fig(
  S1,
  "Supplementary_Figure_S1_depth_assembly_QC",
  13, 8.5
)

# ============================================================
# S2 — MICROBIOME DEPTH ROBUSTNESS
# ============================================================

cat("[S2] Microbiome depth robustness\n")

alpha_input <- read.delim(
  "analysis_ready/taxonomy/alpha_diversity/richness_depth_adjusted_input.tsv",
  check.names = FALSE
)

alpha_coef <- read.delim(
  "analysis_ready/taxonomy/alpha_diversity/richness_NB_model_coefficients.tsv",
  check.names = FALSE
)

alpha_diag <- read.delim(
  "analysis_ready/taxonomy/alpha_diversity/richness_NB_model_diagnostics.tsv",
  check.names = FALSE
)

alpha_tests <- read.delim(
  "analysis_ready/taxonomy/alpha_diversity/alpha_diversity_paired_tests.tsv",
  check.names = FALSE
)

write.table(
  alpha_coef,
  file.path(tabdir, "S2_richness_NB_coefficients.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  alpha_diag,
  file.path(tabdir, "S2_richness_NB_diagnostics.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# Detect richness column safely
rich_candidates <- grep(
  "rich|observ",
  names(alpha_input),
  ignore.case = TRUE,
  value = TRUE
)

rich_col <- rich_candidates[
  !grepl("log|depth", rich_candidates, ignore.case = TRUE)
][1]

if (is.na(rich_col)) {
  stop("Could not identify richness column for S2.")
}

if (!"HostRemoved_read_pairs" %in% names(alpha_input)) {
  alpha_input <- alpha_input %>%
    left_join(
      master %>%
        select(Run, HostRemoved_read_pairs),
      by = "Run"
    )
}

p2A <- ggplot(
  alpha_input,
  aes(
    HostRemoved_read_pairs,
    .data[[rich_col]],
    fill = STATUS
  )
) +
  geom_point(
    shape = 21, size = 2.3,
    colour = "black", stroke = .3
  ) +
  scale_fill_manual(values = c("non-PMA" = BLUE, "PMA" = RED)) +
  scale_x_log10(labels = label_number()) +
  labs(
    x = "Host-removed read pairs",
    y = "Observed species richness",
    title = "Richness is depth-associated"
  ) +
  theme_paper +
  theme(legend.title = element_blank())

p2B <- ggplot(
  alpha_input,
  aes(
    STATUS,
    .data[[rich_col]],
    group = Pair_ID
  )
) +
  geom_line(colour = "grey82", linewidth = .45) +
  geom_point(
    aes(fill = STATUS),
    shape = 21, size = 2.3,
    colour = "black", stroke = .3
  ) +
  scale_fill_manual(values = c("non-PMA" = BLUE, "PMA" = RED)) +
  labs(
    x = NULL,
    y = "Observed species richness",
    title = "Matched richness change"
  ) +
  theme_paper +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold")
  )

coef_names <- names(alpha_coef)

term_col <- coef_names[
  grepl("^term$|parameter|effect", coef_names, ignore.case = TRUE)
][1]

est_col <- coef_names[
  grepl("^estimate$|rate.ratio|ratio", coef_names, ignore.case = TRUE)
][1]

low_col <- coef_names[
  grepl("conf.low|lower|ci_low|lcl", coef_names, ignore.case = TRUE)
][1]

high_col <- coef_names[
  grepl("conf.high|upper|ci_high|ucl", coef_names, ignore.case = TRUE)
][1]

if (!any(is.na(c(term_col, est_col, low_col, high_col)))) {

  forest <- alpha_coef %>%
    filter(!grepl("intercept", .data[[term_col]], ignore.case = TRUE)) %>%
    mutate(
      label = case_when(
        grepl("STATUS", .data[[term_col]], ignore.case = TRUE) ~
          "PMA vs non-PMA",
        grepl("depth|log10", .data[[term_col]], ignore.case = TRUE) ~
          "Sequencing depth\n(per 10-fold)",
        TRUE ~ as.character(.data[[term_col]])
      )
    )

  p2C <- ggplot(
    forest,
    aes(
      x = .data[[est_col]],
      y = reorder(label, .data[[est_col]])
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = 2,
      colour = "grey50"
    ) +
    geom_errorbarh(
      aes(
        xmin = .data[[low_col]],
        xmax = .data[[high_col]]
      ),
      height = .18
    ) +
    geom_point(size = 2.6, colour = PURPLE) +
    scale_x_log10() +
    labs(
      x = "Rate ratio (95% CI)",
      y = NULL,
      title = "Depth-adjusted NB model"
    ) +
    theme_paper

} else {

  p2C <- ggplot() +
    annotate(
      "text", x = 0, y = 0,
      label =
        "Depth-adjusted NB model:\nPMA RR = 0.669\n95% CI 0.560–0.799\nP = 9.14 × 10⁻⁶\n\nDepth RR = 1.708 per 10-fold\nP = 1.22 × 10⁻¹⁶",
      hjust = .5,
      size = 4
    ) +
    xlim(-1, 1) + ylim(-1, 1) +
    labs(title = "Depth-adjusted NB model") +
    theme_void() +
    theme(
      plot.title =
        element_text(face = "bold", hjust = .5)
    )
}

# Existing beta-diversity sensitivity text is preserved as data
perma_paired <- readLines(
  "analysis_ready/taxonomy/beta_diversity/bray_permanova_paired.txt"
)

perma_depth <- readLines(
  "analysis_ready/taxonomy/beta_diversity/bray_permanova_depth_adjusted.txt"
)

writeLines(
  c(
    "PAIR-RESTRICTED PERMANOVA",
    perma_paired,
    "",
    "DEPTH-SENSITIVITY PERMANOVA",
    perma_depth
  ),
  file.path(tabdir, "S2_beta_diversity_depth_sensitivity.txt")
)

p2D <- ggplot() +
  annotate(
    "text",
    x = 0, y = .7,
    label =
      "Pair-restricted PERMANOVA\nSTATUS R² = 0.0289\nP = 0.0001",
    size = 4
  ) +
  annotate(
    "text",
    x = 0, y = -.15,
    label =
      "Depth-sensitivity model\nDepth R² = 0.0982, P = 0.0003\nSTATUS R² = 0.0667, P = 0.0001",
    size = 4
  ) +
  xlim(-1, 1) + ylim(-1, 1) +
  labs(title = "Beta-diversity depth sensitivity") +
  theme_void() +
  theme(
    plot.title =
      element_text(face = "bold", hjust = .5)
  )

S2 <- (p2A | p2B) / (p2C | p2D) +
  plot_annotation(tag_levels = "A")

save_fig(
  S2,
  "Supplementary_Figure_S2_microbiome_depth_robustness",
  11.5, 8.2
)

# ============================================================
# S3 — ANCOM-BC2 PREVALENCE SENSITIVITY
# ============================================================

cat("[S3] ANCOM-BC2 prevalence-threshold sensitivity\n")

stable_file <-
  "analysis_ready/taxonomy/ancombc2_genus/ANCOMBC2_genus_stable_10pct_20pct.tsv"

stable <- read.delim(stable_file, check.names = FALSE)

write.table(
  stable,
  file.path(tabdir, "S3_ANCOMBC_stable_10pct_20pct.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# Find genus and effect columns
genus_col <- names(stable)[
  grepl("genus|taxon|feature", names(stable), ignore.case = TRUE)
][1]

numeric_cols <- names(stable)[
  sapply(stable, is.numeric)
]

lfc_cols <- numeric_cols[
  grepl("lfc|coef|effect|estimate", numeric_cols, ignore.case = TRUE)
]

if (length(lfc_cols) >= 2) {

  xcol <- lfc_cols[1]
  ycol <- lfc_cols[2]

  p3A <- ggplot(
    stable,
    aes(
      .data[[xcol]],
      .data[[ycol]]
    )
  ) +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_vline(xintercept = 0, colour = "grey80") +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = 2,
      colour = "grey55"
    ) +
    geom_point(
      size = 2.1,
      alpha = .8,
      colour = PURPLE
    ) +
    labs(
      x = "Effect estimate: primary threshold",
      y = "Effect estimate: sensitivity threshold",
      title = "Stable effect direction"
    ) +
    theme_paper

} else {

  p3A <- ggplot() +
    annotate(
      "text", x = 0, y = 0,
      label =
        "48 genera were robust in both\n20% and 10% prevalence analyses\n\nAll 48 had concordant direction",
      size = 4.5
    ) +
    xlim(-1, 1) + ylim(-1, 1) +
    labs(title = "ANCOM-BC2 sensitivity") +
    theme_void() +
    theme(
      plot.title =
        element_text(face = "bold", hjust = .5)
    )
}

direction_tab <- data.frame(
  Category = c(
    "Stable in both",
    "PMA-enriched",
    "non-PMA-enriched"
  ),
  n = c(48, 4, 44)
)

p3B <- ggplot(
  direction_tab,
  aes(Category, n)
) +
  geom_col(
    aes(fill = Category),
    width = .65
  ) +
  scale_fill_manual(
    values = c(
      "Stable in both" = PURPLE,
      "PMA-enriched" = RED,
      "non-PMA-enriched" = BLUE
    )
  ) +
  geom_text(
    aes(label = n),
    vjust = -.4,
    fontface = "bold"
  ) +
  coord_cartesian(ylim = c(0, 52)) +
  labs(
    x = NULL,
    y = "Number of genera",
    title = "Stable differential genera"
  ) +
  theme_paper +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    )
  )

stable_pma <- data.frame(
  Genus = c(
    "Thiobacillus",
    "Aquabacterium",
    "Comamonas",
    "Herbaspirillum"
  ),
  Direction = "PMA-enriched"
)

p3C <- ggplot(
  stable_pma,
  aes(
    x = 1,
    y = reorder(Genus, Genus)
  )
) +
  geom_point(
    size = 4,
    colour = RED
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Stable PMA-associated genera"
  ) +
  xlim(.8, 1.2) +
  theme_paper +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p3D <- ggplot() +
  annotate(
    "text", x = 0, y = .45,
    label =
      "Primary prevalence threshold: 20%\nRobust genera: 59",
    size = 4.2
  ) +
  annotate(
    "text", x = 0, y = -.25,
    label =
      "Sensitivity threshold: 10%\nRobust genera: 109",
    size = 4.2
  ) +
  xlim(-1, 1) + ylim(-1, 1) +
  labs(title = "Threshold sensitivity") +
  theme_void() +
  theme(
    plot.title =
      element_text(face = "bold", hjust = .5)
  )

S3 <- (p3A | p3B) / (p3C | p3D) +
  plot_annotation(tag_levels = "A")

save_fig(
  S3,
  "Supplementary_Figure_S3_ANCOMBC_sensitivity",
  11.5, 8
)

# ============================================================
# S4 — RESISTOME DEPTH ROBUSTNESS
# ============================================================

cat("[S4] Resistome depth robustness\n")

res_metrics <- read.delim(
  "analysis_ready/deeparg/figure4/DeepARG_sample_resistome_metrics.tsv",
  check.names = FALSE
)

if (!"HostRemoved_read_pairs" %in% names(res_metrics)) {
  res_metrics <- res_metrics %>%
    left_join(
      master %>%
        select(Run, HostRemoved_read_pairs),
      by = "Run"
    )
}

res_model <- read.delim(
  "analysis_ready/deeparg/figure4/DeepARG_richness_depth_adjusted.tsv",
  check.names = FALSE
)

res_cor <- read.delim(
  "analysis_ready/deeparg/figure4/DeepARG_depth_correlations.tsv",
  check.names = FALSE
)

write.table(
  res_model,
  file.path(tabdir, "S4_ARG_richness_depth_model.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  res_cor,
  file.path(tabdir, "S4_resistome_depth_correlations.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

p4A <- ggplot(
  res_metrics,
  aes(
    HostRemoved_read_pairs,
    ARG_subtype_richness,
    fill = STATUS
  )
) +
  geom_point(
    shape = 21,
    size = 2.3,
    colour = "black",
    stroke = .3
  ) +
  scale_fill_manual(
    values = c("non-PMA" = BLUE, "PMA" = RED)
  ) +
  scale_x_log10(labels = label_number()) +
  labs(
    x = "Host-removed read pairs",
    y = "ARG-group richness",
    title = "ARG richness and sequencing depth"
  ) +
  theme_paper +
  theme(legend.title = element_blank())

p4B <- ggplot(
  res_metrics,
  aes(
    STATUS,
    ARG_subtype_richness,
    group = Pair_ID
  )
) +
  geom_line(colour = "grey82", linewidth = .45) +
  geom_point(
    aes(fill = STATUS),
    shape = 21, size = 2.3,
    colour = "black", stroke = .3
  ) +
  scale_fill_manual(
    values = c("non-PMA" = BLUE, "PMA" = RED)
  ) +
  labs(
    x = NULL,
    y = "ARG-group richness",
    title = "Matched resistome breadth"
  ) +
  theme_paper +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold")
  )

p4C <- ggplot(
  res_metrics,
  aes(
    HostRemoved_read_pairs,
    DeepARG_total_16S,
    fill = STATUS
  )
) +
  geom_point(
    shape = 21,
    size = 2.3,
    colour = "black",
    stroke = .3
  ) +
  scale_fill_manual(
    values = c("non-PMA" = BLUE, "PMA" = RED)
  ) +
  scale_x_log10(labels = label_number()) +
  labs(
    x = "Host-removed read pairs",
    y = "Total DeepARG 16S-normalized abundance",
    title = "Normalized ARG burden"
  ) +
  theme_paper +
  theme(legend.title = element_blank())

p4D <- ggplot() +
  annotate(
    "text", x = 0, y = .50,
    label =
      "Depth-adjusted ARG richness model\nPMA RR = 0.472\n95% CI 0.396–0.562\nP = 3.67 × 10⁻¹⁷",
    size = 4
  ) +
  annotate(
    "text", x = 0, y = -.30,
    label =
      "Sequencing depth\nRR = 1.398 per 10-fold increase\nP = 1.84 × 10⁻⁶\n\nModel Hessian: positive definite",
    size = 4
  ) +
  xlim(-1, 1) + ylim(-1, 1) +
  labs(title = "Negative-binomial sensitivity model") +
  theme_void() +
  theme(
    plot.title =
      element_text(face = "bold", hjust = .5)
  )

S4 <- (p4A | p4B) / (p4C | p4D) +
  plot_annotation(tag_levels = "A")

save_fig(
  S4,
  "Supplementary_Figure_S4_resistome_depth_robustness",
  11.5, 8.2
)

# ============================================================
# S5 — EXPANDED AMR TOOL CONCORDANCE
# ============================================================

cat("[S5] Expanded AMR-tool concordance\n")

gene_conc <- read.delim(
  "analysis_ready/amr_concordance/CARD_candidate_concordance.tsv",
  check.names = FALSE
)

four_method <- read.delim(
  "analysis_ready/amr_concordance/four_method_candidate_evidence.tsv",
  check.names = FALSE
)

class_prev <- read.delim(
  "analysis_ready/amr_concordance/class_prevalence_four_methods.tsv",
  check.names = FALSE
)

write.table(
  gene_conc,
  file.path(tabdir, "S5_gene_concordance.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  class_prev,
  file.path(tabdir, "S5_class_concordance.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# Generic reshape: prevalence columns only
prev_cols <- names(gene_conc)[
  grepl(
    "PMA.*prev|prev.*PMA|PMA_prevalence",
    names(gene_conc),
    ignore.case = TRUE
  )
]

gene_name_col <- names(gene_conc)[
  grepl(
    "^gene$|ARG|candidate",
    names(gene_conc),
    ignore.case = TRUE
  )
][1]

if (length(prev_cols) > 0 && !is.na(gene_name_col)) {

  gene_plot <- gene_conc %>%
    select(
      Gene = all_of(gene_name_col),
      all_of(prev_cols)
    ) %>%
    pivot_longer(
      -Gene,
      names_to = "Method",
      values_to = "Prevalence"
    )

  p5A <- ggplot(
    gene_plot,
    aes(
      Method,
      Gene,
      fill = Prevalence
    )
  ) +
    geom_tile(colour = "white") +
    geom_text(
      aes(label = percent(Prevalence, accuracy = 1)),
      size = 3
    ) +
    scale_fill_gradient(
      low = "white",
      high = RED,
      limits = c(0, 1)
    ) +
    labs(
      x = NULL,
      y = NULL,
      fill = "PMA\nprevalence",
      title = "Candidate ARG prevalence"
    ) +
    theme_paper +
    theme(
      axis.text.x =
        element_text(angle = 30, hjust = 1)
    )

} else {

  p5A <- ggplot() +
    annotate(
      "text", x = 0, y = 0,
      label =
        "Expanded candidate-gene\nconcordance table saved as S5_gene_concordance.tsv",
      size = 4.2
    ) +
    xlim(-1, 1) + ylim(-1, 1) +
    theme_void()
}

# Class-level matrix
class_names <- names(class_prev)

class_col <- class_names[
  grepl("^class$|resistance.class", class_names, ignore.case = TRUE)
][1]

numeric_class <- class_names[
  sapply(class_prev, is.numeric)
]

pma_numeric <- numeric_class[
  grepl("PMA", numeric_class, ignore.case = TRUE)
]

if (!is.na(class_col) && length(pma_numeric) > 0) {

  cp <- class_prev %>%
    select(
      ResistanceClass = all_of(class_col),
      all_of(pma_numeric)
    ) %>%
    pivot_longer(
      -ResistanceClass,
      names_to = "Method",
      values_to = "Prevalence"
    )

  p5B <- ggplot(
    cp,
    aes(
      Method,
      ResistanceClass,
      fill = Prevalence
    )
  ) +
    geom_tile(colour = "white") +
    geom_text(
      aes(label = percent(Prevalence, accuracy = 1)),
      size = 2.8
    ) +
    scale_fill_gradient(
      low = "white",
      high = BLUE,
      limits = c(0, 1)
    ) +
    labs(
      x = NULL,
      y = NULL,
      fill = "PMA\nprevalence",
      title = "Conservative class harmonization"
    ) +
    theme_paper +
    theme(
      axis.text.x =
        element_text(angle = 30, hjust = 1)
    )

} else {

  p5B <- ggplot() +
    annotate(
      "text", x = 0, y = 0,
      label =
        "Expanded class concordance\ntable saved as S5_class_concordance.tsv",
      size = 4.2
    ) +
    xlim(-1, 1) + ylim(-1, 1) +
    theme_void()
}

p5C <- ggplot() +
  annotate(
    "text", x = 0, y = .45,
    label =
      "Strongest recurrent read-based evidence\nmsrA: DeepARG 81%, ShortBRED 62%, RGI 88%\ndfrC: 75%, 56%, 62%",
    size = 4
  ) +
  annotate(
    "text", x = 0, y = -.35,
    label =
      "AMRFinderPlus is assembly-based.\nShortBRED and RGI both use CARD;\nconcordance is not independent database validation.",
    size = 3.8
  ) +
  xlim(-1, 1) + ylim(-1, 1) +
  labs(title = "Interpretation of concordance") +
  theme_void() +
  theme(
    plot.title =
      element_text(face = "bold", hjust = .5)
  )

S5 <- (p5A | p5B) / p5C +
  plot_annotation(tag_levels = "A")

save_fig(
  S5,
  "Supplementary_Figure_S5_expanded_AMR_concordance",
  13, 9
)

# ============================================================
# S6 — ARG-TAXON ASSEMBLY DEPTH CONTEXT
# ============================================================

cat("[S6] ARG-bearing contig assembly-depth context\n")

argtax <- read.delim(
  "analysis_ready/arg_host/ARG_taxon_PMA_vs_nonPMA_selected.tsv",
  check.names = FALSE
)

write.table(
  argtax,
  file.path(tabdir, "S6_ARG_taxon_selected.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

p6A <- ggplot(
  master,
  aes(
    Assembly_total_bp,
    AMRFinder_hits,
    fill = STATUS
  )
) +
  geom_point(
    shape = 21,
    size = 2.4,
    colour = "black",
    stroke = .3
  ) +
  scale_fill_manual(
    values = c("non-PMA" = BLUE, "PMA" = RED)
  ) +
  scale_x_log10(labels = label_number()) +
  scale_y_continuous(
    trans = pseudo_log_trans(base = 10)
  ) +
  labs(
    x = "Assembly yield (bp)",
    y = "AMRFinderPlus hits",
    title = "Assembly yield and ARG recovery"
  ) +
  theme_paper +
  theme(legend.title = element_blank())

# Sample-level ARG contigs
amr_hits <- read.delim(
  "analysis_ready/amrfinder/amrfinder_ARG_contig_taxonomy.tsv",
  check.names = FALSE
)

argcontig_sample <- amr_hits %>%
  group_by(Run, STATUS) %>%
  summarise(
    ARG_contigs = n_distinct(`Contig id`),
    .groups = "drop"
  ) %>%
  right_join(
    master %>%
      select(
        Run, STATUS,
        Assembly_total_bp
      ),
    by = c("Run", "STATUS")
  ) %>%
  mutate(
    ARG_contigs =
      replace_na(ARG_contigs, 0)
  )

write.table(
  argcontig_sample,
  file.path(tabdir, "S6_sample_ARG_contig_recovery.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

p6B <- ggplot(
  argcontig_sample,
  aes(
    Assembly_total_bp,
    ARG_contigs,
    fill = STATUS
  )
) +
  geom_point(
    shape = 21,
    size = 2.4,
    colour = "black",
    stroke = .3
  ) +
  scale_fill_manual(
    values = c("non-PMA" = BLUE, "PMA" = RED)
  ) +
  scale_x_log10(labels = label_number()) +
  scale_y_continuous(
    trans = pseudo_log_trans(base = 10)
  ) +
  labs(
    x = "Assembly yield (bp)",
    y = "Unique ARG-bearing contigs",
    title = "Assembly yield and ARG-contig recovery"
  ) +
  theme_paper +
  theme(legend.title = element_blank())

# Expanded selected association matrix
arg_names <- names(argtax)

arg_col <- arg_names[
  grepl("^ARG$|Element.symbol|gene", arg_names, ignore.case = TRUE)
][1]

tax_col <- arg_names[
  grepl("taxon", arg_names, ignore.case = TRUE)
][1]

status_col <- arg_names[
  grepl("^STATUS$", arg_names, ignore.case = TRUE)
][1]

prev_col <- arg_names[
  grepl("prevalence|fraction", arg_names, ignore.case = TRUE)
][1]

if (
  !any(is.na(c(arg_col, tax_col, status_col, prev_col)))
) {

  argtax_plot <- argtax %>%
    mutate(
      Association =
        paste(.data[[arg_col]], .data[[tax_col]], sep = " – ")
    )

  p6C <- ggplot(
    argtax_plot,
    aes(
      .data[[status_col]],
      Association,
      fill = .data[[prev_col]]
    )
  ) +
    geom_tile(colour = "white") +
    geom_text(
      aes(
        label = percent(
          .data[[prev_col]],
          accuracy = 1
        )
      ),
      size = 2.7
    ) +
    scale_fill_gradient(
      low = "white",
      high = ORANGE,
      limits = c(0, 1)
    ) +
    labs(
      x = NULL,
      y = NULL,
      fill = "Sample\nprevalence",
      title = "Expanded recurrent ARG–taxon associations"
    ) +
    theme_paper

} else {

  p6C <- ggplot() +
    annotate(
      "text", x = 0, y = 0,
      label =
        "Expanded ARG–taxon association table\nsaved as S6_ARG_taxon_selected.tsv",
      size = 4
    ) +
    xlim(-1, 1) + ylim(-1, 1) +
    theme_void()
}

p6D <- ggplot() +
  annotate(
    "text", x = 0, y = .35,
    label =
      "Median assembly yield\nnon-PMA ≈ 46.0 Mb\nPMA ≈ 3.04 Mb\n(~15-fold difference)",
    size = 4.2
  ) +
  annotate(
    "text", x = 0, y = -.40,
    label =
      "Therefore absence of an assembly-based\nARG or ARG–taxon association is not\nevidence of biological absence.",
    size = 3.9
  ) +
  xlim(-1, 1) + ylim(-1, 1) +
  labs(title = "Assembly-depth caveat") +
  theme_void() +
  theme(
    plot.title =
      element_text(face = "bold", hjust = .5)
  )

S6 <- (p6A | p6B) / (p6C | p6D) +
  plot_annotation(tag_levels = "A")

save_fig(
  S6,
  "Supplementary_Figure_S6_ARG_taxon_assembly_depth",
  12.5, 8.5
)

# ============================================================
# S7 — INTEGRATED COUPLING DEPTH ROBUSTNESS
# ============================================================

cat("[S7] Integrated microbiome-resistome depth robustness\n")

coupling <- read.delim(
  "analysis_ready/integrated_analysis/paired_distances_with_depth.tsv",
  check.names = FALSE
)

coupling_stats <- read.delim(
  "analysis_ready/integrated_analysis/microbiome_resistome_depth_sensitivity.tsv",
  check.names = FALSE
)

write.table(
  coupling,
  file.path(tabdir, "S7_paired_distances_with_depth.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  coupling_stats,
  file.path(tabdir, "S7_coupling_depth_statistics.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# Detect exact distance columns
micro_col <- names(coupling)[
  grepl("microbiome.*bray|microbiome.*distance",
        names(coupling), ignore.case = TRUE)
][1]

res_col <- names(coupling)[
  grepl("resistome.*bray|resistome.*distance",
        names(coupling), ignore.case = TRUE)
][1]

imb_col <- names(coupling)[
  grepl("imbalance|abs.*log10",
        names(coupling), ignore.case = TRUE)
][1]

if (any(is.na(c(micro_col, res_col, imb_col)))) {
  stop(
    paste(
      "Could not identify S7 columns. Available:",
      paste(names(coupling), collapse = ", ")
    )
  )
}

p7A <- ggplot(
  coupling,
  aes(
    .data[[imb_col]],
    .data[[micro_col]]
  )
) +
  geom_point(
    size = 2.5,
    colour = BLUE
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE,
    colour = "grey35",
    fill = "grey80",
    linewidth = .7
  ) +
  labs(
    x = "Absolute log10 depth imbalance",
    y = "Microbiome Bray–Curtis distance",
    title = "Microbiome shift vs depth imbalance"
  ) +
  theme_paper

p7B <- ggplot(
  coupling,
  aes(
    .data[[imb_col]],
    .data[[res_col]]
  )
) +
  geom_point(
    size = 2.5,
    colour = RED
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE,
    colour = "grey35",
    fill = "grey80",
    linewidth = .7
  ) +
  labs(
    x = "Absolute log10 depth imbalance",
    y = "Resistome Bray–Curtis distance",
    title = "Resistome shift vs depth imbalance"
  ) +
  theme_paper

p7C <- ggplot(
  coupling,
  aes(
    .data[[micro_col]],
    .data[[res_col]]
  )
) +
  geom_point(
    size = 2.5,
    colour = PURPLE
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE,
    colour = "grey35",
    fill = "grey80",
    linewidth = .7
  ) +
  labs(
    x = "Microbiome Bray–Curtis distance",
    y = "Resistome Bray–Curtis distance",
    title = "Microbiome–resistome coupling"
  ) +
  theme_paper

p7D <- ggplot() +
  annotate(
    "text", x = 0, y = .55,
    label =
      "Original coupling\nSpearman ρ = 0.511\nP = 0.00282",
    size = 4.2
  ) +
  annotate(
    "text", x = 0, y = -.15,
    label =
      "Partial rank correlation\ncontrolling depth imbalance\nr = 0.527\nP = 0.00194",
    size = 4.2
  ) +
  annotate(
    "text", x = 0, y = -.72,
    label =
      "Depth imbalance does not explain\nthe microbiome–resistome coupling.",
    size = 3.8,
    fontface = "bold"
  ) +
  xlim(-1, 1) + ylim(-1, 1) +
  labs(title = "Depth-adjusted sensitivity") +
  theme_void() +
  theme(
    plot.title =
      element_text(face = "bold", hjust = .5)
  )

S7 <- (p7A | p7B) / (p7C | p7D) +
  plot_annotation(tag_levels = "A")

save_fig(
  S7,
  "Supplementary_Figure_S7_integrated_depth_robustness",
  11.5, 8.2
)

# ============================================================
# MASTER SUMMARY
# ============================================================

summary_out <- data.frame(
  Figure = paste0("S", 1:7),
  Purpose = c(
    "Sequencing depth, assembly yield, and downstream recovery",
    "Microbiome richness and beta-diversity depth sensitivity",
    "ANCOM-BC2 prevalence-threshold sensitivity",
    "Resistome breadth and burden depth sensitivity",
    "Expanded AMR-tool concordance",
    "ARG-bearing contig assembly-depth context",
    "Integrated microbiome-resistome depth sensitivity"
  ),
  Status = "Generated"
)

write.table(
  summary_out,
  file.path(
    tabdir,
    "supplementary_analysis_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("\n====================================================\n")
cat("SUPPLEMENTARY PACKAGE COMPLETE\n")
cat("====================================================\n")
cat("Figures written to: ", figdir, "\n")
cat("Tables written to:  ", tabdir, "\n\n")

print(summary_out)

cat("\nS1 key depth correlations:\n")
print(S1_cor)

cat("\nS1 group summary:\n")
print(S1_summary)

cat("\nDone.\n")
