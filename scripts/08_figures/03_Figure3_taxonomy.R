
# ============================================================
# Figure 3: Taxonomic composition and differential abundance
# NASA MT2 PMA vs non-PMA metagenomics
#
# Fig 3A: Top-genus relative abundance
# Fig 3B: Heatmap of stable differential genera
# Fig 3C: ANCOM-BC2 forest plot
#
# Primary inference:
#   ANCOM-BC2, >=20% prevalence, Pair_ID random intercept
# Sensitivity:
#   >=10% prevalence
# Stable taxa:
#   robust in BOTH analyses with concordant direction
# ============================================================

.libPaths(c(path.expand("~/R/4.5"), .libPaths()))

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(pheatmap)
    library(grid)
})

# ------------------------------------------------------------
# Colors
# ------------------------------------------------------------

COL_NONPMA <- "#1479E8"
COL_PMA    <- "#FF4B4B"
COL_PAIR   <- "grey80"

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

genus_ra_file <-
    "analysis_ready/taxonomy/bracken_genus_relative_abundance.tsv"

meta_file <-
    "analysis_ready/metadata/master_sample_table.tsv"

primary_file <-
    "analysis_ready/taxonomy/ancombc2_genus/ANCOMBC2_genus_20pct_primary_STATUS.tsv"

stable_file <-
    "analysis_ready/taxonomy/ancombc2_genus/ANCOMBC2_genus_stable_10pct_20pct.tsv"

outdir <- "figures/figure3"
tabdir <- "analysis_ready/taxonomy/figure3"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Read data
# ------------------------------------------------------------

ra <- read.delim(
    genus_ra_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

meta <- read.delim(
    meta_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

primary <- read.delim(
    primary_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

stable <- read.delim(
    stable_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

meta$STATUS <- factor(
    meta$STATUS,
    levels = c("non-PMA", "PMA")
)

# Genus matrix:
# taxonomy_id | name | Run1 | Run2 ...
genus_names <- ra$name

X <- as.matrix(
    ra[, setdiff(colnames(ra), c("taxonomy_id", "name"))]
)

rownames(X) <- genus_names

storage.mode(X) <- "numeric"

# Ensure sample order agrees with metadata
common_runs <- intersect(colnames(X), meta$Run)

X <- X[, common_runs, drop = FALSE]

meta <- meta[
    match(common_runs, meta$Run),
    ,
    drop = FALSE
]

stopifnot(all(colnames(X) == meta$Run))

cat("\n===== FIGURE 3 INPUT =====\n")
cat("Samples:", ncol(X), "\n")
cat("Genera:", nrow(X), "\n")
cat("Stable taxa:", nrow(stable), "\n")
cat("PMA:", sum(meta$STATUS == "PMA"), "\n")
cat("non-PMA:", sum(meta$STATUS == "non-PMA"), "\n")

# ============================================================
# FIGURE 3A
# Top-genus composition
# ============================================================

mean_abundance <- rowMeans(X, na.rm = TRUE)

TOP_N <- 12

top_genera <- names(
    sort(mean_abundance, decreasing = TRUE)
)[seq_len(min(TOP_N, length(mean_abundance)))]

comp <- as.data.frame(t(X[top_genera, , drop = FALSE]))

comp$Run <- rownames(comp)

comp <- comp %>%
    left_join(
        meta %>%
            select(Run, STATUS, Pair_ID),
        by = "Run"
    )

comp_long <- comp %>%
    pivot_longer(
        cols = all_of(top_genera),
        names_to = "Genus",
        values_to = "Abundance"
    )

# Relative-abundance tables may be fraction or percentage.
# Detect scale.
if (max(comp_long$Abundance, na.rm = TRUE) <= 1.01) {
    comp_long$Abundance_pct <- comp_long$Abundance * 100
} else {
    comp_long$Abundance_pct <- comp_long$Abundance
}

# Average composition by treatment.
comp_summary <- comp_long %>%
    group_by(STATUS, Genus) %>%
    summarise(
        Mean_abundance = mean(Abundance_pct, na.rm = TRUE),
        .groups = "drop"
    )

# Add "Other" so bars sum to ~100%.
other <- comp_summary %>%
    group_by(STATUS) %>%
    summarise(
        Mean_abundance =
            max(0, 100 - sum(Mean_abundance)),
        .groups = "drop"
    ) %>%
    mutate(Genus = "Other")

comp_summary <- bind_rows(comp_summary, other)

# Order genera by overall abundance.
genus_order <- comp_summary %>%
    filter(Genus != "Other") %>%
    group_by(Genus) %>%
    summarise(x = mean(Mean_abundance), .groups = "drop") %>%
    arrange(desc(x)) %>%
    pull(Genus)

comp_summary$Genus <- factor(
    comp_summary$Genus,
    levels = c(genus_order, "Other")
)

write.table(
    comp_summary,
    file.path(tabdir, "Fig3A_top_genus_composition.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

pA <- ggplot(
    comp_summary,
    aes(
        x = STATUS,
        y = Mean_abundance,
        fill = Genus
    )
) +
    geom_col(
        width = 0.72,
        colour = "white",
        linewidth = 0.15
    ) +
    labs(
        x = NULL,
        y = "Mean relative abundance (%)",
        fill = "Genus"
    ) +
    theme_classic(base_size = 11) +
    theme(
        axis.text.x = element_text(
            face = "bold",
            colour = c(COL_NONPMA, COL_PMA)
        ),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(
            face = "italic",
            size = 8
        ),
        plot.margin = margin(6, 6, 6, 6)
    )

ggsave(
    file.path(outdir, "Fig3A_TopGenusComposition.pdf"),
    pA,
    width = 5.2,
    height = 5.0
)

ggsave(
    file.path(outdir, "Fig3A_TopGenusComposition.png"),
    pA,
    width = 5.2,
    height = 5.0,
    dpi = 600
)

# ============================================================
# FIGURE 3B
# Heatmap of stable differential genera
# ============================================================

# Use:
#   top 8 stable non-PMA taxa by primary FDR
#   all 4 stable PMA taxa
#
# This avoids subjective selection and gives balanced biology.

stable2 <- stable %>%
    filter(same_direction %in% TRUE)

nonpma_heat <- stable2 %>%
    filter(lfc20 < 0) %>%
    arrange(q20) %>%
    slice_head(n = 8)

pma_heat <- stable2 %>%
    filter(lfc20 > 0) %>%
    arrange(q20)

heat_taxa <- unique(
    c(nonpma_heat$genus, pma_heat$genus)
)

heat_taxa <- intersect(heat_taxa, rownames(X))

H <- X[heat_taxa, , drop = FALSE]

# log transform relative abundance
Hlog <- log10(H + 1e-6)

# Row Z-score to visualize within-genus pattern
Hz <- t(scale(t(Hlog)))

Hz[!is.finite(Hz)] <- 0

# Arrange samples by Pair_ID and treatment.
sample_order <- meta %>%
  arrange(STATUS, Pair_ID) %>%
  pull(Run)

Hz <- Hz[, sample_order, drop = FALSE]

ann <- meta %>%
    select(Run, STATUS, Pair_ID)

rownames(ann) <- ann$Run

ann <- ann[
    sample_order,
    c("STATUS"),
    drop = FALSE
]

ann_colors <- list(
    STATUS = c(
        "non-PMA" = COL_NONPMA,
        "PMA" = COL_PMA
    )
)

# Keep genus names clean in the final heatmap.
# Differential direction is communicated in Fig. 3C.
rownames(Hz) <- rownames(Hz)

write.table(
    Hz,
    file.path(tabdir, "Fig3B_heatmap_Zscores.tsv"),
    sep = "\t",
    quote = FALSE,
    col.names = NA
)

pdf(
    file.path(outdir, "Fig3B_StableGeneraHeatmap.pdf"),
    width = 10.5,
    height = 5.8
)

pheatmap(
    Hz,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = ann,
    annotation_colors = ann_colors,
    show_colnames = FALSE,
    fontsize_row = 9,
    border_color = NA,
    color = colorRampPalette(
        c("#2166AC", "white", "#B2182B")
    )(100),
    main = "Stable differential genera"
)

dev.off()

png(
    file.path(outdir, "Fig3B_StableGeneraHeatmap.png"),
    width = 10.5,
    height = 5.8,
    units = "in",
    res = 600
)

pheatmap(
    Hz,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = ann,
    annotation_colors = ann_colors,
    show_colnames = FALSE,
    fontsize_row = 9,
    border_color = NA,
    color = colorRampPalette(
        c("#2166AC", "white", "#B2182B")
    )(100),
    main = "Stable differential genera"
)

dev.off()

# ============================================================
# FIGURE 3C
# ANCOM-BC2 forest plot
# ============================================================

# Reproducible main-figure selection:
# top 12 stable non-PMA taxa by primary q-value
# plus ALL stable PMA taxa.

nonpma_forest <- stable2 %>%
    filter(lfc20 < 0) %>%
    arrange(q20) %>%
    slice_head(n = 12)

pma_forest <- stable2 %>%
    filter(lfc20 > 0) %>%
    arrange(q20)

forest_taxa <- unique(
    c(nonpma_forest$genus, pma_forest$genus)
)

forest <- primary %>%
    filter(genus %in% forest_taxa) %>%
    mutate(
        Direction = ifelse(
            lfc_STATUSPMA > 0,
            "Higher in PMA",
            "Higher in non-PMA"
        ),
        CI_low =
            lfc_STATUSPMA -
            1.96 * se_STATUSPMA,
        CI_high =
            lfc_STATUSPMA +
            1.96 * se_STATUSPMA,
        FDR_label = case_when(
            q_STATUSPMA < 0.001 ~ "FDR < 0.001",
            q_STATUSPMA < 0.01  ~ "FDR < 0.01",
            q_STATUSPMA < 0.05  ~ "FDR < 0.05",
            TRUE ~ "NS"
        )
    ) %>%
    arrange(lfc_STATUSPMA)

forest$genus <- factor(
    forest$genus,
    levels = forest$genus
)

write.table(
    forest,
    file.path(tabdir, "Fig3C_ANCOMBC2_forest_data.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

pC <- ggplot(
    forest,
    aes(
        x = lfc_STATUSPMA,
        y = genus,
        colour = Direction
    )
) +
    geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.5,
        colour = "grey45"
    ) +
    geom_errorbarh(
        aes(
            xmin = CI_low,
            xmax = CI_high
        ),
        height = 0,
        linewidth = 0.65
    ) +
    geom_point(size = 2.7) +
    scale_colour_manual(
        values = c(
            "Higher in non-PMA" = COL_NONPMA,
            "Higher in PMA" = COL_PMA
        )
    ) +
    labs(
        x = "ANCOM-BC2 log-fold change\n(PMA vs non-PMA)",
        y = NULL,
        colour = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
        axis.text.y = element_text(
            face = "italic",
            colour = "black"
        ),
        legend.position = "top",
        plot.margin = margin(6, 10, 6, 6)
    )

ggsave(
    file.path(outdir, "Fig3C_ANCOMBC2_Forest.pdf"),
    pC,
    width = 6.2,
    height = 5.5
)

ggsave(
    file.path(outdir, "Fig3C_ANCOMBC2_Forest.png"),
    pC,
    width = 6.2,
    height = 5.5,
    dpi = 600
)

# ============================================================
# Save selected taxa
# ============================================================

selected <- stable2 %>%
    filter(genus %in% unique(c(heat_taxa, forest_taxa))) %>%
    arrange(lfc20)

write.table(
    selected,
    file.path(
        tabdir,
        "Figure3_selected_stable_genera.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# Summary
# ============================================================

cat("\n===== FIGURE 3 COMPLETE =====\n")

cat(
    "Stable genera in both prevalence analyses:",
    nrow(stable2),
    "\n"
)

cat(
    "Stable higher in non-PMA:",
    sum(stable2$lfc20 < 0),
    "\n"
)

cat(
    "Stable higher in PMA:",
    sum(stable2$lfc20 > 0),
    "\n"
)

cat(
    "Heatmap taxa:",
    length(heat_taxa),
    "\n"
)

cat(
    "Forest taxa:",
    length(forest_taxa),
    "\n"
)

cat("\nFiles written to:\n")
cat(outdir, "\n")
cat(tabdir, "\n")

