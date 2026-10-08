#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vegan)
  library(grid)
  library(ggplot2)
  library(dplyr)
  library(permute)
})

set.seed(20261005)

ROOT <- getwd()

REL_FILE <- file.path(
  ROOT,
  "analysis_ready/taxonomy/bracken_species_relative_abundance.tsv"
)

MASTER_FILE <- file.path(
  ROOT,
  "analysis_ready/metadata/master_sample_table.tsv"
)

OUTDIR <- file.path(
  ROOT,
  "analysis_ready/taxonomy/beta_diversity"
)

FIGDIR <- file.path(
  ROOT,
  "figures/figure2"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

COL_NONPMA <- "#1479E8"
COL_PMA    <- "#FF4B4B"
COL_PAIR   <- "grey75"

cat("========================================\n")
cat("BETA DIVERSITY ANALYSIS\n")
cat("========================================\n")

# ============================================================
# Load Bracken relative abundance
# ============================================================

raw <- read.delim(
  REL_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Raw relative-abundance dimensions:",
    nrow(raw), "x", ncol(raw), "\n")

# Matrix should have taxonomy_id + name + 64 samples
id_cols <- intersect(
  c("taxonomy_id", "name"),
  colnames(raw)
)

sample_cols <- setdiff(
  colnames(raw),
  id_cols
)

comm <- t(
  as.matrix(
    raw[, sample_cols, drop = FALSE]
  )
)

storage.mode(comm) <- "numeric"

cat("Community matrix:",
    nrow(comm), "samples x",
    ncol(comm), "species\n")

if (anyNA(comm)) {
  stop("NA values detected in community matrix.")
}

if (any(comm < 0)) {
  stop("Negative abundance values detected.")
}

# ============================================================
# Metadata
# ============================================================

meta <- read.delim(
  MASTER_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

meta <- meta[
  match(rownames(comm), meta$Run),
]

if (anyNA(meta$Run)) {
  stop("Some community samples are missing metadata.")
}

if (!all(meta$Run == rownames(comm))) {
  stop("Community matrix and metadata are not aligned.")
}

meta$STATUS <- factor(
  meta$STATUS,
  levels = c("non-PMA", "PMA")
)

meta$Pair_ID <- factor(meta$Pair_ID)

cat("Samples:", nrow(meta), "\n")
cat("Pairs:", nlevels(meta$Pair_ID), "\n")
cat("PMA:", sum(meta$STATUS == "PMA"), "\n")
cat("non-PMA:", sum(meta$STATUS == "non-PMA"), "\n")

pair_check <- table(meta$Pair_ID, meta$STATUS)

if (!all(pair_check[, "PMA"] == 1) ||
    !all(pair_check[, "non-PMA"] == 1)) {
  stop("Pairing structure is not exactly one PMA + one non-PMA.")
}

# ============================================================
# Remove species absent from all samples
# ============================================================

comm <- comm[, colSums(comm) > 0, drop = FALSE]

cat("Species retained:", ncol(comm), "\n")

# ============================================================
# Bray-Curtis
# ============================================================

bray <- vegdist(
  comm,
  method = "bray"
)

# ============================================================
# PCoA
# ============================================================

pcoa <- cmdscale(
  bray,
  eig = TRUE,
  k = 2,
  add = TRUE
)

coords <- as.data.frame(pcoa$points)

colnames(coords) <- c("PCoA1", "PCoA2")

coords$Run <- rownames(coords)

coords <- coords %>%
  left_join(
    meta %>%
      select(
        Run,
        SampleName,
        Pair_ID,
        STATUS,
        Flight,
        LocationCode,
        Location,
        HostRemoved_read_pairs
      ),
    by = "Run"
  )

positive_eigs <- pcoa$eig[pcoa$eig > 0]

var1 <- 100 * positive_eigs[1] /
  sum(positive_eigs)

var2 <- 100 * positive_eigs[2] /
  sum(positive_eigs)

cat("\n========================================\n")
cat("PCoA\n")
cat("========================================\n")

cat("PCoA1 variance:", round(var1, 2), "%\n")
cat("PCoA2 variance:", round(var2, 2), "%\n")

write.table(
  coords,
  file.path(OUTDIR, "bray_pcoa_coordinates.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ============================================================
# Restricted permutation design
# ============================================================

perm_control <- how(
  nperm = 9999,
  blocks = meta$Pair_ID
)

cat("\nPermutation design:\n")
print(perm_control)

# ============================================================
# Primary paired PERMANOVA
# ============================================================

permanova <- adonis2(
  bray ~ STATUS,
  data = meta,
  permutations = perm_control,
  by = "margin"
)

cat("\n========================================\n")
cat("PAIRED / RESTRICTED PERMANOVA\n")
cat("Bray-Curtis ~ STATUS\n")
cat("Permutations restricted within Pair_ID\n")
cat("========================================\n")

print(permanova)

capture.output(
  permanova,
  file = file.path(
    OUTDIR,
    "bray_permanova_paired.txt"
  )
)

# ============================================================
# Depth-adjusted PERMANOVA sensitivity analysis
#
# We include log10 host-removed read depth and test STATUS
# while preserving within-pair permutations.
# ============================================================

meta$log10_depth <- log10(
  meta$HostRemoved_read_pairs
)

permanova_depth <- adonis2(
  bray ~ log10_depth + STATUS,
  data = meta,
  permutations = perm_control,
  by = "margin"
)

cat("\n========================================\n")
cat("DEPTH-ADJUSTED PAIRED PERMANOVA\n")
cat("Bray-Curtis ~ log10(depth) + STATUS\n")
cat("========================================\n")

print(permanova_depth)

capture.output(
  permanova_depth,
  file = file.path(
    OUTDIR,
    "bray_permanova_depth_adjusted.txt"
  )
)

# ============================================================
# PERMDISP
# ============================================================

bd <- betadisper(
  bray,
  group = meta$STATUS,
  type = "median"
)

bd_anova <- anova(bd)

# PERMDISP is a diagnostic of group dispersion, not the paired
# location test used by PERMANOVA. Use unrestricted permutations.
# Pair_ID restriction remains appropriate for the primary PERMANOVA.
set.seed(20261005)
bd_perm <- permutest(
  bd,
  permutations = 9999,
  pairwise = FALSE
)

cat("\n========================================\n")
cat("PERMDISP\n")
cat("========================================\n")

cat("\nParametric ANOVA:\n")
print(bd_anova)

cat("\nRestricted permutation test:\n")
print(bd_perm)

capture.output(
  bd_anova,
  file = file.path(
    OUTDIR,
    "bray_permdisp_anova.txt"
  )
)

capture.output(
  bd_perm,
  file = file.path(
    OUTDIR,
    "bray_permdisp_restricted.txt"
  )
)

# ============================================================
# Distances within each matched pair
# ============================================================

bray_matrix <- as.matrix(bray)

pair_distances <- lapply(
  levels(meta$Pair_ID),
  function(pair) {

    idx <- which(meta$Pair_ID == pair)

    data.frame(
      Pair_ID = pair,
      Run1 = meta$Run[idx[1]],
      Run2 = meta$Run[idx[2]],
      Bray_distance =
        bray_matrix[idx[1], idx[2]]
    )
  }
)

pair_distances <- bind_rows(pair_distances)

cat("\n========================================\n")
cat("WITHIN-PAIR BRAY-CURTIS DISTANCES\n")
cat("========================================\n")

print(summary(pair_distances$Bray_distance))

write.table(
  pair_distances,
  file.path(
    OUTDIR,
    "within_pair_bray_distances.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ============================================================
# Figure 2B
# ============================================================

p <- ggplot(
  coords,
  aes(
    x = PCoA1,
    y = PCoA2
  )
) +

  geom_segment(
    data = coords,
    aes(
      group = Pair_ID,
      xend = ave(
        PCoA1,
        Pair_ID,
        FUN = function(x) rev(x)
      ),
      yend = ave(
        PCoA2,
        Pair_ID,
        FUN = function(x) rev(x)
      )
    ),
    color = COL_PAIR,
    linewidth = 0.45,
    alpha = 0.55
  ) +

  stat_ellipse(
    aes(
      color = STATUS,
      group = STATUS
    ),
    type = "t",
    level = 0.95,
    linewidth = 0.8,
    show.legend = FALSE
  ) +

  geom_point(
    aes(
      color = STATUS,
      shape = STATUS
    ),
    size = 3,
    alpha = 0.90
  ) +

  scale_color_manual(
    values = c(
      "non-PMA" = COL_NONPMA,
      "PMA" = COL_PMA
    )
  ) +

  scale_shape_manual(
    values = c(
      "non-PMA" = 16,
      "PMA" = 17
    )
  ) +

  labs(
    x = paste0(
      "PCoA1 (",
      round(var1, 1),
      "%)"
    ),
    y = paste0(
      "PCoA2 (",
      round(var2, 1),
      "%)"
    ),
    color = NULL,
    shape = NULL
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(
    axis.text = element_text(
      color = "black"
    ),
    axis.title = element_text(
      color = "black",
      face = "bold"
    ),
    legend.position = c(0.82, 0.82),
    legend.direction = "vertical",
    legend.background = element_rect(
      fill = "white",
      color = NA
    ),
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size = 9)
  )

ggsave(
  file.path(
    FIGDIR,
    "Fig2B_Bray_PCoA.pdf"
  ),
  p,
  width = 5.5,
  height = 4.8
)

ggsave(
  file.path(
    FIGDIR,
    "Fig2B_Bray_PCoA.png"
  ),
  p,
  width = 5.5,
  height = 4.8,
  dpi = 600
)

cat("\n========================================\n")
cat("FILES SAVED\n")
cat("========================================\n")

cat(
  file.path(
    FIGDIR,
    "Fig2B_Bray_PCoA.pdf"
  ),
  "\n"
)

cat(
  file.path(
    OUTDIR,
    "bray_permanova_paired.txt"
  ),
  "\n"
)

cat(
  file.path(
    OUTDIR,
    "bray_permanova_depth_adjusted.txt"
  ),
  "\n"
)

cat(
  file.path(
    OUTDIR,
    "bray_permdisp_restricted.txt"
  ),
  "\n"
)

cat("\nDONE\n")
