# ============================================================
# Build Supplementary Tables S1-S6
# ISS MT-2 matched PMA vs non-PMA metagenomics study
#
# Outputs:
#   1. One formatted Excel workbook containing S1-S6
#   2. Repository-friendly TSV files for every component table
#
# Primary statistical analyses are NOT recomputed.
# Existing frozen analysis outputs are compiled.
# ============================================================

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages(
    "openxlsx",
    repos = "https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
})

options(stringsAsFactors = FALSE)

outdir <- "tables/supplementary"
tsvdir <- file.path(outdir, "tsv")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tsvdir, recursive = TRUE, showWarnings = FALSE)

cat("\n===============================================\n")
cat("BUILDING SUPPLEMENTARY TABLES S1-S6\n")
cat("===============================================\n\n")

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

read_tsv <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required file: ", path)
  }

  read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

write_tsv <- function(x, filename) {
  write.table(
    x,
    file.path(tsvdir, filename),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
}

# ------------------------------------------------------------
# S1 — COMPLETE SAMPLE METADATA
# ------------------------------------------------------------

cat("[S1] Complete MT-2 metadata\n")

meta <- read.csv(
  "analysis_ready/metadata/mt2_metadata_enriched.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_meta <- c(
  "Run",
  "SampleName",
  "Flight",
  "LocationCode",
  "STATUS",
  "Location"
)

if (!all(required_meta %in% names(meta))) {
  stop(
    "S1 metadata missing columns: ",
    paste(
      setdiff(required_meta, names(meta)),
      collapse = ", "
    )
  )
}

# Generate Pair_ID only if absent
if (!"Pair_ID" %in% names(meta)) {
  meta$Pair_ID <- sub("_P$", "", meta$SampleName)
}

S1 <- meta %>%
  select(
    Run,
    SampleName,
    STATUS,
    Pair_ID,
    Flight,
    LocationCode,
    Location,
    everything()
  ) %>%
  distinct(Run, .keep_all = TRUE) %>%
  arrange(Pair_ID, STATUS)

if (nrow(S1) != 64) {
  warning(
    "S1 contains ", nrow(S1),
    " rows; expected 64."
  )
}

write_tsv(
  S1,
  "Supplementary_Table_S1_complete_metadata.tsv"
)

# ------------------------------------------------------------
# S2 — SAMPLE-LEVEL QC
# ------------------------------------------------------------

cat("[S2] Sample-level multi-domain QC\n")

master <- read_tsv(
  "analysis_ready/metadata/master_sample_table.tsv"
)

alpha <- read_tsv(
  paste0(
    "analysis_ready/taxonomy/alpha_diversity/",
    "richness_depth_adjusted_input.tsv"
  )
)

bracken_qc <- read_tsv(
  "analysis_ready/taxonomy/bracken_sample_QC.tsv"
)

res <- read_tsv(
  paste0(
    "analysis_ready/deeparg/figure4/",
    "DeepARG_sample_resistome_metrics.tsv"
  )
)

# Keep only new alpha-diversity fields
alpha_keep <- alpha %>%
  select(
    Run,
    Observed,
    Shannon,
    Simpson,
    Bracken_depth
  ) %>%
  distinct(Run, .keep_all = TRUE)

# Keep nonduplicated DeepARG sample metrics
res_keep <- res %>%
  select(
    Run,
    DeepARG_total_16S,
    DeepARG_total_reads,
    ARG_subtype_richness,
    ARG_category_richness
  ) %>%
  distinct(Run, .keep_all = TRUE)

# Bracken QC: keep fields not already represented
bracken_names <- setdiff(
  names(bracken_qc),
  c(
    "SampleName", "Pair_ID", "STATUS",
    "Flight", "LocationCode", "Location"
  )
)

if ("Run" %in% bracken_names) {
  bracken_keep <- bracken_qc %>%
    select(all_of(bracken_names)) %>%
    distinct(Run, .keep_all = TRUE)
} else {
  bracken_keep <- NULL
}

S2 <- master %>%
  left_join(alpha_keep, by = "Run") %>%
  left_join(res_keep, by = "Run")

if (!is.null(bracken_keep)) {

  new_bcols <- setdiff(
    names(bracken_keep),
    names(S2)
  )

  if (length(new_bcols) > 0) {
    S2 <- S2 %>%
      left_join(
        bracken_keep %>%
          select(Run, all_of(new_bcols)),
        by = "Run"
      )
  }
}

# Remove analysis-only transformed columns from publication QC table
S2 <- S2 %>%
  select(
    -any_of(
      c(
        "log10_HostRemoved_pairs",
        "log10_Assembly_bp",
        "log10_Assembly_contigs",
        "log10_KEGG_assignments",
        "log10_AMRFinder_hits",
        "log10_Bracken_depth"
      )
    )
  ) %>%
  arrange(Pair_ID, STATUS)

if (nrow(S2) != 64) {
  warning(
    "S2 contains ", nrow(S2),
    " rows; expected 64."
  )
}

write_tsv(
  S2,
  "Supplementary_Table_S2_sample_QC.tsv"
)

# ------------------------------------------------------------
# S3 — COMPLETE ANCOM-BC2 RESULTS
# ------------------------------------------------------------

cat("[S3] Complete ANCOM-BC2 results\n")

S3_primary <- read_tsv(
  paste0(
    "analysis_ready/taxonomy/ancombc2_genus/",
    "ANCOMBC2_genus_20pct_primary_STATUS.tsv"
  )
)

# NOTE:
# Directory is sensitivity10 even though historical filename
# incorrectly contains "20pct".
S3_sensitivity <- read_tsv(
  paste0(
    "analysis_ready/taxonomy/",
    "ancombc2_genus_sensitivity10/",
    "ANCOMBC2_genus_20pct_primary_STATUS.tsv"
  )
)

S3_stable <- read_tsv(
  paste0(
    "analysis_ready/taxonomy/ancombc2_genus/",
    "ANCOMBC2_genus_stable_10pct_20pct.tsv"
  )
)

write_tsv(
  S3_primary,
  "Supplementary_Table_S3A_ANCOMBC2_primary_20pct.tsv"
)

write_tsv(
  S3_sensitivity,
  "Supplementary_Table_S3B_ANCOMBC2_sensitivity_10pct.tsv"
)

write_tsv(
  S3_stable,
  "Supplementary_Table_S3C_ANCOMBC2_stable_both.tsv"
)

# ------------------------------------------------------------
# S4 — DEEPARG PAIRED RESULTS
# ------------------------------------------------------------

cat("[S4] DeepARG paired ARG and resistance-category results\n")

S4_ARG <- read_tsv(
  paste0(
    "analysis_ready/deeparg/figure4/",
    "DeepARG_ARG_group_persistence.tsv"
  )
)

S4_category <- read_tsv(
  paste0(
    "analysis_ready/deeparg/figure4/",
    "DeepARG_category_paired_results.tsv"
  )
)

S4_sample <- read_tsv(
  paste0(
    "analysis_ready/deeparg/figure4/",
    "DeepARG_sample_resistome_metrics.tsv"
  )
)

write_tsv(
  S4_ARG,
  "Supplementary_Table_S4A_ARG_group_paired.tsv"
)

write_tsv(
  S4_category,
  "Supplementary_Table_S4B_resistance_category_paired.tsv"
)

write_tsv(
  S4_sample,
  "Supplementary_Table_S4C_sample_resistome_metrics.tsv"
)

# ------------------------------------------------------------
# S5 — CROSS-METHOD AMR CONCORDANCE
# ------------------------------------------------------------

cat("[S5] Cross-method AMR concordance\n")

S5_gene <- read_tsv(
  paste0(
    "analysis_ready/amr_concordance/",
    "CARD_candidate_concordance.tsv"
  )
)

S5_four <- read_tsv(
  paste0(
    "analysis_ready/amr_concordance/",
    "four_method_candidate_evidence.tsv"
  )
)

S5_class <- read_tsv(
  paste0(
    "analysis_ready/amr_concordance/",
    "class_prevalence_four_methods.tsv"
  )
)

S5_class_conc <- read_tsv(
  paste0(
    "analysis_ready/amr_concordance/",
    "class_concordance_PMA.tsv"
  )
)

S5_amrfinder <- read_tsv(
  paste0(
    "analysis_ready/amr_concordance/",
    "AMRFinder_gene_prevalence.tsv"
  )
)

write_tsv(
  S5_gene,
  "Supplementary_Table_S5A_gene_concordance.tsv"
)

write_tsv(
  S5_four,
  "Supplementary_Table_S5B_four_method_evidence.tsv"
)

write_tsv(
  S5_class,
  "Supplementary_Table_S5C_class_prevalence.tsv"
)

write_tsv(
  S5_class_conc,
  "Supplementary_Table_S5D_class_concordance_PMA.tsv"
)

write_tsv(
  S5_amrfinder,
  "Supplementary_Table_S5E_AMRFinder_gene_prevalence.tsv"
)

# ------------------------------------------------------------
# S6 — ARG-BEARING CONTIG TAXONOMY
# ------------------------------------------------------------

cat("[S6] ARG-bearing contig taxonomy and associations\n")

S6_contigs <- read_tsv(
  paste0(
    "analysis_ready/amrfinder/",
    "amrfinder_ARG_contig_taxonomy.tsv"
  )
)

S6_assoc <- read_tsv(
  paste0(
    "analysis_ready/arg_host/",
    "ARG_taxon_associations.tsv"
  )
)

S6_recurrent <- read_tsv(
  paste0(
    "analysis_ready/arg_host/",
    "ARG_taxon_recurrent_PMA.tsv"
  )
)

S6_selected <- read_tsv(
  paste0(
    "analysis_ready/arg_host/",
    "ARG_taxon_PMA_vs_nonPMA_selected.tsv"
  )
)

write_tsv(
  S6_contigs,
  "Supplementary_Table_S6A_ARG_contig_taxonomy.tsv"
)

write_tsv(
  S6_assoc,
  "Supplementary_Table_S6B_ARG_taxon_associations.tsv"
)

write_tsv(
  S6_recurrent,
  "Supplementary_Table_S6C_recurrent_PMA_ARG_taxon.tsv"
)

write_tsv(
  S6_selected,
  "Supplementary_Table_S6D_selected_PMA_vs_nonPMA.tsv"
)

# ============================================================
# EXCEL WORKBOOK
# ============================================================

cat("\nFormatting Excel workbook...\n")

wb <- createWorkbook(
  creator = "ISS MT-2 metagenomics study"
)

# ------------------------------------------------------------
# Styles
# ------------------------------------------------------------

titleStyle <- createStyle(
  fontSize = 14,
  fontColour = "#FFFFFF",
  fgFill = "#5B4B8A",
  textDecoration = "bold",
  halign = "left",
  valign = "center"
)

headerStyle <- createStyle(
  fontColour = "#FFFFFF",
  fgFill = "#404040",
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  wrapText = TRUE,
  border = "Bottom",
  borderColour = "#000000"
)

pmaStyle <- createStyle(
  fgFill = "#FFE2E2"
)

nonPmaStyle <- createStyle(
  fgFill = "#DDEEFF"
)

percentStyle <- createStyle(
  numFmt = "0.0%"
)

decimalStyle <- createStyle(
  numFmt = "0.000"
)

scientificStyle <- createStyle(
  numFmt = "0.00E+00"
)

wrapStyle <- createStyle(
  wrapText = TRUE,
  valign = "top"
)

# ------------------------------------------------------------
# README
# ------------------------------------------------------------

addWorksheet(wb, "README")

readme <- data.frame(
  Supplement = c(
    "S1",
    "S2",
    "S3",
    "S4",
    "S5",
    "S6"
  ),
  Description = c(
    paste(
      "Complete metadata for the 64 public MT-2 metagenomes,",
      "including Run, SampleName, PMA status, Pair_ID, flight,",
      "location code, and location."
    ),
    paste(
      "Sample-level sequencing-depth, assembly-yield, taxonomic,",
      "resistome, AMRFinderPlus, and KEGG QC metrics."
    ),
    paste(
      "Complete ANCOM-BC2 genus-level results for the primary",
      "20% prevalence analysis, 10% sensitivity analysis, and",
      "stable results shared by both analyses."
    ),
    paste(
      "DeepARG ARG-group and resistance-category paired results,",
      "prevalence, effect summaries, Wilcoxon statistics and BH-FDR."
    ),
    paste(
      "Cross-method gene- and resistance-class-level AMR",
      "concordance across DeepARG, ShortBRED, RGI and AMRFinderPlus."
    ),
    paste(
      "Complete ARG-bearing contig taxonomy and recurrent",
      "ARG-taxon association results."
    )
  ),
  stringsAsFactors = FALSE
)

writeData(
  wb,
  "README",
  "Supplementary Tables S1-S6",
  startRow = 1,
  startCol = 1
)

mergeCells(
  wb,
  "README",
  cols = 1:2,
  rows = 1
)

addStyle(
  wb,
  "README",
  titleStyle,
  rows = 1,
  cols = 1:2,
  gridExpand = TRUE
)

writeData(
  wb,
  "README",
  readme,
  startRow = 3,
  headerStyle = headerStyle
)

setColWidths(
  wb,
  "README",
  cols = 1,
  widths = 15
)

setColWidths(
  wb,
  "README",
  cols = 2,
  widths = 100
)

addStyle(
  wb,
  "README",
  wrapStyle,
  rows = 4:(nrow(readme) + 3),
  cols = 2,
  gridExpand = TRUE
)

freezePane(
  wb,
  "README",
  firstActiveRow = 4
)

# ------------------------------------------------------------
# Generic sheet writer
# ------------------------------------------------------------

add_table_sheet <- function(
  wb,
  sheet,
  title,
  dat,
  status_col = NULL
) {

  addWorksheet(wb, sheet)

  ncols <- ncol(dat)

  writeData(
    wb,
    sheet,
    title,
    startRow = 1,
    startCol = 1
  )

  mergeCells(
    wb,
    sheet,
    cols = 1:ncols,
    rows = 1
  )

  addStyle(
    wb,
    sheet,
    titleStyle,
    rows = 1,
    cols = 1:ncols,
    gridExpand = TRUE
  )

  writeData(
    wb,
    sheet,
    dat,
    startRow = 3,
    startCol = 1,
    headerStyle = headerStyle,
    withFilter = TRUE
  )

  freezePane(
    wb,
    sheet,
    firstActiveRow = 4
  )

  # Sensible widths
  setColWidths(
    wb,
    sheet,
    cols = 1:ncols,
    widths = "auto"
  )

  # Prevent enormous widths from long descriptions
  for (j in seq_len(ncols)) {

    vals <- as.character(dat[[j]])

    maxlen <- max(
      nchar(
        c(names(dat)[j], vals),
        type = "width"
      ),
      na.rm = TRUE
    )

    width <- min(
      max(maxlen + 2, 10),
      35
    )

    setColWidths(
      wb,
      sheet,
      cols = j,
      widths = width
    )
  }

  addStyle(
    wb,
    sheet,
    wrapStyle,
    rows = 4:(nrow(dat) + 3),
    cols = 1:ncols,
    gridExpand = TRUE,
    stack = TRUE
  )

  # Highlight PMA/non-PMA rows if STATUS exists
  if (!is.null(status_col) &&
      status_col %in% names(dat)) {

    sc <- which(names(dat) == status_col)

    pma_rows <- which(
      dat[[status_col]] == "PMA"
    ) + 3

    nonpma_rows <- which(
      dat[[status_col]] == "non-PMA"
    ) + 3

    if (length(pma_rows) > 0) {
      addStyle(
        wb,
        sheet,
        pmaStyle,
        rows = pma_rows,
        cols = 1:ncols,
        gridExpand = TRUE,
        stack = TRUE
      )
    }

    if (length(nonpma_rows) > 0) {
      addStyle(
        wb,
        sheet,
        nonPmaStyle,
        rows = nonpma_rows,
        cols = 1:ncols,
        gridExpand = TRUE,
        stack = TRUE
      )
    }
  }

  invisible(TRUE)
}

# ------------------------------------------------------------
# Write all workbook sheets
# ------------------------------------------------------------

add_table_sheet(
  wb,
  "S1_Metadata",
  "Supplementary Table S1. Complete metadata for 64 MT-2 metagenomes",
  S1,
  "STATUS"
)

add_table_sheet(
  wb,
  "S2_Sample_QC",
  paste(
    "Supplementary Table S2. Sample-level sequencing, assembly,",
    "taxonomy, resistome, AMRFinderPlus and KEGG QC metrics"
  ),
  S2,
  "STATUS"
)

add_table_sheet(
  wb,
  "S3A_ANCOM20",
  "Supplementary Table S3A. ANCOM-BC2 primary analysis (20% prevalence threshold)",
  S3_primary
)

add_table_sheet(
  wb,
  "S3B_ANCOM10",
  "Supplementary Table S3B. ANCOM-BC2 sensitivity analysis (10% prevalence threshold)",
  S3_sensitivity
)

add_table_sheet(
  wb,
  "S3C_Stable",
  "Supplementary Table S3C. Genera stable across 20% and 10% prevalence analyses",
  S3_stable
)

add_table_sheet(
  wb,
  "S4A_ARG_groups",
  "Supplementary Table S4A. DeepARG ARG-group paired results and prevalence",
  S4_ARG
)

add_table_sheet(
  wb,
  "S4B_ARG_classes",
  "Supplementary Table S4B. DeepARG resistance-category paired results and prevalence",
  S4_category
)

add_table_sheet(
  wb,
  "S4C_Sample_metrics",
  "Supplementary Table S4C. Sample-level DeepARG resistome metrics",
  S4_sample,
  "STATUS"
)

add_table_sheet(
  wb,
  "S5A_Gene",
  "Supplementary Table S5A. Cross-method candidate-gene concordance",
  S5_gene
)

add_table_sheet(
  wb,
  "S5B_4method",
  "Supplementary Table S5B. Four-method candidate ARG evidence",
  S5_four
)

add_table_sheet(
  wb,
  "S5C_Class_prev",
  "Supplementary Table S5C. Resistance-class prevalence across AMR methods",
  S5_class,
  "STATUS"
)

add_table_sheet(
  wb,
  "S5D_Class_PMA",
  "Supplementary Table S5D. PMA resistance-class concordance",
  S5_class_conc
)

add_table_sheet(
  wb,
  "S5E_AMRFinder",
  "Supplementary Table S5E. AMRFinderPlus gene prevalence",
  S5_amrfinder
)

add_table_sheet(
  wb,
  "S6A_ARG_contigs",
  "Supplementary Table S6A. Complete taxonomy of ARG-bearing contigs",
  S6_contigs,
  "STATUS"
)

add_table_sheet(
  wb,
  "S6B_Associations",
  "Supplementary Table S6B. Complete ARG-taxon association results",
  S6_assoc,
  "STATUS"
)

add_table_sheet(
  wb,
  "S6C_Recurrent_PMA",
  "Supplementary Table S6C. Recurrent PMA ARG-taxon associations",
  S6_recurrent,
  "STATUS"
)

add_table_sheet(
  wb,
  "S6D_Selected",
  "Supplementary Table S6D. Selected recurrent ARG-taxon associations in PMA and non-PMA samples",
  S6_selected,
  "STATUS"
)

# ------------------------------------------------------------
# Save workbook
# ------------------------------------------------------------

outfile <- file.path(
  outdir,
  "Supplementary_Tables_S1-S6.xlsx"
)

saveWorkbook(
  wb,
  outfile,
  overwrite = TRUE
)

# ------------------------------------------------------------
# Verification
# ------------------------------------------------------------

cat("\n===============================================\n")
cat("SUPPLEMENTARY TABLE PACKAGE COMPLETE\n")
cat("===============================================\n")

cat("\nExcel workbook:\n")
cat(outfile, "\n")

cat("\nRow counts:\n")
cat("S1 metadata:                 ", nrow(S1), "\n")
cat("S2 sample QC:                ", nrow(S2), "\n")
cat("S3 primary 20%:              ", nrow(S3_primary), "\n")
cat("S3 sensitivity 10%:          ", nrow(S3_sensitivity), "\n")
cat("S3 stable both:              ", nrow(S3_stable), "\n")
cat("S4 ARG groups:               ", nrow(S4_ARG), "\n")
cat("S4 resistance categories:    ", nrow(S4_category), "\n")
cat("S5 gene concordance:         ", nrow(S5_gene), "\n")
cat("S5 four-method evidence:     ", nrow(S5_four), "\n")
cat("S5 class prevalence:         ", nrow(S5_class), "\n")
cat("S6 ARG-bearing contigs:      ", nrow(S6_contigs), "\n")
cat("S6 ARG-taxon associations:   ", nrow(S6_assoc), "\n")
cat("S6 recurrent PMA:            ", nrow(S6_recurrent), "\n")

cat("\nWorkbook sheets:\n")
print(names(wb))

cat("\nTSV directory:\n")
cat(tsvdir, "\n")

cat("\nDone.\n")
