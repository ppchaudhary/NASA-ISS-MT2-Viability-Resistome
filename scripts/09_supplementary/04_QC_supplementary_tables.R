suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
})

xlsx <- "tables/supplementary/Supplementary_Tables_S1-S6.xlsx"
out  <- "tables/supplementary/Supplementary_Tables_QC_report.txt"

sink(out, split = TRUE)

cat("====================================================\n")
cat("SUPPLEMENTARY TABLES S1-S6 — FINAL QC\n")
cat("====================================================\n\n")

PASS <- function(x) ifelse(isTRUE(x), "PASS", "FAIL")

# ------------------------------------------------------------
# Read sheets
# Data begin at Excel row 3; skip title + blank row
# ------------------------------------------------------------

read_sheet <- function(sheet) {
  read.xlsx(
    xlsx,
    sheet = sheet,
    startRow = 3,
    check.names = FALSE
  )
}

S1  <- read_sheet("S1_Metadata")
S2  <- read_sheet("S2_Sample_QC")
S3A <- read_sheet("S3A_ANCOM20")
S3B <- read_sheet("S3B_ANCOM10")
S3C <- read_sheet("S3C_Stable")
S4A <- read_sheet("S4A_ARG_groups")
S4B <- read_sheet("S4B_ARG_classes")
S5A <- read_sheet("S5A_Gene")
S5B <- read_sheet("S5B_4method")
S5C <- read_sheet("S5C_Class_prev")
S6A <- read_sheet("S6A_ARG_contigs")
S6B <- read_sheet("S6B_Associations")
S6C <- read_sheet("S6C_Recurrent_PMA")

# ------------------------------------------------------------
# S1
# ------------------------------------------------------------

cat("S1 — METADATA\n")
cat("Rows = ", nrow(S1), " [", PASS(nrow(S1) == 64), "]\n", sep="")
cat("Unique Runs = ", n_distinct(S1$Run),
    " [", PASS(n_distinct(S1$Run) == 64), "]\n", sep="")
cat("Unique Pair_IDs = ", n_distinct(S1$Pair_ID),
    " [", PASS(n_distinct(S1$Pair_ID) == 32), "]\n", sep="")

status1 <- table(S1$STATUS)
print(status1)

cat(
  "32 PMA + 32 non-PMA = ",
  PASS(
    all(c("PMA", "non-PMA") %in% names(status1)) &&
    status1["PMA"] == 32 &&
    status1["non-PMA"] == 32
  ),
  "\n"
)

paircheck <- S1 %>%
  count(Pair_ID, STATUS) %>%
  tidyr::pivot_wider(
    names_from = STATUS,
    values_from = n,
    values_fill = 0
  )

cat(
  "All 32 pairs complete = ",
  PASS(
    nrow(paircheck) == 32 &&
    all(paircheck$PMA == 1) &&
    all(paircheck$`non-PMA` == 1)
  ),
  "\n\n"
)

# ------------------------------------------------------------
# S2
# ------------------------------------------------------------

cat("S2 — SAMPLE QC\n")
cat("Rows = ", nrow(S2), " [", PASS(nrow(S2) == 64), "]\n", sep="")
cat("Unique Runs = ", n_distinct(S2$Run),
    " [", PASS(n_distinct(S2$Run) == 64), "]\n", sep="")

required_S2 <- c(
  "HostRemoved_read_pairs",
  "Assembly_total_bp",
  "AMRFinder_hits",
  "KEGG_valid_assignments",
  "Observed",
  "Shannon",
  "Simpson",
  "DeepARG_total_16S",
  "ARG_subtype_richness",
  "ARG_category_richness"
)

cat(
  "Required QC columns present = ",
  PASS(all(required_S2 %in% names(S2))),
  "\n"
)

for (v in intersect(required_S2, names(S2))) {
  cat(
    sprintf(
      "%-28s missing=%d\n",
      v,
      sum(is.na(S2[[v]]))
    )
  )
}

cat("\n")

# ------------------------------------------------------------
# S3
# ------------------------------------------------------------

cat("S3 — ANCOM-BC2\n")
cat("Primary rows = ", nrow(S3A),
    " [", PASS(nrow(S3A) == 235), "]\n", sep="")
cat("Sensitivity rows = ", nrow(S3B),
    " [", PASS(nrow(S3B) == 254), "]\n", sep="")
cat("Stable rows = ", nrow(S3C),
    " [", PASS(nrow(S3C) == 48), "]\n", sep="")

if ("q_STATUSPMA" %in% names(S3A)) {
  cat(
    "Primary q-values in [0,1] = ",
    PASS(
      all(
        S3A$q_STATUSPMA >= 0 &
        S3A$q_STATUSPMA <= 1,
        na.rm = TRUE
      )
    ),
    "\n"
  )
}

if ("q_STATUSPMA" %in% names(S3B)) {
  cat(
    "Sensitivity q-values in [0,1] = ",
    PASS(
      all(
        S3B$q_STATUSPMA >= 0 &
        S3B$q_STATUSPMA <= 1,
        na.rm = TRUE
      )
    ),
    "\n"
  )
}

if ("same_direction" %in% names(S3C)) {
  cat(
    "All stable genera same direction = ",
    PASS(all(S3C$same_direction %in% c(TRUE, "TRUE"))),
    "\n"
  )
}

cat("\n")

# ------------------------------------------------------------
# S4
# ------------------------------------------------------------

cat("S4 — DEEPARG\n")
cat("ARG groups = ", nrow(S4A),
    " [", PASS(nrow(S4A) == 412), "]\n", sep="")
cat("Resistance categories = ", nrow(S4B),
    " [", PASS(nrow(S4B) == 39), "]\n", sep="")

for (dat in list(S4A, S4B)) {

  if ("FDR_BH" %in% names(dat)) {
    cat(
      "FDR range valid = ",
      PASS(
        all(
          dat$FDR_BH >= 0 &
          dat$FDR_BH <= 1,
          na.rm = TRUE
        )
      ),
      "\n"
    )
  }

  for (v in c("Prevalence_PMA", "Prevalence_non-PMA")) {
    if (v %in% names(dat)) {
      cat(
        v, " range valid = ",
        PASS(
          all(
            dat[[v]] >= 0 &
            dat[[v]] <= 1,
            na.rm = TRUE
          )
        ),
        "\n",
        sep=""
      )
    }
  }
}

cat("\n")

# ------------------------------------------------------------
# S5
# ------------------------------------------------------------

cat("S5 — AMR CONCORDANCE\n")
cat("Candidate genes = ", nrow(S5A),
    " [", PASS(nrow(S5A) == 7), "]\n", sep="")
cat("Four-method candidates = ", nrow(S5B),
    " [", PASS(nrow(S5B) == 7), "]\n", sep="")
cat("Class prevalence rows = ", nrow(S5C),
    " [", PASS(nrow(S5C) == 80), "]\n", sep="")

if ("Prevalence" %in% names(S5C)) {
  cat(
    "Class prevalence in [0,1] = ",
    PASS(
      all(
        S5C$Prevalence >= 0 &
        S5C$Prevalence <= 1,
        na.rm = TRUE
      )
    ),
    "\n"
  )
}

cat("\n")

# ------------------------------------------------------------
# S6
# ------------------------------------------------------------

cat("S6 — ARG-BEARING CONTIG TAXONOMY\n")
cat("ARG records = ", nrow(S6A),
    " [", PASS(nrow(S6A) == 811), "]\n", sep="")
cat("ARG-taxon associations = ", nrow(S6B),
    " [", PASS(nrow(S6B) == 224), "]\n", sep="")
cat("Recurrent PMA associations = ", nrow(S6C),
    " [", PASS(nrow(S6C) == 10), "]\n", sep="")

cat("Unique ARG contigs = ",
    n_distinct(S6A$`Contig id`), "\n")

cat(
  "Missing Element symbol = ",
  sum(is.na(S6A$`Element symbol`) |
      S6A$`Element symbol` == ""),
  "\n"
)

cat(
  "Missing Taxon_name = ",
  sum(is.na(S6A$Taxon_name) |
      S6A$Taxon_name == ""),
  "\n"
)

if ("Kraken_classified" %in% names(S6A)) {
  cat("\nKraken classification:\n")
  print(table(S6A$Kraken_classified, useNA = "ifany"))
}

if ("Prevalence" %in% names(S6B)) {
  cat(
    "Association prevalence in [0,1] = ",
    PASS(
      all(
        S6B$Prevalence >= 0 &
        S6B$Prevalence <= 1,
        na.rm = TRUE
      )
    ),
    "\n"
  )
}

cat("\n====================================================\n")
cat("QC COMPLETE\n")
cat("Review any line marked FAIL before submission.\n")
cat("====================================================\n")

sink()
