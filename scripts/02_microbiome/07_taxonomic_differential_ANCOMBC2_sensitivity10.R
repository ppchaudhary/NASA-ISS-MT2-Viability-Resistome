
.libPaths(c(path.expand("~/R/4.5"), .libPaths()))

suppressPackageStartupMessages({
    library(ANCOMBC)
    library(dplyr)
})

set.seed(20261005)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
count_file <- "analysis_ready/taxonomy/bracken_genus_counts.tsv"
meta_file  <- "analysis_ready/metadata/master_sample_table.tsv"

outdir <- "analysis_ready/taxonomy/ancombc2_genus_sensitivity10"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Read data
# ------------------------------------------------------------
counts <- read.delim(
    count_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

meta <- read.delim(
    meta_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

# Taxon identifier avoids duplicate-name problems
taxon_id <- paste0(counts$taxonomy_id, "|", counts$name)

X <- as.matrix(counts[, -(1:2)])
rownames(X) <- taxon_id

storage.mode(X) <- "numeric"

# ------------------------------------------------------------
# Metadata alignment
# ------------------------------------------------------------
meta <- meta[match(colnames(X), meta$Run), ]

if (any(is.na(meta$Run))) {
    stop("Some count-matrix samples are missing from metadata.")
}

if (!all(meta$Run == colnames(X))) {
    stop("Metadata/count matrix alignment failed.")
}

meta$STATUS <- factor(
    meta$STATUS,
    levels = c("non-PMA", "PMA")
)

meta$Pair_ID <- factor(meta$Pair_ID)

rownames(meta) <- meta$Run

cat("===== DATA CHECK =====\n")
cat("Samples:", ncol(X), "\n")
cat("Genera before filtering:", nrow(X), "\n")
cat("PMA:", sum(meta$STATUS == "PMA"), "\n")
cat("non-PMA:", sum(meta$STATUS == "non-PMA"), "\n")
cat("Pairs:", length(unique(meta$Pair_ID)), "\n")

# ------------------------------------------------------------
# Prevalence filter
# Primary threshold = 20%
# ------------------------------------------------------------
prev <- rowMeans(X > 0)

keep <- prev >= 0.10
X20 <- X[keep, , drop = FALSE]

cat("\nGenera retained at >=10% prevalence:",
    nrow(X20), "\n")

# Save prevalence table
prev_table <- data.frame(
    taxon = rownames(X),
    taxonomy_id = counts$taxonomy_id,
    genus = counts$name,
    prevalence = prev,
    samples_present = rowSums(X > 0),
    retained_20pct = keep
)

write.table(
    prev_table,
    file = file.path(outdir, "genus_prevalence.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ------------------------------------------------------------
# ANCOM-BC2
#
# Matched design:
# STATUS = fixed effect
# Pair_ID = random intercept
#
# non-PMA is reference.
# Positive STATUSPMA LFC = higher in PMA.
# Negative STATUSPMA LFC = higher in non-PMA.
# ------------------------------------------------------------

cat("\n===== RUNNING ANCOM-BC2 =====\n")

fit <- ancombc2(
    data = X20,
    taxa_are_rows = TRUE,
    meta_data = meta,
    fix_formula = "STATUS",
    rand_formula = "(1 | Pair_ID)",
    p_adj_method = "BH",
    pseudo = 0,
    pseudo_sens = TRUE,
    prv_cut = 0,
    lib_cut = 0,
    s0_perc = 0.05,
    group = "STATUS",
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    n_cl = 4,
    verbose = TRUE,
    global = FALSE,
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE
)

# Save complete R object for reproducibility
saveRDS(
    fit,
    file.path(outdir, "ANCOMBC2_genus_20pct_fit.rds")
)

# ------------------------------------------------------------
# Inspect/save primary result
# ------------------------------------------------------------
res <- fit$res

write.table(
    res,
    file = file.path(outdir, "ANCOMBC2_genus_20pct_full_results.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("\n===== RESULT COLUMNS =====\n")
print(colnames(res))

cat("\n===== RESULT DIMENSIONS =====\n")
print(dim(res))

cat("\n===== FIRST RESULTS =====\n")
print(head(res))

# ------------------------------------------------------------
# Identify STATUS coefficient automatically
# ------------------------------------------------------------
lfc_cols <- grep("^lfc_", colnames(res), value = TRUE)

cat("\nLFC columns:\n")
print(lfc_cols)

status_lfc <- grep("STATUSPMA", lfc_cols, value = TRUE)

if (length(status_lfc) != 1) {
    warning(
        "Could not uniquely identify STATUSPMA coefficient. ",
        "Inspect result columns before downstream plotting."
    )
} else {

    suffix <- sub("^lfc_", "", status_lfc)

    q_col <- paste0("q_", suffix)
    p_col <- paste0("p_", suffix)
    se_col <- paste0("se_", suffix)
    diff_col <- paste0("diff_", suffix)

    primary <- res

    # Attach clean taxonomy names
    split_tax <- strsplit(primary$taxon, "\\|", fixed = FALSE)

    primary$taxonomy_id <- sapply(split_tax, `[`, 1)
    primary$genus <- sapply(
        split_tax,
        function(z) paste(z[-1], collapse = "|")
    )

    primary$prevalence <- prev[primary$taxon]

    # 95% CI if SE available
    if (se_col %in% colnames(primary)) {
        primary$CI_low <-
            primary[[status_lfc]] - 1.96 * primary[[se_col]]

        primary$CI_high <-
            primary[[status_lfc]] + 1.96 * primary[[se_col]]
    }

    # Direction
    primary$Direction <- ifelse(
        primary[[status_lfc]] > 0,
        "Higher in PMA",
        "Higher in non-PMA"
    )

    # FDR significance
    if (q_col %in% colnames(primary)) {
        primary$FDR_significant <-
            !is.na(primary[[q_col]]) &
            primary[[q_col]] < 0.05
    }

    primary <- primary %>%
        arrange(.data[[q_col]])

    write.table(
        primary,
        file = file.path(
            outdir,
            "ANCOMBC2_genus_20pct_primary_STATUS.tsv"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    if ("FDR_significant" %in% colnames(primary)) {

        sig <- primary %>%
            filter(FDR_significant)

        write.table(
            sig,
            file = file.path(
                outdir,
                "ANCOMBC2_genus_20pct_FDR05.tsv"
            ),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )

        cat("\n===== FDR < 0.05 =====\n")
        cat("Significant genera:", nrow(sig), "\n")

        if (nrow(sig) > 0) {
            cat(
                "Higher in PMA:",
                sum(sig$Direction == "Higher in PMA"),
                "\n"
            )

            cat(
                "Higher in non-PMA:",
                sum(sig$Direction == "Higher in non-PMA"),
                "\n"
            )

            cat("\nTop significant results:\n")

            showcols <- intersect(
                c(
                    "genus",
                    status_lfc,
                    se_col,
                    p_col,
                    q_col,
                    "CI_low",
                    "CI_high",
                    "prevalence",
                    "Direction"
                ),
                colnames(sig)
            )

            print(
                head(sig[, showcols, drop = FALSE], 20),
                row.names = FALSE
            )
        }
    }
}

cat("\n===== ANALYSIS COMPLETE =====\n")
cat("Output directory:", outdir, "\n")

