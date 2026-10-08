
#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

INFILE = "analysis_ready/amrfinder/amrfinder_hits_long.tsv"
OUTDIR = Path("analysis_ready/amr_concordance")
OUTDIR.mkdir(parents=True, exist_ok=True)

x = pd.read_csv(INFILE, sep="\t")

# ------------------------------------------------------------
# Master sample universe from project metadata
# Important because AMRFinder has zero-hit samples.
# ------------------------------------------------------------

meta = pd.read_csv(
    "analysis_ready/metadata/master_sample_table.tsv",
    sep="\t"
)

samples = (
    meta[["Run", "STATUS", "Pair_ID"]]
    .drop_duplicates()
)

print("\n===== SAMPLE UNIVERSE =====")
print(samples.groupby("STATUS")["Run"].nunique())

# ------------------------------------------------------------
# Overall AMRFinder detection
# ------------------------------------------------------------

hit_runs = set(x["Run"].unique())

overall = (
    samples.assign(
        AMRFinder_detected=samples["Run"].isin(hit_runs)
    )
    .groupby("STATUS")
    .agg(
        N_samples=("Run", "nunique"),
        N_positive=("AMRFinder_detected", "sum")
    )
    .reset_index()
)

overall["Prevalence"] = (
    overall["N_positive"] /
    overall["N_samples"]
)

print("\n===== OVERALL AMRFINDER DETECTION =====")
print(overall.to_string(index=False))

# ------------------------------------------------------------
# Class prevalence
#
# A class is detected if >=1 AMRFinder hit belonging to that
# class occurs in a sample.
# ------------------------------------------------------------

class_detect = (
    x[["Run", "Class"]]
    .dropna()
    .drop_duplicates()
)

classes = sorted(class_detect["Class"].unique())

grid = (
    samples.assign(key=1)
    .merge(
        pd.DataFrame({"Class": classes, "key": 1}),
        on="key"
    )
    .drop(columns="key")
)

grid = grid.merge(
    class_detect.assign(Detected=1),
    on=["Run", "Class"],
    how="left"
)

grid["Detected"] = grid["Detected"].fillna(0).astype(int)

class_summary = (
    grid.groupby(["Class", "STATUS"])
    .agg(
        N_samples=("Run", "nunique"),
        N_detected=("Detected", "sum")
    )
    .reset_index()
)

class_summary["Prevalence"] = (
    class_summary["N_detected"] /
    class_summary["N_samples"]
)

class_wide = (
    class_summary
    .pivot(
        index="Class",
        columns="STATUS",
        values=["N_detected", "Prevalence"]
    )
)

class_wide.columns = [
    "_".join(map(str, c))
    for c in class_wide.columns
]

class_wide = (
    class_wide
    .reset_index()
    .sort_values(
        ["Prevalence_PMA", "Prevalence_non-PMA"],
        ascending=False
    )
)

# ------------------------------------------------------------
# Gene-symbol prevalence
# ------------------------------------------------------------

gene_detect = (
    x[["Run", "STATUS", "Element symbol", "Class"]]
    .dropna(subset=["Element symbol"])
    .drop_duplicates()
)

gene_summary = (
    gene_detect
    .groupby(["Element symbol", "Class", "STATUS"])
    .agg(
        N_detected=("Run", "nunique")
    )
    .reset_index()
)

gene_summary["Prevalence"] = (
    gene_summary["N_detected"] / 32.0
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

overall.to_csv(
    OUTDIR / "AMRFinder_overall_detection.tsv",
    sep="\t",
    index=False
)

class_summary.to_csv(
    OUTDIR / "AMRFinder_class_prevalence_long.tsv",
    sep="\t",
    index=False
)

class_wide.to_csv(
    OUTDIR / "AMRFinder_class_prevalence_wide.tsv",
    sep="\t",
    index=False
)

gene_summary.to_csv(
    OUTDIR / "AMRFinder_gene_prevalence.tsv",
    sep="\t",
    index=False
)

# ------------------------------------------------------------
# Console summaries
# ------------------------------------------------------------

print("\n===== AMRFINDER CLASS PREVALENCE =====")
print(
    class_wide.head(25).to_string(index=False)
)

print("\n===== TOP PMA AMRFINDER GENES =====")

pma_genes = (
    gene_summary[
        gene_summary["STATUS"] == "PMA"
    ]
    .sort_values(
        ["N_detected", "Element symbol"],
        ascending=[False, True]
    )
)

print(
    pma_genes.head(30).to_string(index=False)
)

print("\nSaved AMRFinder concordance summaries.")

