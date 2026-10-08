#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

ROOT = Path.cwd()

COUNT = ROOT / "analysis_ready/taxonomy/bracken_species_counts.tsv"
MASTER = ROOT / "analysis_ready/metadata/master_sample_table.tsv"
OUT = ROOT / "analysis_ready/taxonomy/bracken_sample_QC.tsv"

# --------------------------------------------------
# Load Bracken matrix
# --------------------------------------------------
raw = pd.read_csv(COUNT, sep="\t")

print("========================================")
print("BRACKEN RAW MATRIX")
print("========================================")
print("Raw shape:", raw.shape)
print("First column:", raw.columns[0])
print("Number of sample columns:", len(raw.columns) - 1)

# Bracken matrix contains two taxonomy columns:
# taxonomy_id and name. All remaining columns are samples.
id_cols = ["taxonomy_id", "name"]

missing_id_cols = [c for c in id_cols if c not in raw.columns]
if missing_id_cols:
    raise ValueError(f"Missing taxonomy columns: {missing_id_cols}")

sample_cols = [c for c in raw.columns if c not in id_cols]

raw[sample_cols] = raw[sample_cols].apply(
    pd.to_numeric, errors="coerce"
).fillna(0)

# Use taxonomy_id as unique feature identifier.
# rows = samples; columns = species
x = raw.set_index("taxonomy_id")[sample_cols].T

x.index.name = "Run"

print("\n========================================")
print("BRACKEN TRANSPOSED MATRIX")
print("========================================")
print("Shape:", x.shape)
print("Samples:", x.shape[0])
print("Species:", x.shape[1])

# --------------------------------------------------
# Load master metadata
# --------------------------------------------------
m = pd.read_csv(MASTER, sep="\t")

# --------------------------------------------------
# Sample-level QC
# --------------------------------------------------
qc = pd.DataFrame({
    "Run": x.index.astype(str),
    "Bracken_total_estimated_reads": x.sum(axis=1).values,
    "Bracken_species_richness": (x > 0).sum(axis=1).values
})

qc = qc.merge(
    m[[
        "Run",
        "SampleName",
        "Pair_ID",
        "STATUS",
        "HostRemoved_read_pairs",
        "Assembly_total_bp"
    ]],
    on="Run",
    how="left",
    validate="one_to_one"
)

qc.to_csv(OUT, sep="\t", index=False)

# --------------------------------------------------
# QC checks
# --------------------------------------------------
print("\n========================================")
print("BRACKEN SAMPLE QC")
print("========================================")

print("Rows:", len(qc))
print("Unique runs:", qc["Run"].nunique())

print("\nSTATUS:")
print(qc["STATUS"].value_counts(dropna=False))

print("\nMissing metadata:")
print(qc["STATUS"].isna().sum())

expected = set(m["Run"])
observed = set(qc["Run"])

print("\nMetadata samples missing Bracken:")
print(sorted(expected - observed))

print("\nUnexpected Bracken samples:")
print(sorted(observed - expected))

# Pair check
paircheck = (
    qc.dropna(subset=["Pair_ID"])
      .groupby("Pair_ID")["STATUS"]
      .nunique()
)

print("\nUnique Pair_IDs:", qc["Pair_ID"].nunique())
print("Complete PMA/non-PMA pairs:", (paircheck == 2).sum())

# --------------------------------------------------
# Group summaries
# --------------------------------------------------
print("\n========================================")
print("GROUP MEDIANS")
print("========================================")

print(
    qc.groupby("STATUS")[
        [
            "HostRemoved_read_pairs",
            "Bracken_total_estimated_reads",
            "Bracken_species_richness"
        ]
    ]
    .median()
    .round(1)
    .T
)

# --------------------------------------------------
# Depth correlations
# --------------------------------------------------
print("\n========================================")
print("SPEARMAN CORRELATIONS")
print("========================================")

cols = [
    "HostRemoved_read_pairs",
    "Bracken_total_estimated_reads",
    "Bracken_species_richness"
]

print(
    qc[cols]
    .corr(method="spearman")
    .round(3)
    .to_string()
)

# --------------------------------------------------
# Lowest-depth samples
# --------------------------------------------------
print("\n========================================")
print("LOWEST BRACKEN DEPTH")
print("========================================")

print(
    qc.sort_values("Bracken_total_estimated_reads")[
        [
            "Run",
            "SampleName",
            "Pair_ID",
            "STATUS",
            "HostRemoved_read_pairs",
            "Bracken_total_estimated_reads",
            "Bracken_species_richness"
        ]
    ]
    .head(15)
    .to_string(index=False)
)

print("\nSaved:", OUT)
print("\nDONE")
