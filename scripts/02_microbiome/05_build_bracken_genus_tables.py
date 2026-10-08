#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

ROOT = Path(".")
BRACKEN = ROOT / "bracken"
OUT = ROOT / "analysis_ready" / "taxonomy"
META = ROOT / "analysis_ready" / "metadata" / "master_sample_table.tsv"

OUT.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(META, sep="\t")

required = {"Run", "STATUS", "Pair_ID"}
missing = required - set(meta.columns)
if missing:
    raise ValueError(f"Missing metadata columns: {missing}")

meta = meta.drop_duplicates("Run").copy()

records = []

files = sorted(BRACKEN.rglob("*.S.kreport"))
print(f"Found {len(files)} Bracken kreport files")

if len(files) != 64:
    raise ValueError(f"Expected 64 kreport files, found {len(files)}")

for f in files:

    run = f.name.replace(".bracken.S.kreport", "")

    dat = pd.read_csv(
        f,
        sep="\t",
        header=None,
        names=[
            "percent",
            "clade_reads",
            "taxon_reads",
            "rank",
            "taxonomy_id",
            "name"
        ]
    )

    # Remove indentation used to display taxonomy hierarchy
    dat["name"] = dat["name"].astype(str).str.strip()

    # Exact genus rank only.
    # Do NOT include G1, G2, etc.
    genus = dat.loc[dat["rank"] == "G"].copy()

    genus["Run"] = run

    records.append(
        genus[
            [
                "Run",
                "taxonomy_id",
                "name",
                "clade_reads",
                "percent"
            ]
        ]
    )

long = pd.concat(records, ignore_index=True)

# Attach sample metadata
long = long.merge(
    meta,
    on="Run",
    how="left",
    validate="many_to_one"
)

if long["STATUS"].isna().any():
    bad = long.loc[long["STATUS"].isna(), "Run"].unique()
    raise ValueError(f"Metadata missing for runs: {bad}")

# Save long table
long_file = OUT / "bracken_genus_long.tsv"
long.to_csv(long_file, sep="\t", index=False)

# Genus x sample count matrix
counts = (
    long.pivot_table(
        index=["taxonomy_id", "name"],
        columns="Run",
        values="clade_reads",
        aggfunc="sum",
        fill_value=0
    )
    .reset_index()
)

counts_file = OUT / "bracken_genus_counts.tsv"
counts.to_csv(counts_file, sep="\t", index=False)

# Relative abundance matrix calculated within genus-level abundance
mat = counts.set_index(["taxonomy_id", "name"]).astype(float)

sample_totals = mat.sum(axis=0)

rel = mat.div(sample_totals, axis=1)

rel = rel.reset_index()

rel_file = OUT / "bracken_genus_relative_abundance.tsv"
rel.to_csv(rel_file, sep="\t", index=False)

# QC summary
qc = (
    long.groupby("Run")
        .agg(
            Genus_total_reads=("clade_reads", "sum"),
            Genus_richness=("taxonomy_id", "nunique")
        )
        .reset_index()
        .merge(
            meta[["Run", "STATUS", "Pair_ID"]],
            on="Run",
            how="left"
        )
)

qc_file = OUT / "bracken_genus_QC.tsv"
qc.to_csv(qc_file, sep="\t", index=False)

print("\n===== GENUS TABLE BUILD COMPLETE =====")
print("Long table:", long.shape)
print("Count matrix:", counts.shape)
print("Unique genera:", long["taxonomy_id"].nunique())
print("Runs:", long["Run"].nunique())

print("\nSTATUS:")
print(qc["STATUS"].value_counts())

print("\nComplete Pair_IDs:")
paircheck = (
    qc.groupby("Pair_ID")["STATUS"]
      .nunique()
)
print((paircheck == 2).sum(), "of", paircheck.size)

print("\nGenus richness by STATUS:")
print(
    qc.groupby("STATUS")["Genus_richness"]
      .describe()[["min", "25%", "50%", "mean", "75%", "max"]]
)

print("\nFiles written:")
print(long_file)
print(counts_file)
print(rel_file)
print(qc_file)
