#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path

INFILE = "analysis_ready/amrfinder/amrfinder_ARG_contig_taxonomy.tsv"
OUT = Path("analysis_ready/arg_host")
OUT.mkdir(parents=True, exist_ok=True)

# ============================================================
# Load
# ============================================================

x = pd.read_csv(INFILE, sep="\t")

print("\n===== INPUT =====")
print("Rows:", len(x))
print("Columns:")
print("\n".join(x.columns))

# ============================================================
# Identify important columns
# ============================================================

def pick(candidates):
    for c in candidates:
        if c in x.columns:
            return c
    return None

arg_col = pick([
    "Element symbol",
    "ARG",
    "Gene",
    "gene"
])

run_col = pick(["Run"])
status_col = pick(["STATUS"])
contig_col = pick([
    "Contig id",
    "Contig",
    "contig"
])

taxid_col = pick([
    "TaxID",
    "taxid",
    "taxonomy_id"
])

taxname_col = pick([
    "TaxName",
    "Taxon",
    "taxon",
    "Name",
    "name",
    "ScientificName",
    "scientific_name",
    "Taxon_name"
])

rank_col = pick([
    "Rank",
    "rank",
    "TaxRank",
    "tax_rank",
    "Taxon_rank"
])

print("\n===== DETECTED COLUMNS =====")
print("ARG:", arg_col)
print("Run:", run_col)
print("STATUS:", status_col)
print("Contig:", contig_col)
print("TaxID:", taxid_col)
print("Taxon:", taxname_col)
print("Rank:", rank_col)

required = [arg_col, run_col, status_col, contig_col, taxname_col]

if any(v is None for v in required):
    raise SystemExit(
        "\nRequired column not detected. "
        "See column listing above."
    )

# ============================================================
# Basic classified dataset
# ============================================================

z = x.copy()

# Unique contig identifier across independently assembled runs
z["Run_Contig_ID"] = (
    z[run_col].astype(str) + "::" +
    z[contig_col].astype(str)
)

z[taxname_col] = (
    z[taxname_col]
    .astype(str)
    .str.strip()
)

bad_taxa = {
    "",
    "nan",
    "NA",
    "N/A",
    "unclassified",
    "Unclassified",
    "unclassified sequences",
    "root"
}

z["Classified"] = ~z[taxname_col].isin(bad_taxa)

classified = z[z["Classified"]].copy()

print("\n===== CLASSIFICATION QC =====")

qc = (
    z.groupby(status_col)
    .agg(
        ARG_hits=(arg_col, "size"),
        ARG_contigs=("Run_Contig_ID", "nunique"),
        ARG_symbols=(arg_col, "nunique")
    )
)

qc["Classified_hits"] = (
    classified.groupby(status_col)
    .size()
)

qc["Classified_fraction"] = (
    qc["Classified_hits"] /
    qc["ARG_hits"]
)

print(qc.to_string())

qc.to_csv(
    OUT / "ARG_contig_taxonomy_QC.tsv",
    sep="\t"
)

# ============================================================
# Top taxonomic assignments by STATUS
# Count unique ARG-bearing contigs, not duplicated hit rows
# ============================================================

contig_tax = (
    classified[
        [run_col, status_col, contig_col, "Run_Contig_ID", taxname_col]
    ]
    .drop_duplicates()
)

tax_summary = (
    contig_tax
    .groupby([status_col, taxname_col])
    .agg(
        ARG_contigs=("Run_Contig_ID", "nunique"),
        Samples=(run_col, "nunique")
    )
    .reset_index()
)

tax_summary["Total_ARG_contigs"] = (
    tax_summary.groupby(status_col)["ARG_contigs"]
    .transform("sum")
)

tax_summary["Proportion"] = (
    tax_summary["ARG_contigs"] /
    tax_summary["Total_ARG_contigs"]
)

tax_summary = tax_summary.sort_values(
    [status_col, "ARG_contigs"],
    ascending=[True, False]
)

tax_summary.to_csv(
    OUT / "ARG_contig_taxon_composition.tsv",
    sep="\t",
    index=False
)

print("\n===== TOP TAXA: PMA =====")
print(
    tax_summary[
        tax_summary[status_col] == "PMA"
    ].head(20).to_string(index=False)
)

print("\n===== TOP TAXA: non-PMA =====")
print(
    tax_summary[
        tax_summary[status_col] == "non-PMA"
    ].head(20).to_string(index=False)
)

# ============================================================
# ARG-taxon associations
#
# Count recurrence by independent samples, not raw hit count.
# ============================================================

pairs = (
    classified[
        [
            run_col,
            status_col,
            contig_col,
            "Run_Contig_ID",
            arg_col,
            taxname_col
        ]
    ]
    .drop_duplicates()
)

pair_summary = (
    pairs
    .groupby(
        [status_col, arg_col, taxname_col]
    )
    .agg(
        Samples=(run_col, "nunique"),
        Contigs=("Run_Contig_ID", "nunique")
    )
    .reset_index()
)

pair_summary["Prevalence"] = (
    pair_summary["Samples"] / 32
)

pair_summary = pair_summary.sort_values(
    [status_col, "Samples", "Contigs"],
    ascending=[True, False, False]
)

pair_summary.to_csv(
    OUT / "ARG_taxon_associations.tsv",
    sep="\t",
    index=False
)

print("\n===== RECURRENT PMA ARG-TAXON PAIRS =====")

pma_pairs = (
    pair_summary[
        pair_summary[status_col] == "PMA"
    ]
    .sort_values(
        ["Samples", "Contigs"],
        ascending=False
    )
)

print(
    pma_pairs.head(40).to_string(index=False)
)

# ============================================================
# Recurrent PMA subset
# Require >=2 independent PMA samples
# ============================================================

recurrent = pma_pairs[
    pma_pairs["Samples"] >= 2
].copy()

recurrent.to_csv(
    OUT / "ARG_taxon_recurrent_PMA.tsv",
    sep="\t",
    index=False
)

print("\n===== PMA RECURRENCE SUMMARY =====")
print("PMA ARG-taxon pairs total:", len(pma_pairs))
print("Pairs in >=2 PMA samples:", len(recurrent))
print("Unique recurrent ARGs:", recurrent[arg_col].nunique())
print("Unique recurrent taxa:", recurrent[taxname_col].nunique())

print("\n===== TOP PMA ARG SYMBOLS =====")

arg_pma = (
    pairs[
        pairs[status_col] == "PMA"
    ]
    .groupby(arg_col)
    .agg(
        Samples=(run_col, "nunique"),
        Contigs=("Run_Contig_ID", "nunique"),
        Taxa=(taxname_col, "nunique")
    )
    .reset_index()
    .sort_values(
        ["Samples", "Contigs"],
        ascending=False
    )
)

print(arg_pma.head(25).to_string(index=False))

arg_pma.to_csv(
    OUT / "PMA_ARG_summary.tsv",
    sep="\t",
    index=False
)

print("\nSaved outputs to:", OUT)
