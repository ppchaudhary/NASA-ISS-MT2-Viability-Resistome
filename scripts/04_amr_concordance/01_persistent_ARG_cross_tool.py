#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

OUT = Path("analysis_ready/amr_concordance")
OUT.mkdir(parents=True, exist_ok=True)

# Exact/clean mappings identified from CARD annotations.
mapping = pd.DataFrame([
    ["RPOB2", "ARO_3000501", "rpoB2"],
    ["BACA",  "ARO_3002986", "bacA"],
    ["NORA",  "ARO_3000391", "norA"],
    ["MSRA",  "ARO_3000251", "msrA"],
    ["UGD",   "ARO_3003577", "ugd"],
    ["DFRC",  "ARO_3002865", "dfrC"],
    ["ACRB",  "ARO_3000216", "acrB"],
], columns=["DeepARG_ARG", "ARO", "CARD_ARG"])

mapping.to_csv(
    OUT / "persistent_ARG_exact_mapping.tsv",
    sep="\t", index=False
)

s = pd.read_csv(
    "analysis_ready/shortbred/shortbred_long.tsv",
    sep="\t"
)

r = pd.read_csv(
    "analysis_ready/rgi/rgi_gene_long.tsv",
    sep="\t"
)

# ----------------------------------------------------------
# ShortBRED
# Detection: Hits > 0
# Quantification: Count
# ----------------------------------------------------------

sb = s.merge(mapping, on="ARO", how="inner")

sb_summary = (
    sb.groupby(["ARO", "CARD_ARG", "STATUS"])
      .agg(
          N_samples=("Run", "nunique"),
          N_detected=("Hits", lambda x: (x > 0).sum()),
          Median_Count=("Count", "median"),
          Mean_Count=("Count", "mean"),
          Median_Hits=("Hits", "median")
      )
      .reset_index()
)

# Ensure prevalence denominator = 32 samples/status.
sb_summary["Prevalence"] = sb_summary["N_detected"] / 32.0
sb_summary["Tool"] = "ShortBRED"

# ----------------------------------------------------------
# RGI
# Each row represents an ARO/sample mapping observation.
# Complete missing sample×ARO combinations as zero.
# ----------------------------------------------------------

sample_meta = (
    r[["Run", "STATUS", "Pair_ID"]]
    .drop_duplicates()
)

aros = mapping[["ARO", "CARD_ARG"]].drop_duplicates()

grid = (
    sample_meta.assign(key=1)
    .merge(aros.assign(key=1), on="key")
    .drop(columns="key")
)

rsub = (
    r[r["ARO"].isin(mapping["ARO"])]
    .groupby(["Run", "ARO"], as_index=False)
    .agg(
        All_Mapped_Reads=("All Mapped Reads", "sum"),
        Percent_Coverage=("Average Percent Coverage", "max")
    )
)

rg = grid.merge(
    rsub,
    on=["Run", "ARO"],
    how="left"
)

rg[["All_Mapped_Reads", "Percent_Coverage"]] = (
    rg[["All_Mapped_Reads", "Percent_Coverage"]].fillna(0)
)

rgi_summary = (
    rg.groupby(["ARO", "CARD_ARG", "STATUS"])
      .agg(
          N_samples=("Run", "nunique"),
          N_detected=("All_Mapped_Reads", lambda x: (x > 0).sum()),
          Median_Mapped_Reads=("All_Mapped_Reads", "median"),
          Mean_Mapped_Reads=("All_Mapped_Reads", "mean"),
          Median_Percent_Coverage=("Percent_Coverage", "median")
      )
      .reset_index()
)

rgi_summary["Prevalence"] = rgi_summary["N_detected"] / 32.0
rgi_summary["Tool"] = "RGI"

sb_summary.to_csv(
    OUT / "ShortBRED_persistent_ARG_summary.tsv",
    sep="\t", index=False
)

rgi_summary.to_csv(
    OUT / "RGI_persistent_ARG_summary.tsv",
    sep="\t", index=False
)

print("\n===== EXACT ARG MAPPING =====")
print(mapping.to_string(index=False))

print("\n===== SHORTBRED PREVALENCE =====")
print(
    sb_summary[
        ["CARD_ARG", "STATUS", "N_detected",
         "Prevalence", "Median_Count"]
    ].to_string(index=False)
)

print("\n===== RGI PREVALENCE =====")
print(
    rgi_summary[
        ["CARD_ARG", "STATUS", "N_detected",
         "Prevalence", "Median_Mapped_Reads",
         "Median_Percent_Coverage"]
    ].to_string(index=False)
)

