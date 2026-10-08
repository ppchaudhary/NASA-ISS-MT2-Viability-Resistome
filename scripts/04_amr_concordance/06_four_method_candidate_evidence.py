#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

OUT = Path("analysis_ready/amr_concordance")
OUT.mkdir(parents=True, exist_ok=True)

# ============================================================
# Read harmonized three-method table
# ============================================================

card = pd.read_csv(
    OUT / "CARD_candidate_concordance.tsv",
    sep="\t"
)

amr = pd.read_csv(
    "analysis_ready/amrfinder/amrfinder_hits_long.tsv",
    sep="\t"
)

# ============================================================
# AMRFinder candidate aliases
#
# Only biologically equivalent/recognizable gene-symbol aliases
# are used. We do NOT force unmatched CARD genes to AMRFinder
# calls.
# ============================================================

aliases = {
    "msrA": ["msr(A)", "msrA"],
    "dfrC": ["dfrC"],
    "norA": ["norA"],
    "rpoB2": ["rpoB2"],
    "acrB": ["acrB"],
    "bacA": ["bacA"],
    "ugd": ["ugd"]
}

rows = []

for _, r in card.iterrows():

    gene = r["ARG"]

    possible = {
        str(x).strip().lower()
        for x in aliases.get(gene, [])
    }

    z = amr[
        amr["Element symbol"]
        .astype(str)
        .str.strip()
        .str.lower()
        .isin(possible)
    ].copy()

    pma_runs = z.loc[
        z["STATUS"] == "PMA",
        "Run"
    ].nunique()

    non_runs = z.loc[
        z["STATUS"] == "non-PMA",
        "Run"
    ].nunique()

    rows.append({
        "ARG": gene,

        "DeepARG_PMA_prevalence":
            r["DeepARG_PMA"],

        "ShortBRED_PMA_prevalence":
            r["ShortBRED_PMA"],

        "RGI_PMA_prevalence":
            r["RGI_PMA"],

        "AMRFinder_PMA_samples":
            pma_runs,

        "AMRFinder_nonPMA_samples":
            non_runs,

        "AMRFinder_PMA_prevalence":
            pma_runs / 32.0,

        # Recurrent support for read-based methods
        "DeepARG_recurrent":
            r["DeepARG_PMA"] >= 0.50,

        "ShortBRED_recurrent":
            r["ShortBRED_PMA"] >= 0.50,

        "RGI_recurrent":
            r["RGI_PMA"] >= 0.50,

        # For assembly-based AMRFinder, any PMA contig detection
        # is treated as complementary evidence, NOT as an
        # equivalent >=50% recurrence criterion.
        "AMRFinder_PMA_detected":
            pma_runs > 0
    })

res = pd.DataFrame(rows)

res["Read_methods_recurrent"] = (
    res[
        [
            "DeepARG_recurrent",
            "ShortBRED_recurrent",
            "RGI_recurrent"
        ]
    ]
    .sum(axis=1)
)

res["Four_method_support"] = (
    res["Read_methods_recurrent"] +
    res["AMRFinder_PMA_detected"].astype(int)
)

res = res.sort_values(
    [
        "Four_method_support",
        "Read_methods_recurrent",
        "DeepARG_PMA_prevalence"
    ],
    ascending=False
)

outfile = OUT / "four_method_candidate_evidence.tsv"

res.to_csv(
    outfile,
    sep="\t",
    index=False
)

print("\n===== FOUR-METHOD CANDIDATE EVIDENCE =====")
print(res.to_string(index=False))

print("\n===== AMRFINDER PMA CANDIDATE DETECTION =====")
print(
    res[
        [
            "ARG",
            "AMRFinder_PMA_samples",
            "AMRFinder_nonPMA_samples",
            "AMRFinder_PMA_prevalence"
        ]
    ].to_string(index=False)
)

print("\nSaved:", outfile)
