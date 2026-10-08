#!/usr/bin/env python3

import pandas as pd
import numpy as np
import re
from pathlib import Path

OUT = Path("analysis_ready/amr_concordance")
OUT.mkdir(parents=True, exist_ok=True)

# ============================================================
# Shared conservative class vocabulary
# ============================================================

CLASSES = [
    "Aminoglycoside",
    "Beta-lactam",
    "Fosfomycin",
    "Fusidic acid",
    "Glycopeptide",
    "Phenicol",
    "Tetracycline",
    "Diaminopyrimidine",
    "Rifamycin",
    "Sulfonamide"
]

# ============================================================
# Harmonization functions
# ============================================================

def harmonize_deeparg(x):
    x = str(x).lower().strip()

    exact = {
        "aminoglycoside": "Aminoglycoside",
        "beta-lactam": "Beta-lactam",
        "fosfomycin": "Fosfomycin",
        "fusidic-acid": "Fusidic acid",
        "glycopeptide": "Glycopeptide",
        "phenicol": "Phenicol",
        "tetracycline": "Tetracycline",
        "diaminopyrimidine": "Diaminopyrimidine",
        "rifamycin": "Rifamycin",
        "sulfonamide": "Sulfonamide"
    }

    return exact.get(x)


def harmonize_card_piece(x):
    x = str(x).lower().strip()

    if x == "aminoglycoside antibiotic":
        return "Aminoglycoside"

    if x in {
        "penicillin beta-lactam",
        "cephalosporin",
        "carbapenem",
        "monobactam"
    }:
        return "Beta-lactam"

    if x == "phosphonic acid antibiotic":
        return "Fosfomycin"

    if x == "fusidane antibiotic":
        return "Fusidic acid"

    if x == "glycopeptide antibiotic":
        return "Glycopeptide"

    if x == "phenicol antibiotic":
        return "Phenicol"

    if x in {
        "tetracycline antibiotic",
        "glycylcycline"
    }:
        return "Tetracycline"

    if x == "diaminopyrimidine antibiotic":
        return "Diaminopyrimidine"

    if x == "rifamycin antibiotic":
        return "Rifamycin"

    if x == "sulfonamide antibiotic":
        return "Sulfonamide"

    return None


def harmonize_amrfinder(x):
    x = str(x).upper().strip()

    exact = {
        "AMINOGLYCOSIDE": "Aminoglycoside",
        "BETA-LACTAM": "Beta-lactam",
        "FOSFOMYCIN": "Fosfomycin",
        "FUSIDIC ACID": "Fusidic acid",
        "GLYCOPEPTIDE": "Glycopeptide",
        "PHENICOL": "Phenicol",
        "TETRACYCLINE": "Tetracycline",
        "TRIMETHOPRIM": "Diaminopyrimidine"
    }

    return exact.get(x)


# ============================================================
# Metadata denominator
# ============================================================

meta = pd.read_csv(
    "analysis_ready/metadata/master_sample_table.tsv",
    sep="\t"
)

denom = (
    meta.groupby("STATUS")["Run"]
    .nunique()
    .to_dict()
)

print("Sample denominators:", denom)

# ============================================================
# DeepARG
# Detection = category ReadCount > 0
# ============================================================

deep = pd.read_csv(
    "analysis_ready/deeparg/deeparg_category_long.tsv",
    sep="\t"
)

deep["Class"] = deep["ARG-category"].map(harmonize_deeparg)

deep = deep[
    deep["Class"].notna() &
    (deep["ReadCount"] > 0)
].copy()

deep_det = (
    deep[["Run", "STATUS", "Class"]]
    .drop_duplicates()
)

# ============================================================
# RGI
# Detection = All Mapped Reads > 0
# CARD multi-class strings split conservatively
# ============================================================

rgi = pd.read_csv(
    "analysis_ready/rgi/rgi_gene_long.tsv",
    sep="\t"
)

rgi = rgi[
    rgi["All Mapped Reads"] > 0
].copy()

rgi_rows = []

for _, r in rgi.iterrows():

    for piece in str(r["Drug Class"]).split(";"):

        cls = harmonize_card_piece(piece)

        if cls is not None:
            rgi_rows.append({
                "Run": r["Run"],
                "STATUS": r["STATUS"],
                "Class": cls
            })

rgi_det = pd.DataFrame(rgi_rows).drop_duplicates()

# ============================================================
# ShortBRED
# Detection = Count > 0
#
# Map ARO -> CARD Drug Class using RGI annotation lookup.
# ============================================================

short = pd.read_csv(
    "analysis_ready/shortbred/shortbred_long.tsv",
    sep="\t"
)

lookup = pd.read_csv(
    OUT / "CARD_ARO_drugclass_lookup.tsv",
    sep="\t"
)

short = short[
    short["Count"] > 0
].copy()

short = short.merge(
    lookup,
    on="ARO",
    how="left"
)

short_rows = []

for _, r in short.dropna(subset=["Drug Class"]).iterrows():

    for piece in str(r["Drug Class"]).split(";"):

        cls = harmonize_card_piece(piece)

        if cls is not None:
            short_rows.append({
                "Run": r["Run"],
                "STATUS": r["STATUS"],
                "Class": cls
            })

short_det = pd.DataFrame(short_rows).drop_duplicates()

# ============================================================
# AMRFinderPlus
# Detection = assembled AMR hit
# ============================================================

amf = pd.read_csv(
    "analysis_ready/amrfinder/amrfinder_hits_long.tsv",
    sep="\t"
)

amf["Class_harmonized"] = amf["Class"].map(
    harmonize_amrfinder
)

amf_det = (
    amf[
        amf["Class_harmonized"].notna()
    ][
        ["Run", "STATUS", "Class_harmonized"]
    ]
    .rename(
        columns={"Class_harmonized": "Class"}
    )
    .drop_duplicates()
)

# ============================================================
# Combine
# ============================================================

datasets = {
    "DeepARG": deep_det,
    "ShortBRED": short_det,
    "RGI": rgi_det,
    "AMRFinderPlus": amf_det
}

rows = []

for method, dat in datasets.items():

    for status in ["PMA", "non-PMA"]:

        z = dat[
            dat["STATUS"] == status
        ]

        for cls in CLASSES:

            n = z.loc[
                z["Class"] == cls,
                "Run"
            ].nunique()

            rows.append({
                "Method": method,
                "STATUS": status,
                "Class": cls,
                "Detected_samples": n,
                "Total_samples": denom.get(status, 32),
                "Prevalence": n / denom.get(status, 32)
            })

res = pd.DataFrame(rows)

# ============================================================
# PMA-focused concordance summary
# ============================================================

pma = res[
    res["STATUS"] == "PMA"
].copy()

wide = pma.pivot(
    index="Class",
    columns="Method",
    values="Prevalence"
).reset_index()

for m in [
    "DeepARG",
    "ShortBRED",
    "RGI",
    "AMRFinderPlus"
]:
    if m not in wide.columns:
        wide[m] = 0

wide["Methods_PMA_detected"] = (
    wide[
        [
            "DeepARG",
            "ShortBRED",
            "RGI",
            "AMRFinderPlus"
        ]
    ] > 0
).sum(axis=1)

wide["Methods_PMA_50"] = (
    wide[
        [
            "DeepARG",
            "ShortBRED",
            "RGI"
        ]
    ] >= 0.50
).sum(axis=1)

wide["Mean_read_method_PMA_prevalence"] = (
    wide[
        [
            "DeepARG",
            "ShortBRED",
            "RGI"
        ]
    ].mean(axis=1)
)

wide = wide.sort_values(
    [
        "Methods_PMA_detected",
        "Methods_PMA_50",
        "Mean_read_method_PMA_prevalence"
    ],
    ascending=False
)

# ============================================================
# Save
# ============================================================

res.to_csv(
    OUT / "class_prevalence_four_methods.tsv",
    sep="\t",
    index=False
)

wide.to_csv(
    OUT / "class_concordance_PMA.tsv",
    sep="\t",
    index=False
)

print("\n===== PMA CLASS CONCORDANCE =====")
print(wide.to_string(index=False))

print("\n===== PMA PREVALENCE LONG =====")
print(
    pma.sort_values(
        ["Class", "Method"]
    ).to_string(index=False)
)

print(
    "\nSaved:",
    OUT / "class_prevalence_four_methods.tsv"
)

print(
    "Saved:",
    OUT / "class_concordance_PMA.tsv"
)
