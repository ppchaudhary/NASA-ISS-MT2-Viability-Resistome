#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

OUT = Path("analysis_ready/amr_concordance")
OUT.mkdir(parents=True, exist_ok=True)

# ============================================================
# Load tables
# ============================================================

deep = pd.read_csv(
    "analysis_ready/deeparg/deeparg_category_long.tsv",
    sep="\t"
)

short = pd.read_csv(
    "analysis_ready/shortbred/shortbred_long.tsv",
    sep="\t"
)

rgi = pd.read_csv(
    "analysis_ready/rgi/rgi_gene_long.tsv",
    sep="\t"
)

amf = pd.read_csv(
    "analysis_ready/amrfinder/amrfinder_hits_long.tsv",
    sep="\t"
)

# ============================================================
# 1. DeepARG categories
# ============================================================

deep_classes = (
    deep[["ARG-category"]]
    .drop_duplicates()
    .sort_values("ARG-category")
)

deep_classes.to_csv(
    OUT / "DeepARG_class_vocabulary.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 2. RGI CARD ARO -> Drug Class lookup
# ============================================================

rgi_lookup = (
    rgi[["ARO", "Drug Class"]]
    .dropna()
    .drop_duplicates()
)

# Some CARD entries may contain multiple drug classes.
# Keep the original CARD string; do not force-collapse yet.

rgi_lookup.to_csv(
    OUT / "CARD_ARO_drugclass_lookup.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 3. Check how many ShortBRED AROs can be assigned CARD class
# ============================================================

short_aro = (
    short[["ARO"]]
    .dropna()
    .drop_duplicates()
)

short_map = short_aro.merge(
    rgi_lookup,
    on="ARO",
    how="left"
)

short_map["Mapped_to_CARD_class"] = (
    short_map["Drug Class"].notna()
)

short_map.to_csv(
    OUT / "ShortBRED_ARO_class_mapping.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 4. RGI drug-class vocabulary
# ============================================================

rgi_classes = (
    rgi[["Drug Class"]]
    .dropna()
    .drop_duplicates()
    .sort_values("Drug Class")
)

rgi_classes.to_csv(
    OUT / "RGI_drugclass_vocabulary.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 5. AMRFinder class vocabulary
# ============================================================

amf_classes = (
    amf[["Class"]]
    .dropna()
    .drop_duplicates()
    .sort_values("Class")
)

amf_classes.to_csv(
    OUT / "AMRFinder_class_vocabulary.tsv",
    sep="\t",
    index=False
)

# ============================================================
# Report
# ============================================================

print("\n===== DEEPARG CLASS VOCABULARY =====")
print(deep_classes.to_string(index=False))

print("\n===== SHORTBRED -> CARD CLASS MAPPING =====")
print("Unique ShortBRED AROs:", len(short_aro))
print(
    "AROs mapped to >=1 CARD drug class:",
    short_map.loc[
        short_map["Mapped_to_CARD_class"],
        "ARO"
    ].nunique()
)

print(
    "AROs without CARD class mapping:",
    short_map.loc[
        ~short_map["Mapped_to_CARD_class"],
        "ARO"
    ].nunique()
)

print("\n===== RGI CARD DRUG-CLASS VOCABULARY =====")
print(rgi_classes.to_string(index=False))

print("\n===== AMRFINDER CLASS VOCABULARY =====")
print(amf_classes.to_string(index=False))

print("\nSaved class-harmonization diagnostic tables to:")
print(OUT)
