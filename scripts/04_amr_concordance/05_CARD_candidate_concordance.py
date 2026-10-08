#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

OUT = Path("analysis_ready/amr_concordance")
OUT.mkdir(parents=True, exist_ok=True)

# ============================================================
# Exact DeepARG -> CARD mappings established previously
# ============================================================
# Generic transporter groups are intentionally excluded because
# they cannot be mapped uniquely to a single CARD ARO.

mapping = pd.DataFrame({
    "DeepARG": [
        "RPOB2",
        "BACA",
        "NORA",
        "MSRA",
        "UGD",
        "DFRC",
        "ACRB"
    ],
    "ARO": [
        "ARO_3000501",
        "ARO_3002986",
        "ARO_3000391",
        "ARO_3000251",
        "ARO_3003577",
        "ARO_3002865",
        "ARO_3000216"
    ],
    "Display": [
        "rpoB2",
        "bacA",
        "norA",
        "msrA",
        "ugd",
        "dfrC",
        "acrB"
    ]
})

# ============================================================
# DeepARG persistence results from Figure 4
# ============================================================

d = pd.read_csv(
    "analysis_ready/deeparg/figure4/DeepARG_ARG_group_persistence.tsv",
    sep="\t"
)

print("\n===== DEEPARG COLUMNS =====")
print(list(d.columns))

# ============================================================
# ShortBRED
# ============================================================

sb = pd.read_csv(
    "analysis_ready/shortbred/shortbred_long.tsv",
    sep="\t"
)

print("\n===== SHORTBRED COLUMNS =====")
print(list(sb.columns))

# ============================================================
# RGI
# ============================================================

rgi = pd.read_csv(
    "analysis_ready/rgi/rgi_gene_long.tsv",
    sep="\t"
)

print("\n===== RGI COLUMNS =====")
print(list(rgi.columns))

# ============================================================
# Normalize ARO identifiers
# ============================================================

def norm_aro(x):

    if pd.isna(x):
        return None

    x = str(x).strip()

    if x.startswith("ARO:"):
        return "ARO_" + x.split(":")[-1]

    if x.startswith("ARO_"):
        return x

    if x.isdigit():
        return "ARO_" + x

    return x


sb["ARO_norm"] = sb["ARO"].map(norm_aro)
rgi["ARO_norm"] = rgi["ARO"].map(norm_aro)

print("\n===== EXAMPLE NORMALIZED ARO IDs =====")

print(
    "ShortBRED:",
    sb["ARO_norm"]
    .dropna()
    .drop_duplicates()
    .head(10)
    .tolist()
)

print(
    "RGI:",
    rgi["ARO_norm"]
    .dropna()
    .drop_duplicates()
    .head(10)
    .tolist()
)

# ============================================================
# Detection criteria
# ============================================================
# ShortBRED:
# Count > 0 indicates detection.
#
# RGI:
# All Mapped Reads > 0 indicates detection.

sb["Detected"] = (
    pd.to_numeric(
        sb["Count"],
        errors="coerce"
    )
    .fillna(0)
    > 0
)

rgi["Detected"] = (
    pd.to_numeric(
        rgi["All Mapped Reads"],
        errors="coerce"
    )
    .fillna(0)
    > 0
)

# ============================================================
# Build concordance table
# ============================================================

rows = []

for _, m in mapping.iterrows():

    # --------------------------------------------------------
    # DeepARG
    # --------------------------------------------------------

    dz = d[
        d["ARG_group"] == m["DeepARG"]
    ]

    if len(dz) == 1:

        deep_pma = float(
            dz.iloc[0]["Prevalence_PMA"]
        )

        deep_non = float(
            dz.iloc[0]["Prevalence_non-PMA"]
        )

    else:

        deep_pma = float("nan")
        deep_non = float("nan")

    # --------------------------------------------------------
    # ShortBRED
    # --------------------------------------------------------

    sz = sb[
        (sb["ARO_norm"] == m["ARO"]) &
        (sb["Detected"])
    ]

    sb_pma = (
        sz.loc[
            sz["STATUS"] == "PMA",
            "Run"
        ].nunique()
        / 32.0
    )

    sb_non = (
        sz.loc[
            sz["STATUS"] == "non-PMA",
            "Run"
        ].nunique()
        / 32.0
    )

    # --------------------------------------------------------
    # RGI
    # --------------------------------------------------------

    rz = rgi[
        (rgi["ARO_norm"] == m["ARO"]) &
        (rgi["Detected"])
    ]

    rgi_pma = (
        rz.loc[
            rz["STATUS"] == "PMA",
            "Run"
        ].nunique()
        / 32.0
    )

    rgi_non = (
        rz.loc[
            rz["STATUS"] == "non-PMA",
            "Run"
        ].nunique()
        / 32.0
    )

    # --------------------------------------------------------
    # Save row
    # --------------------------------------------------------

    rows.append({
        "ARG": m["Display"],
        "DeepARG_ID": m["DeepARG"],
        "ARO": m["ARO"],

        "DeepARG_PMA": deep_pma,
        "DeepARG_nonPMA": deep_non,

        "ShortBRED_PMA": sb_pma,
        "ShortBRED_nonPMA": sb_non,

        "RGI_PMA": rgi_pma,
        "RGI_nonPMA": rgi_non
    })


res = pd.DataFrame(rows)

# ============================================================
# Number of methods showing recurrent PMA detection
# Threshold is descriptive: >=50% PMA prevalence.
# ============================================================

res["Methods_PMA50"] = (
    (res["DeepARG_PMA"] >= 0.50).astype(int)
    +
    (res["ShortBRED_PMA"] >= 0.50).astype(int)
    +
    (res["RGI_PMA"] >= 0.50).astype(int)
)

# Also calculate mean PMA prevalence across the three methods
res["Mean_PMA_prevalence"] = res[
    [
        "DeepARG_PMA",
        "ShortBRED_PMA",
        "RGI_PMA"
    ]
].mean(axis=1)

res = res.sort_values(
    [
        "Methods_PMA50",
        "Mean_PMA_prevalence"
    ],
    ascending=False
)

# ============================================================
# Save
# ============================================================

outfile = (
    OUT /
    "CARD_candidate_concordance.tsv"
)

res.to_csv(
    outfile,
    sep="\t",
    index=False
)

# ============================================================
# Console output
# ============================================================

print("\n===== CARD CANDIDATE CONCORDANCE =====")

print(
    res.to_string(index=False)
)

print("\n===== PMA >=50% SUPPORT =====")

print(
    res[
        [
            "ARG",
            "DeepARG_PMA",
            "ShortBRED_PMA",
            "RGI_PMA",
            "Methods_PMA50",
            "Mean_PMA_prevalence"
        ]
    ].to_string(index=False)
)

print(
    "\nSaved:",
    outfile
)

