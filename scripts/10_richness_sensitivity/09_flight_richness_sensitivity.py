from pathlib import Path
import pandas as pd

out = Path("analysis_ready/integrated_analysis")
source = out / "richness_location_paired_differences.tsv"

df = pd.read_csv(source, sep="\t")

assert len(df) == 64
assert df.groupby("Metric").size().eq(32).all()

results = (
    df.groupby(["Metric", "Flight"])
    .agg(
        N_locations=("LocationCode", "nunique"),
        Median_PMA_minus_nonPMA=("PMA_minus_nonPMA", "median"),
        Mean_PMA_minus_nonPMA=("PMA_minus_nonPMA", "mean"),
        Locations_lower_after_PMA=(
            "PMA_minus_nonPMA",
            lambda x: (x < 0).sum()
        ),
        Locations_higher_after_PMA=(
            "PMA_minus_nonPMA",
            lambda x: (x > 0).sum()
        )
    )
    .reset_index()
)

assert len(results) == 8
assert results["N_locations"].eq(8).all()

path = out / "richness_flight_sensitivity.tsv"
results.to_csv(path, sep="\t", index=False)

print("===== FLIGHT SENSITIVITY RESULTS =====")
print(results.to_string(index=False))

print("\nSaved:", path)
print("===== STEP 61 COMPLETE =====")
