# Viability-informed metagenomics of the ISS microbiome and resistome

## Overview

This repository contains processed metagenomic data, sample metadata,
analysis scripts, and supporting inputs for the manuscript:

**Viability-informed metagenomics characterizes microbiome and
resistome restructuring on the International Space Station.**

The analysis examines 64 shotgun metagenomes representing 32 matched
PMA-treated and untreated sample pairs from NASA Microbial Tracking-2.

PMA treatment provides a membrane-integrity-associated DNA proxy.
It does not directly establish microbial viability.

## Repository structure

- `metadata/` — sample metadata and matched-pair information
- `counts/` — processed microbiome and resistome count matrices
- `analysis_inputs/` — additional processed analysis inputs
- `analysis_ready/` — compatibility paths used by analysis scripts
- `scripts/02_microbiome/` — microbiome analyses
- `scripts/03_resistome/` — resistome analyses
- `scripts/04_amr_concordance/` — cross-method AMR comparisons
- `scripts/05_arg_host/` — ARG-associated taxonomic analyses
- `scripts/07_integrated_analysis/` — integrated analyses
- `scripts/08_figures/` — figure-generation scripts
- `scripts/09_supplementary/` — supplementary analyses
- `scripts/10_richness_sensitivity/` — location and flight sensitivity analyses

## Reproducibility

Run analysis scripts from the repository root.

The included processed matrices support downstream analyses without
repeating raw FASTQ processing.

The upstream Bracken genus-construction script additionally requires
64 Bracken `.S.kreport` files in `bracken/`. These upstream files are
not included in the processed-data release.

Some scripts generate intermediate outputs needed by subsequent scripts.
A complete execution-order guide is being finalized.


### Preparing processed inputs

From the repository root, run:

    python3 scripts/00_stage_inputs.py

This stages 51 included processed input files into the
git-ignored `analysis_ready/` directory.

## Verified analyses

The following analyses were independently rerun in an isolated
workspace and numerically compared against reference results:

- Bracken sample QC
- Alpha-diversity metrics and paired tests
- Depth-adjusted microbiome richness
- Depth-adjusted ARG richness
- Microbiome–resistome paired distances
- Coupling depth-sensitivity analysis
- ARG burden coupling analysis

## Scientific interpretation

PMA-associated richness reductions should not be interpreted as
direct measurements of viable microorganisms.

Sequencing depth, sampling location, flight, and other technical
factors may influence observed associations.

The primary depth-adjusted ARG subtype richness model had a
positive-definite Hessian. The secondary ARG category richness model
showed convergence problems and requires caution.

The microbiome–resistome association was weaker after adjustment
for location and flight.


### Additional reproducibility limitations

ANCOM-BC2 differential-abundance analyses were not independently
rerun because the ANCOMBC package was unavailable in the verification
environment.

Microbiome-resistome coupling remained significant after sequencing
depth adjustment (P = 0.00231), but was less conclusive after
additional location and flight adjustment (P = 0.0602).
Within-flight permutation sensitivity analysis gave P = 0.0826.

The repository supports verified downstream analyses, but not a
complete end-to-end reproduction of all upstream analyses.

## Data availability

Raw sequencing accessions and final public archive links will be
verified and added before public release.

## Code availability

The repository will be made public only after applicable NIH/NIAID
clearance and release requirements have been satisfied.

## License

License selection is pending applicable NIH/NIAID review.
