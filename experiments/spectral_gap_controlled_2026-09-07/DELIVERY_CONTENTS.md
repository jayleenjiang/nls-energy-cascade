# Delivery contents

## Human-readable result

- `report/controlled_relaxation_report.pdf`
- `report/controlled_relaxation_report.tex`
- `FINAL_VERDICT.md`
- `VALIDATION_REPORT.md`

## Raw and derived numerical data

- Four local raw directories under `raw/`, one per bath/timestep case.
- 64 stream files per case; 100001 rows by 12 observables per stream.
- `analysis/*_autocorrelations.csv`: mean autocovariance, SE, local rate, and
  bootstrap band at every reported lag for every observable.
- `analysis/*_plateaus.csv`: plateau decisions and effective counts.
- `analysis/early_odd_damped_fits.csv`: raw fit parameters and bootstrap CIs.
- `analysis/early_odd_fit_curves.csv`: plotted data and fitted curves.
- `analysis/previous_vs_controlled_rates.csv`: historical/new comparison.
- `analysis/controlled_summary.json`: complete machine-readable summary.

## Reproduction and integrity

- Production code, analysis, figure, table, and audit scripts are in this
  directory.
- Exact commands and seeds are under `provenance/`.
- `provenance/raw_artifact_manifest.sha256` hashes every raw stream and key
  derived artifact.
- Raw trajectories are intentionally not committed to Git because their total
  size is about 1.1 GiB; all analysis/report artifacts are committed.
