# Delivery contents

This compact delivery contains everything needed to inspect and reproduce the
stationary-autocorrelation analysis without duplicating the 342 MB full
trajectory store.

- `report/spectral_gap_autocorr_report.pdf`: visually checked five-page report.
- `analysis/*_autocorrelations.csv`: requested raw autocorrelation and local-rate
  tables for all seven observables.
- `analysis/*_trajectory_correlations.npz`: per-trajectory correlations used by
  the trajectory-level bootstrap.
- `analysis/*_plateaus.csv`, `*_damped_fits.csv`, and `*_diagnostics.csv`:
  machine-readable plateau, frequency, stationarity, sign, and monotonicity
  results.
- `source/NLS_stationary_autocorr.cpp`: exact sampler source.
- `analyze_autocorr.py`, `make_figures.py`, and `audit_results.py`: analysis,
  plotting, and integrity-audit code.
- `PROTOCOL.md`, `FREEZE.md`, and `ANALYSIS_ERRATUM.md`: frozen design and the
  documented analysis-only frequency-identifiability correction.
- `provenance/`: exact production commands, seeds, raw-file SHA-256 manifest,
  integrity audit, and artifact manifest.
- `raw/*/metadata.csv` and `raw/*/run.log`: production metadata and logs.

The 128 production trajectory binaries are retained locally at
`raw/<case>/stream_*.f32`. Their hashes are in
`provenance/raw_manifest.sha256`; they are omitted from the compact ZIP because
the two per-trajectory NPZ files preserve the sufficient inputs for every
reported autocorrelation, bootstrap interval, plateau, and damped fit.
