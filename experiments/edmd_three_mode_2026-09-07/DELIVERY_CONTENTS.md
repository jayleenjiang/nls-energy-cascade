# Delivery contents

- `report/edmd_report.pdf` and `report/edmd_report.tex`: readable five-page
  report and its source.
- `FINAL_VERDICT.md`: strict claim boundary.
- `VALIDATION_REPORT.md`: reproducibility, limitations, and fallacy scan.
- `PROTOCOL.md`: choices and gates frozen before spectrum inspection.
- `run_edmd.py`, `summarize_edmd.py`, `make_figures.py`,
  `audit_outputs.py`: complete analysis pipeline.
- `analysis/*_spectra.csv`: every EDMD eigenvalue for all dictionaries, lags,
  and cutoff checks.
- `analysis/*_bootstrap_raw.csv`: every per-replicate matched rate for all 500
  trajectory bootstraps.
- `analysis/matrices/*.csv.gz`: the 500 by 500 D3 `A` matrix and four `B`
  matrices for each of the four datasets, in compressed CSV format.
- `analysis/leading_modes_summary.csv`, `dictionary_convergence.csv`,
  `lag_convergence.csv`, `timestep_convergence.csv`,
  `cutoff_sensitivity.csv`, `known_oscillation_check.csv`, and
  `trivial_mode_check.csv`: all frozen gate inputs and outcomes.
- `analysis/*_observable_modal_weights.csv`: observable projection weights.
- `analysis/reportable_eigenfunction_slices.csv`: full grids used in the
  eigenfunction figure.
- `figures/`: publication PDF and PNG figures.
- `provenance/`: exact commands, integrity audit, and SHA-256 manifest.

The original 1.1 GiB trajectory data are not duplicated.  They remain at
`../spectral_gap_controlled_2026-09-07/raw`, with their own provenance manifest.

