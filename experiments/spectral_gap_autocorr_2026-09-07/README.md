# Three-mode stationary-autocorrelation rate experiment

This directory replaces the window-dependent conditional-expectation rate
diagnostic with stationary autocorrelations from 64 independent trajectories.

- Frozen design: `PROTOCOL.md`
- Sampler: `source/NLS_stationary_autocorr.cpp`
- Analysis: `analyze_autocorr.py`
- Exact production pipeline: `run_pipeline.sh`
- Raw trajectory observables: `raw/<case>/stream_*.f32`
- Raw per-trajectory autocorrelations and CSV tables: `analysis/`

The experiment estimates rates only and does not read or modify the
eigenfunction network `Q`.

## Result

The run is complete.  The frozen gate finds observable-specific plateaus near
`-1` but no common plateau across the seven observables, in either the driven
or equal-temperature case.  Therefore a generator-wide spectral gap remains
numerically unresolved.  See `FINAL_VERDICT.md`, `VALIDATION_REPORT.md`, and
`report/spectral_gap_autocorr_report.pdf`.

The 128 full trajectory files remain local under `raw/` and are intentionally
not tracked by Git.  Their complete SHA-256 manifest is
`provenance/raw_manifest.sha256`; the requested raw autocorrelation tables are
the two `analysis/*_autocorrelations.csv` files.

Scientific correction recorded before production: Gibbs invariance at equal
temperatures does not make the full Hamiltonian-plus-bath generator
self-adjoint.  Equilibrium monotonicity and sign changes are therefore reported
but are not treated as hard correctness gates.
