# Controlled n=3 stationary-autocorrelation experiment

This directory contains the controlled-integrator remeasurement requested on
7 September 2026.

## Main result

The fixed, double-precision, projection-free Cartesian integrator does not
produce a common late-time plateau under the unchanged rule.  The unique
spectral gap remains unresolved.  The `I1-I3` damped frequency near 5.35 is
reproduced and is timestep-stable, but other exchange-odd probes have different
frequencies.

## Entry points

- `PROTOCOL.md`: frozen design and gates.
- `source/NLS_stationary_autocorr_controlled.cpp`: production source.
- `run_pipeline.sh`: exact four-case production workflow.
- `analysis/*_plateaus.csv`: every observable's raw accepted/unresolved result.
- `analysis/early_odd_damped_fits.csv`: all early-window fit parameters.
- `analysis/previous_vs_controlled_rates.csv`: requested side-by-side table.
- `report/controlled_relaxation_report.pdf`: six-page report.
- `provenance/raw_artifact_manifest.sha256`: hashes of every raw stream and
  report artifact.

Raw trajectories are local under `raw/<case>/stream_*.f32`.  Each file has
shape `100001 x 12`, little-endian float32, with columns declared in its
case-local `FORMAT.txt`.  The simulated state itself is double precision.
