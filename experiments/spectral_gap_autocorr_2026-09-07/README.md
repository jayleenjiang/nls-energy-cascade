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

Scientific correction recorded before production: Gibbs invariance at equal
temperatures does not make the full Hamiltonian-plus-bath generator
self-adjoint.  Equilibrium monotonicity and sign changes are therefore reported
but are not treated as hard correctness gates.
