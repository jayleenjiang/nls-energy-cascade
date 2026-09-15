# Validation report

- Production cases: driven `(2,8)` and equal-temperature `(5,5)`.
- Independent streams: 64 per case; 100,001 saved snapshots per stream.
- Aggregate stationary time: 64,000 per case; maximum lag: 10.
- Raw binary files: 64 per case, 2,800,028 bytes each.
- Non-finite trajectories or values: zero.
- Projection events: 91,674 driven; 97,047 equilibrium.  These are inherited
  from the historical positive-action projection and are a numerical caveat.
- Stationarity: all 14 normalized full-vs-last-three-quarter RMS differences
  are below 0.004, versus the frozen 0.05 gate.
- Plateau bootstrap: all accepted plateaus have 2,000/2,000 valid
  trajectory-bootstrap fits.
- Analysis-only correction: the initial constrained-frequency labels are
  preserved under `analysis_initial_v1/`; `ANALYSIS_ERRATUM.md` documents the
  non-tuned half-cycle/AIC identifiability correction.  Autocorrelations and
  plateau outputs are unchanged.
- Integrity audit: `provenance/integrity_audit.json`.
- Claim boundary: observable-specific relaxation rates are resolved; a common
  spectral gap is not.

