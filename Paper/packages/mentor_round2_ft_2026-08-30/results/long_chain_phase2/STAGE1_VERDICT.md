# Phase-II Stage-I verdict

Stage I completed all 40 prospectively specified runs (40 summaries, 40
timeseries, and 40 timing logs).  A finite-value and completion audit found
zero malformed or nonfinite records.  This is an interim checkpoint: the
conditional `n=30` timestep and selection controls remain to be run.

## Frozen-gate results

| chain | baseline residual `D=psi(0.4)-psi(0.6)` | 95% CI | Stage-I result |
|---:|---:|---:|:---|
| 20 | -0.003389 | [-0.012081, 0.005302] | symmetry, support, time, timestep, and selection gates pass |
| 30 | -0.002336 | [-0.007251, 0.002578] | symmetry, support, time, and `N_c=2048 -> 4096` population gates pass; conditional controls authorized |
| 40 | 0.004498 | [-0.002343, 0.011340] | symmetry, support, time, timestep, and selection gates pass |

For `n=20`, all four individual-member control comparisons and both paired
control comparisons pass.  For `n=40`, the same 4/4 member and 2/2 paired
controls pass.  For `n=30`, both `k=0.4` and `k=0.6` member comparisons and
the paired residual pass the `N_c=2048 -> 4096` population gate at the final
horizon; the final pair also passes the support and `t=40 -> 60` gates.

## Current claim boundary

The complete controlled suite is numerically consistent with the
Gallavotti--Cohen relation over the single resolved pair `(0.4,0.6)` at
`n=20` and `n=40`.  The `n=30` result remains provisional until its
prospectively conditional timestep and selection controls finish.  These
tests do not constitute a proof, an all-tilt test, or an infinite-chain
limit.

The secondary normalized baseline residual magnitudes are approximately
4.0%, 4.9%, and 12.9% for `n=20,30,40`, respectively.  In particular, the
absolute frozen gate passes at `n=40`, while the relative uncertainty remains
non-negligible because the SCGF itself is small.  This effect-size caveat must
be retained in any paper-level interpretation.

## Traceability

- prospective specification: `PROTOCOL.md` and `RUN_MATRIX.csv`;
- analysis-grid correction: `ANALYSIS_ERRATUM.md`;
- raw Stage-I output: `raw/`;
- accepted Stage-I analysis: `analysis/final/`;
- conditional execution queue: `run_n30_controls_when_ac.sh`.
