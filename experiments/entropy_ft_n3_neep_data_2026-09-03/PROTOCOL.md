# Frozen protocol: n=3 consecutive-pair data for a future NEEP test

This task generates and validates trajectory data only.  It does not train,
select, or evaluate a neural entropy-production estimator.

## Immutable simulator

- Validated sampler: `flux/NLS_entropy_ft.cpp`.
- Validated source snapshot commit: `3d89659432fac0e512a1cc86fea2b63f8f849762`.
- Source SHA-256: `9ae5835ed708c8794c8b00ba799b23761482953aaf0ed47cd0b4ba3966d4eaf2`.
- Reused validated binary:
  `experiments/entropy_ft_n3_equilibrium_2026-09-01/bin/entropy_ft_n3_eq`.
- Binary SHA-256:
  `93c4aa6d046f3c35cd0ea1136091fc44ef1b64baf6237e4663ca9e86b995a156`.
- Mode: `sample_n3`.  No C++ source, integrator, force, bath, heat, or entropy
  accumulation code is modified for this experiment.

## Production design

Both cases use `n=3`, `gamma=0.1` (the immutable source constant),
`dt=0.0005`, burn-in `500`, transition duration `delta_t=0.1`, eight SIMD
batches, sixteen lanes per batch, 128 independent streams, 39,063 consecutive
transitions per stream, and eight OpenMP threads.  Thus each case contains
5,000,064 transition pairs and 5,000,192 state snapshots when the unrepeated
initial snapshot of each stream is counted.  Each transition spans exactly
200 integration steps.

Cases and seeds:

| label | `(T_L,T_R)` | seed |
|---|---:|---:|
| driven | `(10,2)` | `2026090310` |
| equilibrium | `(6,6)` | `2026090306` |

The measured action-current bond is `1`, matching the validated n=3 runs.

## Stored rows

The raw archive is the exact unchanged `sample_n3` CSV.  Each row contains
`stream_id`, consecutive interval id, `Q_L`, `Q_R`, `Delta E`,
`Sigma_m=-Q_L/T_L-Q_R/T_R`, action current, first-law residual, and both
endpoint states `(I1,I2,I3,theta1,theta3)`.

A deterministic streaming conversion preserves every raw column and appends
`cos(theta)` and `sin(theta)` for both angles at both endpoints.  It changes no
trajectory value.  Both raw and NEEP-ready CSVs are compressed losslessly with
Zstandard; uncompressed files are never materialized.

## Frozen sanity analysis

For each case report:

1. archive integrity, compressed and decompressed SHA-256, byte size, row and
   stream counts, finite values, contiguous interval ids, endpoint continuity,
   angle range, positive actions, entropy identity, and trigonometric encoding;
2. RMS, mean, standard deviation, maximum absolute value, and selected
   quantiles of `Q_L+Q_R-Delta E` per `delta_t=0.1` interval;
3. decorrelation diagnostics for `I2` and the periodic `theta1` state sequence.

The `I2` correlation is the pooled within-stream normalized autocovariance.
The periodic angle correlation is computed from the centered complex encoding
`exp(i theta1)`, equivalently the summed centered cosine and sine
autocovariances.  For each observable report lag-one correlation, the first
`1/e` crossing, and the Geyer initial-positive-sequence statistical
inefficiency `g=1+2 sum rho(l)`.  The integrated correlation time is
`tau_int=g*delta_t/2`.  Per-stream `1/e` and `tau_int` 5th, 50th, and 95th
percentiles are also reported.  These definitions are frozen before reading
production output.
