# Boundary-profile validation report

Date: 2026-09-06

## Completion and integrity

- Frozen matrix: 72 replicate runs = 4 boundary conditions × 3 temperature pairs × 3 chain lengths × 2 seeds.
- Each replicate contains 16 SIMD batches × 16 trajectories = 256 trajectories; total valid trajectories: 18,432.
- Required outputs present: 72/72 profile files, 72/72 burn-in checkpoint files, 72/72 measurement-checkpoint files, 72/72 trajectory-diagnostic files, and 72/72 summaries.
- Parsed production rows: 4,200 profile rows, 16,800 burn-in checkpoint rows, 16,800 measurement-checkpoint rows, 18,432 trajectory rows, and 72 summary rows.
- Missing files: 0. Non-finite parsed values: 0. Non-finite trajectories: 0. Discarded trajectories: 0.
- Merged logical conditions: 36/36.

## Frozen BC1 reproduction gate

All six endpoint checks passed the predeclared tolerance:

| n | endpoint | existing reference | new mean ± SE |
|---:|:---|---:|---:|
| 25 | left | 0.737907 | 0.7372700 ± 0.0005985 |
| 25 | right | 0.165409 | 0.1652941 ± 0.0001680 |
| 50 | left | 0.514213 | 0.5169271 ± 0.0007108 |
| 50 | right | 0.105599 | 0.1063664 ± 0.0001563 |
| 100 | left | 0.361686 | 0.3582186 ± 0.0008584 |
| 100 | right | 0.0728212 | 0.0720500 ± 0.0001775 |

## Burn-in and stationarity diagnostics

Burn-ins were fixed before production at 2,000, 8,000, and 32,000 for n = 25, 50, and 100. Full profiles at 25%, 50%, 75%, and 100% of burn-in are preserved in every merged `*_burnin_checkpoints.csv` file. Between the final two burn-in snapshots, the largest left/mid/right action displacement among the 36 logical conditions is 2.72 conservative standard errors; these are instantaneous ensemble snapshots rather than time averages. During the subsequent fixed measurement interval, every cumulative endpoint profile changes by less than 0.2% from the 75% checkpoint to the final checkpoint.

The largest pointwise two-seed discrepancy is 4.14 independent-replicate standard errors over the full multiple-comparison family. No BC3 or BC3b trajectory became non-finite, and their endpoint cumulative means did not show continued growth over the measurement interval. Thus the historical finite-time instability was not observed under this frozen protocol; this is not a proof of arbitrarily long-time stability.

## Midpoint-action scaling

Three-point log-log fits give:

| BC | (T1,Tn) | slope ± residual SE | R² |
|:---|:---|---:|---:|
| BC1 | (10,2) | -0.4709 ± 0.0077 | 0.9997 |
| BC1 | (4,6) | -0.4850 ± 0.0120 | 0.9994 |
| BC1 | (6,6) | -0.4728 ± 0.0170 | 0.9987 |
| BC2 | (10,2) | -0.4640 ± 0.0149 | 0.9990 |
| BC2 | (4,6) | -0.4772 ± 0.0116 | 0.9994 |
| BC2 | (6,6) | -0.4794 ± 0.0147 | 0.9991 |
| BC3 | (10,2) | +0.0708 ± 0.0265 | 0.8772 |
| BC3 | (4,6) | +0.0835 ± 0.0233 | 0.9274 |
| BC3 | (6,6) | +0.0750 ± 0.0275 | 0.8814 |
| BC3b | (10,2) | +0.0320 ± 0.0156 | 0.8080 |
| BC3b | (4,6) | +0.0449 ± 0.0094 | 0.9583 |
| BC3b | (6,6) | +0.0455 ± 0.0063 | 0.9810 |

The data support the proposed contrast: with the global M term, BC1/BC2 are close to n^-1/2; without it, BC3/BC3b are close to n^0 over n = 25, 50, 100. Because each exponent uses only three chain lengths, its quoted residual SE is not an asymptotic confidence interval.

## Equal-temperature control

The requested simultaneous zero-sine control fails for all 12 equal-temperature conditions. The largest bondwise |z| is 13.0 for BC1/BC2 and 266 for BC3/BC3b. BC1/BC2 action profiles are comparatively flat at n = 25 but their relative spatial ranges grow to about 9% at n = 100. BC3/BC3b show strong antisymmetric boundary sine profiles.

For BC3/BC3b, equal bath temperatures alone do not guarantee detailed balance because the boundary drift is noncanonical. Nevertheless, the requested gate is failed and must not be reported as passed. For BC1/BC2, the smaller but statistically resolved boundary signal means that a strong equilibrium-profile claim needs a timestep/double-precision or Cartesian cross-check.

## Claim boundary

The profile contrast is suitable as numerical evidence for the boundary-fixed-point mechanism behind the different transport scalings. It is not evidence that every equal-temperature coupling produces a flat, zero-sine equilibrium profile. The inherited single-precision adaptive Euler-Maruyama integrator, positive-action projection, and approximate trigonometric functions are explicit limitations.

## Provenance

- Base `NLS_flux_SIMD_fixed.cpp` SHA-256: `3919ab963e9d94bcb25ae5ef1c30c2bb032636525db7f6df4d6a318dd41f0656`
- Instrumented source SHA-256: `40d13c7547d9c147ca8241ed701e1a1e85540fc0c21db48c4a63a943eccfe3ec`
- Binary SHA-256: `8d3438d3ea893e59216e9784d2104a33ce00d051d4df54c8c19f814fe2c0e6f4`
- Protocol integration commit: `3e6d8e817ee44c8c686f2d80f7c333f3e9d3a973`
- Exact 72 commands and seeds: `COMMANDS.csv`
- File-level manifest: `analysis/FILE_HASHES.csv`
- Final PDF SHA-256 before packaging: `52472156e813e8f2bc238ba9701b2b30355b12349460932f20c0502e181e1922`

