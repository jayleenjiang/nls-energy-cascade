# Frozen protocol: BC1 equal-temperature canonical diagnostic

Frozen on 2026-09-06 before any production output existed.

## Question

Does the systematic nonzero equal-temperature bond-sine profile seen in the inherited single-precision SIMD implementation persist in the double-precision fixed-step canonical implementation?

For BC1 at equal temperatures, the continuum Gibbs measure is even under phase reversal, so every exact stationary expectation `mean(sin(theta_j))` is zero.

## Source audit

The unmodified `flux/NLS_flux_canonical.cpp` is double precision, uses a fixed timestep and Eigen/standard trigonometric functions, and has no adaptive-step rule. It does retain a positive-action projection at `ACTION_FLOOR = 1e-12`; therefore this experiment is a numerical-scheme comparison, not a completely projection-free integration.

The diagnostic copy changes no dynamics. It adds per-trajectory time integrals of every `I_j` and every `sin(theta_j)`, standard errors across independent trajectories, raw per-trajectory profile output, and aggregate projection/near-floor diagnostics. `near_floor_count` counts proposed action values `I <= 10 ACTION_FLOOR`; `minimum_proposed_action` is recorded before projection.

## Frozen runs

- Boundary condition: BC1 canonical only.
- `T1 = Tn = 6`, `gamma = 0.1`.
- Fixed `dt = 5e-4`.
- Measurement duration: 2000.
- 16 batches × 16 trajectories = 256 independent trajectories per seed.
- Two independent seeds per chain length; 512 trajectories after merging.
- `n=25`, burn-in 2000, seeds 2026090625 and 2026090626.
- `n=100`, burn-in 32000, seeds 2026090700 and 2026090701.
- Four OpenMP threads per process. The two seeds of one chain length may run concurrently only on AC power.
- No production parameter may be changed after any production output exists.
- A failed run is preserved and reported; it is not silently retried with altered parameters.

## Frozen statistical rules

For each bond, `z_j = mean(sin(theta_j)) / SE_j`, where the SE is computed across all 512 independent trajectory time averages.

- Pointwise 95% failure: `abs(z_j) > 1.9599639845`.
- Bonferroni simultaneous 95% critical value: `Phi^-1(1 - 0.05/[2(n-1)])`.
- The simultaneous control passes only when the maximum absolute bondwise z-score does not exceed this critical value.
- No bond is excluded and no spatial window is selected.

Action-profile agreement is reported pointwise and through the spatial mean. Before production, “reproduce the SIMD mean action” is operationalized as a relative spatial-mean difference no larger than 5%; the exact difference and full profiles are reported regardless of this descriptive threshold.

## Frozen SIMD references

The comparison uses the already completed BC1 `(6,6)` merged profiles:

- `n=25`: mean action 0.49467016564708577; max `abs(z_sin)=8.678126289279529`; 19/24 pointwise failures.
- `n=100`: mean action 0.24921949993997064; max `abs(z_sin)=12.92852477395248`; 71/99 pointwise failures.

No new SIMD simulation is run.

## Provenance hashes at freeze

- Original canonical source SHA-256: `76f937608280272397a555931b353ba770b06ee87d2f5b0dce08fe1e6bb3727e`
- Instrumented source SHA-256: `c2b897c56a31b0ac38477df31a0900d58d68262f807eba5a38b6cf73e4eac0f5`
- Binary SHA-256: `e950c8b60a2b42e694bed0f19379d485f5b4e33ad160738e3f49fa0c6a944629`

