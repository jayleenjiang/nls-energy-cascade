# Frozen protocol: stationary profiles under four boundary couplings

Protocol date: 2026-09-05. The scientific parameters below are frozen before
the first production output is generated.

## Objective

Measure the stationary site-action profile `<I_j>` and bond-angle profile
`<sin(theta_j)>`, with `theta_j = 2(phi_{j+1}-phi_j)`, for the boundary
couplings used in Sections 3.1--3.2. BC3b is retained as a fourth diagnostic
because it differs from BC3 by a phase drift.

## Frozen dynamics

- Base implementation: `NLS_flux_SIMD_fixed.cpp`, archived unchanged in
  `source_archive/`; its required SHA-256 is
  `3919ab963e9d94bcb25ae5ef1c30c2bb032636525db7f6df4d6a318dd41f0656`.
- Production source: `src/NLS_boundary_profiles.cpp`.
- `gamma = 0.1`, maximum `dt = 5e-4`.
- The inherited adaptive step is `max(1e-5,min(dt,1/max_drift))`.
- Corrected boundary noise is used in every case:
  `sigma_I=2 sqrt(2 gamma T I)`, `sigma_phi=sqrt(2 gamma T/I)`.
- The inherited positive-action floor is retained. Every floor projection is
  counted and reported; it is never silently discarded.

Let `M=sum_j I_j`,
`F=2 M I_1-I_1^2+2 I_1 I_2 cos(delta)`, and
`P=2 I_2 sin(delta)`, with the analogous reflected expression at the right
boundary.

| source id | label | action drift | phase drift |
|---:|---|---|---|
| 0 | BC1 | `2 gamma (2T-F)` | `gamma P` |
| 1 | BC2 | `2 gamma (2T-F)` | `0` |
| 3 | BC3 | `2 gamma (2T-I_1^2)` | `0` |
| 2 | BC3b | `2 gamma (2T-I_1^2)` | `gamma P` |

## Frozen run matrix

- Chain lengths: `n in {25,50,100}`.
- Bath temperatures: `(T1,Tn) in {(10,2),(4,6),(6,6)}`.
- Boundary couplings: BC1, BC2, BC3, BC3b.
- Burn-in: 2000 (n=25), 8000 (n=50), 32000 (n=100).
- Measurement duration after burn-in: 2000 for every run.
- Two independently seeded replicates per logical condition.
- Each replicate has 16 SIMD batches x 16 trajectories = 256 trajectories;
  the merged logical run therefore has 512 trajectories.
- Seeds are generated deterministically by `run_profiles.sh` and recorded in
  `RUN_MATRIX.csv` before execution.

Burn-in snapshots are fixed at 25%, 50%, 75%, and 100% of the requested
burn-in. Post-burn cumulative profiles are fixed at 25%, 50%, 75%, and 100%
of the measurement duration. They diagnose stationarity without selecting a
burn-in or fit window after seeing the final profile.

## Statistics and outputs

For each trajectory, the program forms a time average over the measurement
interval. Means and standard errors are then computed across independent
trajectories. Time steps are not treated as independent samples. Invalid
trajectories are excluded from profile statistics and counted explicitly.

Each replicate writes `_profile.csv`, `_burnin_checkpoints.csv`,
`_checkpoints.csv`, `_trajectory_diagnostics.csv`, and `_summary.csv`. The
analysis merges the two predeclared seed replicates into one logical-run CSV
with columns
`j,mean_I,se_mean_I,mean_sin_theta,se_mean_sin_theta` and metadata listing both
seeds, sample count, discarded trajectories, non-finite trajectories, and
projection count.

The central site statistic averages the one central site for odd n and the two
central sites for even n. The central-bond statistic is defined analogously.

## Predeclared BC1 reproduction gate

BC1 at `(10,2)` is run before every other condition. The disk-backed reference
endpoints are:

| n | left | right | reference file |
|---:|---:|---:|---|
| 25 | 0.737907 | 0.165409 | `experiments/lte/test_profile.txt` |
| 50 | 0.514213 | 0.105599 | `experiments/lte/simd_n50_profile.txt` |
| 100 | 0.361686 | 0.0728212 | `experiments/lte/n100_dt25_profile.txt` |

A boundary endpoint passes when its absolute discrepancy is no larger than
`max(3*SE_new, 0.02*abs(reference), 0.002)`. All six endpoints and all three
zero-nonfinite checks must pass. Failure stops the pipeline before BC2/BC3/BC3b
or other temperature pairs are run. The discrepancy is reported rather than
tuned away.

## Diagnostics (not selection rules)

- Burn-in: show all four fixed snapshots and compare the 75% and 100% profiles.
- Sampling convergence: compare cumulative profiles at 50%, 75%, and 100%.
- Equilibrium: report spatial range/slope of `<I_j>` and every
  `<sin(theta_j)>/SE`; the reference is flat action and zero sine.
- Stability: report invalid/discarded counts, projection counts, maximum
  boundary actions, endpoint evolution across checkpoints, and seed-replicate
  differences. Instability or nonreproducibility is a result, not a reason to
  change parameters.
- Scaling: fit the three-point log--log slope of central `<I>` versus n for
  each boundary condition and temperature pair. The hypotheses are slope
  `-1/2` for BC1/BC2 and `0` for BC3/BC3b; these are tests, not imposed fits.

