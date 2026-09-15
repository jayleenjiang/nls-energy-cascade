# Gibbs-anchored density-ratio validation report

## Frozen design

- Input: 256 existing projection-free float64 trajectories from the controlled
  three-mode sampler.
- Physical split: 32/16/16 independent streams for train/validation/blind test.
- Equilibrium known answer: disjoint halves of each equilibrium split; exact
  normalized log ratio zero.
- Candidates: `linear`, `mlp32`, `mlp64`; seeds 5201--5203.
- Reference density: `rho_0 proportional to exp(-E/5)`.
- No new stochastic simulation.

## Equilibrium validation and blind result

The all-seed validation-admissible set was `linear, mlp64`.  The `mlp32`
candidate failed because seed 5201 had median pointwise residual 0.1227, above
the frozen 0.10 threshold.  The once-opened blind equilibrium gate passed for
all six admissible fits.  Blind centered log-ratio RMS ranged from 0 to 0.00893,
AUC from 0.5000 to 0.5026, and importance ESS fraction from 0.99992 to 1.00000.

## Driven validation result

| Architecture | Seed | ESS fraction | FP median | FP p90 | Marginals | Moments | Overall |
|---|---:|---:|---:|---:|---|---|---|
| linear | 5201 | 0.9613 | 0.8656 | 3.9998 | pass | fail | fail |
| linear | 5202 | 0.9608 | 0.8702 | 4.0111 | pass | fail | fail |
| linear | 5203 | 0.9598 | 0.8535 | 3.9549 | pass | fail | fail |
| mlp64 | 5201 | 0.9550 | 0.9824 | 4.4574 | pass | fail | fail |
| mlp64 | 5202 | 0.9542 | 0.8930 | 4.0756 | pass | fail | fail |
| mlp64 | 5203 | 0.9617 | 0.8849 | 4.0822 | pass | fail | fail |

The largest marginal TV in each fit was the two-angle marginal and ranged from
0.0377 to 0.0420, below the frozen 0.05 gate.  Pairwise seed centered
log-ratio RMS was 0.0154--0.0169 for `linear` and 0.0525--0.0616 for `mlp64`,
below the 0.10 gate.

The decisive current-sensitive failures were reproducible across every seed:

| Architecture | Seed | Observable | target | ratio estimate | z |
|---|---:|---|---:|---:|---:|
| linear | 5201 | I1 I2 sin(theta1) | 0.23046 | 0.10122 | -29.59 |
| linear | 5201 | I2 I3 sin(theta3) | -0.22980 | -0.10395 | 26.17 |
| linear | 5202 | I1 I2 sin(theta1) | 0.23046 | 0.10274 | -26.88 |
| linear | 5202 | I2 I3 sin(theta3) | -0.22980 | -0.09667 | 24.97 |
| linear | 5203 | I1 I2 sin(theta1) | 0.23046 | 0.09978 | -29.39 |
| linear | 5203 | I2 I3 sin(theta3) | -0.22980 | -0.09377 | 26.68 |
| mlp64 | 5201 | I1 I2 sin(theta1) | 0.23046 | 0.15946 | -14.55 |
| mlp64 | 5201 | I2 I3 sin(theta3) | -0.22980 | -0.12449 | 20.33 |
| mlp64 | 5202 | I1 I2 sin(theta1) | 0.23046 | 0.12236 | -23.47 |
| mlp64 | 5202 | I2 I3 sin(theta3) | -0.22980 | -0.10582 | 24.27 |
| mlp64 | 5203 | I1 I2 sin(theta1) | 0.23046 | 0.12775 | -22.64 |
| mlp64 | 5203 | I2 I3 sin(theta3) | -0.22980 | -0.09571 | 27.05 |

## Integrity and stopping rule

The synchronized Documents tree blocked during the first no-output blind read.
The delivered lossless archive was extracted to a local cache and all 256
files matched the frozen byte counts and SHA-256 hashes.  No scientific data
changed.  The failed I/O attempt, amendment and cache audit are retained.

No driven candidate passed validation, so the protocol required stopping
before driven blind evaluation and before the coarse-timestep control.  No
post-hoc architecture, threshold, fit window or seed was selected.

## Methodological interpretation

The density-ratio fit can reproduce coarse marginals while still missing the
stationary probability-current structure.  The consistent failure of the two
current-sensitive moments, together with the large pointwise adjoint residual,
shows why marginal agreement alone is insufficient for a full NESS-density
claim.  A subsequent protocol must change the optimization/selection objective
rather than weaken these validation gates.

