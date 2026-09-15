# Frozen protocol: Gibbs-anchored density-ratio reconstruction for Section 4.1

Protocol date: 2026-09-14  
Protocol version: `section4-density-ratio-v1`

## 1. Question and strict scope

The previous normalized conditional-mixture density failed its predeclared
equilibrium Gibbs/Fokker--Planck calibration.  This new experiment tests a
different representation,

\[
  \rho(X)=\rho_0(X)\exp f_\eta(X)/Z_\eta,
  \qquad \rho_0(X)\propto \exp[-E(X)/5],
\]

where `rho_0` is the exact equal-temperature `(5,5)` Gibbs density and the
network learns only a density ratio.  The reduced state and reference measure
are exactly those audited in `short_chain_section4_2026-09-14`.

No new stochastic simulation is authorized in this phase.  The frozen inputs
are the existing 256 projection-free float64 direct-state trajectories.  The
original density split (32 train, 16 validation, 16 blind-test streams per
physical case) is retained without alteration.

## 2. Independent equilibrium known-answer construction

Within each existing equilibrium split, streams are ranked by

`SHA256("section4-density-ratio-v1|equilibrium-role|split|stream_id")`.

The first half is the pseudo-target sample and the second half is the reference
sample.  Thus equilibrium training/validation/test compare disjoint physical
streams from the same exact Gibbs law.  The exact density ratio is one and the
exact normalized log ratio is zero.  The blind equilibrium roles are written
before any fit and may not be changed.

For the driven problem, all driven streams in a split are target samples and
all equal-temperature streams at the same timestep are reference samples.
No blind stream is read during fitting or selection.

## 3. Frozen model and training candidates

The density-ratio input features are smooth at `I_j=0` and periodic in both
angles.  With action scales `s_j` fixed to the reference-training medians,

`h_j = I_j/(I_j+s_j)`.

The fixed feature vector contains `h_1,h_2,h_3`, all three pairwise products,
`sin(theta_1),cos(theta_1),sin(theta_3),cos(theta_3)`, and
`sin(theta_1-theta_3),cos(theta_1-theta_3)`.  Three candidates are tested:

- `linear`: one linear output on the fixed features;
- `mlp32`: hidden widths `(32,32)`, `tanh` activations;
- `mlp64`: hidden widths `(64,64)`, `tanh` activations.

The final output layer is initialized identically to zero.  Seeds are
`5201,5202,5203`.  Every fit minimizes the balanced logistic density-ratio loss
plus `0.01 * mean((L^dagger rho/rho)^2)` on target-state collocation points.
The pointwise term starts after 20 likelihood-only epochs.  Adam learning rate
is `1e-3`, global gradient norm is clipped at 10, batch sizes are 1024 target,
1024 reference and 256 physics points, with 32 steps per epoch, at most 100
epochs, and validation-loss early stopping after 15 stale epochs.  The initial
zero-ratio model is included as epoch zero in early stopping.

## 4. Frozen equilibrium selection gates

Validation uses fixed, seed-independent subsamples of 100,000 target and
100,000 reference states plus 10,000 reference states for the pointwise
stationarity residual.  A candidate is equilibrium-admissible only when all
three seeds satisfy all of:

1. centered normalized log-ratio RMS `<= 0.05`;
2. held-out target/reference ROC AUC in `[0.48,0.52]`;
3. reference importance-weight ESS fraction `>= 0.50`;
4. median `|L^dagger rho/rho| <= 0.10` and p90 `<= 0.50`;
5. finite values and finite positive normalization estimate.

All equilibrium-admissible architectures advance to the driven validation
stage.  If none pass, the result is `FAIL` and no driven model is trained.
The equilibrium blind split is opened once, after the admissible set is frozen,
and must satisfy the same gates for every seed.

## 5. Frozen driven validation and blind gates

For each equilibrium-admissible architecture, train all three seeds on the
fine-timestep driven/reference training streams.  Validation uses fixed
subsamples and never the blind split.  A candidate is driven-admissible only
when every seed satisfies:

1. reference importance-weight ESS fraction `>= 0.10`;
2. median `|L^dagger rho/rho| <= 0.10` and p90 `<= 0.50` on target states;
3. all 12 predeclared physical moment means agree with target trajectories
   under the equilibrium-calibrated familywise weak-identity threshold;
4. TV distance `<= 0.05` for each 72-bin one-dimensional marginal and the
   fixed `72 x 72` angular marginal;
5. finite validation relative log score and normalization estimate.

The 12 observables are the same as in Section-4 protocol v1:
`I1,I2,I3,M,E,sin(theta1),cos(theta1),sin(theta3),cos(theta3),I1-I3,`
`I1*I2*sin(theta1),I2*I3*sin(theta3)`.

Among admissible candidates select the lowest median relative validation NLL;
ties within across-seed SE go to the smaller architecture.  Freeze the selected
architecture before opening the driven blind split.  Blind gates are identical
and must pass all three seeds.  Importance diagnostics include maximum weight
share and Pareto-tail warning; no clipping or truncation is allowed.

Only after a fine-timestep blind PASS is the selected architecture trained and
tested at `dt=1e-3`.  A numerical-convergence claim additionally requires all
12 coarse/fine moment differences to lie in their combined whole-stream 95%
intervals and pairwise centered blind log-ratio RMS `<=0.10` after fixing the
validation normalization constants.

## 6. Verdicts and claim boundary

- `PASS`: equilibrium validation and blind known-answer gates, driven
  validation and blind gates, all-seed reproducibility, and timestep control
  all pass.  A normalized Gibbs-anchored numerical approximation to the 5D
  NESS may then be reported on the sampled support.
- `PARTIAL`: equilibrium passes but the driven representation, importance
  support, seed, or timestep gate fails.  Only the equilibrium calibration and
  direct trajectory results may be claimed.
- `FAIL`: the equilibrium known-answer stage fails.  NESS remains unopened.

No architecture, threshold, split, feature, physics weight, seed, or fit rule
may be changed after the first non-smoke training output is created.  Software
smoke-test artifacts are stored separately and cannot be used scientifically.

