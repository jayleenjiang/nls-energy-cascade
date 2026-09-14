# Final scientific verdict

Protocol: `section4-v1`  
Frozen protocol SHA-256: `0513a546ba80ea0bd8abbbb990c93b74bc416e83732701941429bea391ae1b70`

## Section 4.1: stationary-density reconstruction

**Verdict: FAIL under the predeclared known-answer gate.**

The reduced generator, normalization implementation, direct-state replay,
trajectory-level split, and effective-sample-size checks passed. The complete
equilibrium validation matrix also ran to completion: six model families,
three independent training seeds per family, for 18 fits in total. No family
was admissible for all three seeds. Consequently
`stationary_density/FROZEN_DENSITY_CHOICE.json` records
`NO_ADMISSIBLE_CANDIDATE`.

The decisive failure is pointwise stationarity at the known Gibbs state. Across
all 18 fits, the validation median of
`|L^dagger rho/rho|` lies in `[1.133,2.624]`, against the frozen selection gate
`0.20`; the 90th percentile lies in `[3.705,9.456]`, against the gate `1.00`.
Centered Gibbs log-density RMSE lies in `[0.302,0.595]`. Several fitted Gibbs
slopes are close to one, but good likelihood or slope alone does not establish
the stationary Fokker--Planck equation.

The blind equilibrium density split was therefore not opened, and the NESS
density model was not trained. The permitted paper-level statement is:

> The tested normalized neural density family did not pass an equilibrium
> Gibbs/Fokker--Planck known-answer calibration; hence the five-dimensional
> NESS density and density-derived stabilization mechanism remain unresolved.

It is not permissible to claim that the full five-dimensional nonequilibrium
steady-state density has been numerically solved.

## Section 4.2: stabilization mechanism

**Verdict: not evaluated by the density branch.**

The controlled trajectories are stationary enough for the prespecified direct
weak identities and have adequate effective sample size, but the density gate
failed before the NESS stage. Density slices, tail geometry, and quantitative
high-/low-energy stabilization claims from this model family are therefore not
used as validated evidence.

## Section 4.3: relaxation modes

**Verdict: held-out EDMD-visible modes are supported; a unique spectral gap is
not established.**

The EDMD dictionary and regularization were fitted on 32 trajectories,
selected on 16 trajectories, frozen, and evaluated once on 16 independent
test trajectories. The real-mode blind autocorrelation rates are:

| Bath condition | timestep | blind rate | 95% stream-bootstrap CI |
|---|---:|---:|---:|
| driven `(2,8)` | `1e-3` | -0.8840 | [-0.9707, -0.8103] |
| driven `(2,8)` | `2.5e-4` | -0.9088 | [-1.0139, -0.8178] |
| equilibrium `(5,5)` | `1e-3` | -0.9507 | [-1.0751, -0.8420] |
| equilibrium `(5,5)` | `2.5e-4` | -0.8264 | [-0.9083, -0.7587] |

All four satisfy the frozen validation-to-test residual-growth rule and the
coarse/fine interval-overlap rule. The corresponding slow-mode absolute modal
weights are stable: approximately `0.47--0.50` for `I2`, `0.28--0.30` for
`I1+I3`, `0.20--0.22` for the cosine sum, and `0.106--0.124` for
`cos(theta3)`.

The oscillatory mode is also recovered on held-out data. At the fine timestep,
the blind frequencies are `6.281 [5.388,6.862]` in the driven case and
`5.088 [4.487,5.558]` at equilibrium, consistent with the independently fitted
values `5.500` and `5.264`, respectively. Coarse-timestep intervals technically
exclude zero but are very broad, so the fine-timestep estimates carry the
useful precision.

The pre-existing common-autocorrelation-plateau gate still fails. Therefore the
real rate may be called a *held-out validated EDMD-visible slow mode*, but not
the unique spectral gap of the full generator.

