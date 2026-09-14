# Final verdict: Gibbs-anchored density-ratio phase

Protocol: `section4-density-ratio-v1`

## Verdict

**PARTIAL: equilibrium calibration passed, but no driven candidate passed the
predeclared validation gates.**

The equilibrium validation stage evaluated three architectures and three seeds.
`linear` and `mlp64` were admissible across all seeds; `mlp32` failed one seed's
median pointwise residual gate.  The frozen equilibrium blind test was then
opened once and all six admissible fits passed the exact zero-density-ratio,
AUC, importance-ESS, finite-support and pointwise-stationarity gates.

Both admissible architectures were trained on the fine-timestep driven data.
All six fits retained excellent equilibrium-reference overlap (importance ESS
fraction 0.954--0.961), seed agreement, and passed every frozen marginal-TV
gate.  Nevertheless, none was driven-admissible:

- median `|L^dagger rho/rho|` was 0.853--0.982, versus the 0.10 gate;
- p90 was 3.955--4.457, versus the 0.50 gate;
- current-sensitive moments `I1*I2*sin(theta1)` and
  `I2*I3*sin(theta3)` missed the held-out trajectory values by 14.6--29.6
  stream-bootstrap standard errors, and additional action-asymmetry moments
  failed in several fits.

Accordingly `FROZEN_DRIVEN_CHOICE.json` records
`NO_ADMISSIBLE_DRIVEN_CANDIDATE`.  The driven blind split was not read and the
coarse-timestep stage was not run.

## Permitted claim

The Gibbs-anchored representation exactly retains the known equilibrium
stationary law and avoids the support/normalization failure of the previous
free density model.  However, the tested linear and MLP density-ratio fits do
not satisfy the driven stationary Fokker--Planck equation or the predeclared
current-sensitive moment gates.  A validated numerical approximation of the
full five-dimensional NESS density therefore remains unresolved.

The marginal agreement may be reported only as a diagnostic; it is not
evidence that the joint stationary density or stationary probability current
has been solved.

