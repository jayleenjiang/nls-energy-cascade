# Frozen protocol: rebuilding the three-mode stabilization diagnostics

Frozen before any diagnostic value from the new or legacy density is computed.
This is analysis only.  No trajectory generation, density training, model
selection, or calibration is permitted.

## Scientific case and inputs

- driven case: `n=3`, `(T1,T3)=(2,8)`, `gamma=0.1`;
- known-answer control: `(T1,T3)=(5,5)`;
- state and density measure:
  `(I1,I2,I3,theta1,theta3)` with respect to
  `dI1 dI2 dI3 dtheta1 dtheta3`;
- primary density: the fine-step R11 ensemble of seeds `8201,8202,8203`;
- primary held-out data: the already-open V6 fine validation and test streams;
- legacy comparison: `KDE/4:15_NN/h5_files/final.keras` and
  `KDE/4:15_NN/h5_files_eq/final.keras`.

The new ensemble is the arithmetic mean of the three separately normalized
density components.  Density expectations are evaluated by self-normalized
importance sampling from the independent equilibrium-reference test streams.
The old density has no absolute normalization, so its expectations use the
same proposal and self-normalization; its effective sample size is reported.

## Support mask

Section 4.1 operationally validated the density on sampled stationary support
and used the central 0.1--99.9% action ranges for its marginal diagnostics.  We
make that operational support explicit here.  For each physical case, the mask
is fixed from its validation target streams before test evaluation:

- `I1` and `I3` share the 0.001 and 0.999 quantiles of their pooled validation
  values, making the mask exactly invariant under endpoint exchange;
- `I2` uses its own validation 0.001 and 0.999 quantiles;
- both angles retain the complete periodic domain `[-pi,pi)`.

All pointwise and integrated density diagnostics use only this mask.  We report
its trajectory probability and each density's probability mass.  No claim is
made outside it.

## 1. Endpoint-exchange asymmetry

The endpoint exchange is
`sigma(I1,I2,I3,theta1,theta3)=(I3,I2,I1,theta3,theta1)`.
For masked test states, define

`A_sigma(x)=|rho(x)-rho(sigma x)|/[rho(x)+rho(sigma x)]`.

Report its median, 90th percentile, maximum, and stationary-sample mean for each
new seed, the new arithmetic ensemble, and the legacy density.  Test states
are the deterministic union of the held-out streams; no point is selected by
its asymmetry value.  Because those points are themselves stationary draws,
they are not weighted by a second factor of the density.

The integrated angular asymmetry is the total variation distance

`TV_sigma = 0.5 sum_ab |P_ab-P_ba|`

on a fixed `60 x 60` equal-width grid over the full angular torus.  Report the
new density, legacy density, and direct held-out trajectory estimate.  The
same estimators and grid are applied to the equilibrium control.

The driven density asymmetry is reportable only when every new seed has a
larger value than the upper 95% whole-stream-bootstrap limit of the equilibrium
floor and the driven ensemble lower 95% limit exceeds three times that floor.
If the pointwise equilibrium value is at floating-point zero, its multiple is
reported as `Inf`/`N/A`, not replaced by an arbitrary epsilon.

## 2. Conditional phase locking

The exact Hamiltonian used by the new density defines low, middle, and high
energy sectors.  Their boundaries are the 1/3 and 2/3 quantiles of the masked
validation target energy and are frozen before test analysis.

Within each sector and for each angle, report:

- the mode of a fixed 72-bin periodic histogram;
- circular mean `arg E exp(i theta)`;
- circular resultant `R=|E exp(i theta)|`;
- probability within `pi/6` of the high-energy stable branch
  `theta0=pi-asin[(sqrt(4 gamma^2+3)-gamma)/(2(1+gamma^2))]`;
- for `theta3`, probability within `pi/6` of the low-energy reference
  `-2 pi/3`.

Also report the joint probability that both angles lie within `pi/6` of
`theta0`.  Compare each density estimate with the independent held-out
trajectory estimate.  A density/trajectory match requires overlapping 95%
whole-stream-bootstrap intervals; failures remain failures.  New-density seed
bands are the minimum and maximum of the three seed estimates.

## 3. Current balance

Use the state observables

`J12=I1*I2*sin(theta1)`, `J23=I2*I3*sin(theta3)`.

Stationarity of the unforced middle mode requires

`E[J12+J23]=0`.

Report `E[J12]`, `E[J23]`, and their sum from the new density, the legacy
density, and direct held-out trajectories.  Use 1000 whole-stream bootstrap
replicates.  The current-balance gate passes only if the new-density and direct
95% intervals for the sum each contain zero and their difference has
`|z|<=1.96`, where the denominator is the bootstrap standard deviation of the
difference.

## Bootstrap and reproducibility

- whole-stream bootstrap replicates: `1000`;
- bootstrap seed: `20260915042`;
- no per-snapshot IID uncertainty is reported;
- trajectory and proposal streams are resampled independently;
- histogram grids, support, energy boundaries, and all gates remain fixed;
- every reported table is written as CSV, with a JSON verdict and exact
  commands and hashes retained.

## Interpretation gate

An asymmetry claim requires a clear excess over the measured equilibrium
floor and consistency across all three new-density seeds.  A stabilization
statistic is not called density-validated when its density and direct
trajectory bootstrap intervals do not overlap.  The legacy result is retained
for comparison but cannot override either gate.
