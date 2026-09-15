# Frozen analysis-only sector follow-up

Frozen before fitting on 7 September 2026. This follow-up starts no simulator
and changes no saved correlations.

## Inputs

- `analysis/driven_T2_T8_autocorrelations.csv`, SHA-256
  `a7c17121d6d0cc77a89162e31c5caa361e14f12a4d8698163f7e4e8ae1164c5e`.
- `analysis/equilibrium_T5_T5_autocorrelations.csv`, SHA-256
  `c36dee1b2da89d8a193cce009ca9ca5775a357ef555c5d7fb3cba9854eb60223`.
- `analysis/driven_T2_T8_trajectory_correlations.npz`, SHA-256
  `9bb7e2c8f9a40a373867ee0bb18d7480178ae884c94d21438a0b43b35c7fff81`.
- `analysis/equilibrium_T5_T5_trajectory_correlations.npz`, SHA-256
  `367cc1e99a2f22a353aebd49e28f212bbdcdebc534e630713a0cc50c079a31f8`.

The NPZ files contain the 64 independent-trajectory correlations used for all
bootstrap intervals. No overlapping time origin is treated as an independent
bootstrap unit.

## Symmetry-label correction

Under endpoint exchange `1 <-> 3`, `theta1 <-> theta3`. Consequently
`sin(theta1)` alone is not an odd eigenobservable and `cos(theta1)` and
`cos(theta3)` alone are not invariant. Exact exchange-sector observables would
use left/right sums and differences. The saved set lacks `sin(theta3)`, so the
requested pair `sin(theta1), I1-I3` is called the **sign-changing pair**, not a
complete odd sector. `I1-I3` is genuinely exchange odd. At unequal bath
temperatures the exchange symmetry is explicitly broken in any case.

The requested five-observable group
`cos(theta1), cos(theta3), cos(theta1-theta3), I2, I1+I3` is retained and
called the **even-candidate group**, with the above limitation stated in every
claim.

## Even-candidate common plateau

Apply the original frozen individual plateau rule unchanged:

- `t >= 1`;
- positive `C` with `C/SE >= 3`;
- at least eight contiguous grid points spanning `t_end/t_start >= 1.5`;
- local-rate range no larger than 20% of the absolute median rate.

A group-level common plateau requires at least three individually accepted
observables, a nonempty intersection of their plateau time intervals, a rate
spread no larger than 20% of the absolute median, and a nonempty intersection
of their 95% trajectory-bootstrap rate intervals. The common estimate, if all
gates pass, is the inverse-variance weighted mean of the accepted rates; its CI
uses a joint 2,000-replicate trajectory bootstrap. Otherwise report `NO COMMON
PLATEAU` and no group rate.

## Sign-changing-pair damped fit

The fit window is fixed at **[0.05, 1.00]** for both observables and both bath
conditions, before inspecting fit results. It brackets the previously reported
first supported negative lags and is not adjusted per observable or case.

Fit by weighted nonlinear least squares

`C(t) = exp(lambda_R t) [A cos(lambda_I t) + B sin(lambda_I t)]`,

with `lambda_R in [-20,0]`, `lambda_I in [0,20]`, and multiple deterministic
starting values. Compare against the nested two-parameter real exponential
`A exp(lambda_R t)` using
`AIC = n log(RSS/n) + 2 p`, with the same trajectory-SE weights.

Uncertainty uses 2,000 bootstrap resamples of the 64 trajectories, seed
`2026090713` plus deterministic case/observable offsets. Report percentile 95%
CIs and acceptance counts. A nonzero oscillatory frequency is resolved only if
all hold:

1. at least 95% of bootstrap fits succeed;
2. the 95% interval for `lambda_I` excludes zero;
3. the damped fit has `Delta AIC <= -10` against the real exponential;
4. the fitted interval contains at least half a cycle,
   `lambda_I >= pi/(1.00-0.05)`.

The first zero of the fitted curve inside the window is calculated from the
fitted phase and compared with the linearly interpolated first zero of the
mean saved correlation. The period alone is not used to predict a first zero,
because the latter also depends on phase.

The two observables agree on a parameter only if their 95% bootstrap intervals
overlap. No fit window, support rule, AIC threshold, or bootstrap setting will
be changed after results are seen.
