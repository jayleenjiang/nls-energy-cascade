# Frozen protocol: stationary-autocorrelation estimate of the three-mode rate

Frozen before any production output is inspected.  This experiment estimates
relaxation rates only.  It does not train, modify, or evaluate the network
surrogate `Q`.

## Historical target and dynamics

The historical conditional-expectation diagnostic used
`KDE/NLS_backward_Y_train.txt` and reported rates `-1.60`, `-0.934`, and
`-0.68` from different fit windows; an order-two Prony diagnostic gave about
`-1.65`.  The new calculation uses stationary autocorrelations instead.

To make the comparison meaningful, the simulation map is copied from
`cpp/backward/NLS_backward.cpp`:

- three modes with state `(I1,I2,I3,phi1,phi2,phi3)` and reported reduced
  angles `theta1=2(phi1-phi2)`, `theta3=2(phi3-phi2)`;
- Euler--Maruyama, single-precision state, `dt=0.001`;
- the same canonical boundary drift and noise convention;
- `gamma=0.1`;
- the same polynomial trigonometric approximation in the dynamics and the
  same positive-action projection at `1e-14`;
- standard trigonometric functions are used only to measure observables, as in
  the historical `cos(theta1)` measurement.

The new sampler changes only stream management: every trajectory receives a
declared independent xoshiro128++ seed, rather than inheriting an OpenMP-thread
RNG state.  This does not change the Euler--Maruyama update.

Cases:

1. driven: `(T1,T3)=(2,8)`, matching the nonequilibrium three-mode section;
2. equilibrium reference: `(T1,T3)=(5,5)`, matching its Gibbs comparison.

## Frozen sampling

- independent trajectories per case: `64`;
- base seeds: `2026090701` (driven), `2026090702` (equilibrium), with stream
  seed `base_seed + stream_id`;
- initial actions: `I1=I2=I3=1`; initial phases are independent uniform draws
  on `[-pi,pi)` from each stream RNG;
- burn-in: `500` time units;
- stationary measurement duration per trajectory: `1000` time units;
- saved-observable interval: `0.01` (`10` integration steps);
- largest reported lag: `10` time units;
- aggregate stationary time per case: `64,000`, which is `6,400` times the
  largest reported lag and exceeds the requested factor `200`.

The seven saved observables, in order, are:

1. `cos(theta1)`;
2. `sin(theta1)`;
3. `cos(theta3)`;
4. `cos(theta1-theta3)`;
5. `I2`;
6. `I1+I3`;
7. `I1-I3`.

Each trajectory is stored separately as row-major float32 binary data.  This
preserves trajectory-level independence for uncertainty estimation and allows
partial-run recovery without pooling seeds.

## Frozen autocorrelation and lag grid

For each case and observable, one global stationary mean is first estimated
from all saved samples.  Each trajectory then contributes

`C_r(l)=mean_s[(g_r(s)-mean_g)(g_r(s+l)-mean_g)]`.

The reported `C(l)` is the unweighted mean of the 64 trajectory values.  Its
standard error and all confidence intervals use trajectory-to-trajectory
scatter only.

The lag indices are the unique rounded values of
`geomspace(1,1000,121)`, corresponding to `t=0.01,...,10`, plus `t=0`.
At an interior grid point `i`, the local rate is the centred finite difference

`lambda_eff(t_i) = [log C(t_{i+2})-log C(t_{i-2})]/[t_{i+2}-t_{i-2}]`.

It is reported only when both endpoint correlations are positive.  Pointwise
uncertainty bands are percentile intervals from `2000` fixed-seed bootstrap
resamples of the 64 trajectories (`bootstrap seed 2026090711`).

## Frozen plateau rule

The rule is declared before production results are viewed.

- `t_min = 1.0`;
- every included point must have `C(t)>0` and `C(t)/SE[C(t)] >= 3`;
- every included local rate must be finite and negative;
- at least `8` consecutive local-rate grid points;
- interval span `t_end/t_start >= 1.5`;
- rate variation
  `max(lambda_eff)-min(lambda_eff) <= 0.20*abs(median(lambda_eff))`.

Among all qualifying intervals, select the one with the largest log-time span;
ties go to the later interval.  The interval is selected once from the full
sample and is not reselected in bootstrap replicates.  `lambda_R` is the
weighted least-squares slope of `log C(t)` on `t` over that interval, using
`SE[C]/C` as the pointwise scale.  Its uncertainty is the percentile 95% CI
from trajectory bootstrap fits on the frozen interval.  If no interval passes,
report `NO PLATEAU`; no threshold may be relaxed.

A common plateau is declared only if at least three observables have accepted
plateaus, their plateau time intervals have a nonempty common intersection,
their rate estimates differ by at most 20% relative to the median magnitude,
and their individual 95% CIs have a nonempty intersection.  Otherwise no
generator-wide common rate is quoted.

## Damped-oscillation check

For each observable, fit

`C(t)=exp(lambda_R*t)*(A*cos(lambda_I*t)+B*sin(lambda_I*t))`

to the predeclared late-time support: all log-grid points from `t=1` through
the last point with `abs(C)/SE >= 3`, provided there are at least 12 points.
The bounds are `-20 <= lambda_R <= 0` and `0 <= lambda_I <= 20`.  Deterministic
multistart values are used; the solution with the lowest weighted residual is
retained.  Uncertainty comes from 500 trajectory bootstrap fits using seed
`2026090712`.  `lambda_I` is called resolved only if its 95% CI excludes zero
and the optimum is not on a bound.

## Stationarity, effective samples, and equilibrium interpretation

The burn-in/stationarity check repeats the correlation analysis after deleting
the first quarter of every measured trajectory.  For each observable, report
the RMS difference between full and last-three-quarter correlations over
`0<=t<=5`, normalized by `C(0)`.  The predeclared stationarity gate is `<=5%`.

The integrated autocorrelation time uses the initial-positive-sequence rule on
the dense correlation grid through lag 10.  At the largest plateau lag report
the requested effective count

`N_eff = total sampled time / (2*tau_int)`.

At equilibrium, monotonicity and sign changes are reported descriptively.
They are **not** a hard correctness gate: the Gibbs-preserving generator
contains a Hamiltonian Liouville part, which is skew-adjoint in Gibbs
`L2`; the full generator is therefore not generally self-adjoint even at equal
temperatures.  Gibbs invariance alone does not force every autocorrelation to
be positive and monotone.

## Frozen outcomes

- `COMMON PLATEAU`: only if the common-plateau rule passes.
- `OBSERVABLE-SPECIFIC PLATEAU`: one or more observables pass, but the common
  rule fails.
- `NO PLATEAU`: no observable passes or all apparent late-time rates remain
  unsupported/drifting.

The historical values are compared only after this classification.  A null or
inconsistent result is retained unchanged.
