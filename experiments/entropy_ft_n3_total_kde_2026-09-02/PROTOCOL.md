# Frozen protocol: n=3 finite-time total entropy from a held-out periodic KDE

This protocol is frozen before the driven total-entropy values or any FT fit
are inspected.  It does not change the simulator, the existing driven sample,
or the two equilibrium known-answer samples.

## Scientific question and claim boundary

For the driven `n=3` chain at `(T_L,T_R)=(10,2)`, estimate the stationary
system-entropy boundary term from the saved reduced endpoints and ask whether
ordinary direct sampling resolves the finite-time total-entropy fluctuation
theorem at `t=20`.

The density estimate is an approximation.  A total-entropy slope is a
validated result only if the same estimator passes the exact equilibrium Gibbs
density audit below and the driven data pass the frozen two-sided-support gate.
Otherwise all raw diagnostics are still reported, but the total-entropy FT is
declared not reliably resolved.

## Immutable source data

- Driven sample: `n=3`, `(T_L,T_R)=(10,2)`, `gamma=0.1`, `dt=0.0005`,
  burn-in `500`, block time `20`, `128` independent streams, `7813` blocks per
  stream, `1,000,064` blocks, seed `2026083133`.
- Driven decompressed CSV SHA-256 expected from the completed integrity audit:
  `4f728b3d0e007d704d90734b0888c00ec05b60f09385c6cfd079f3417d7a088f`.
- Equilibrium controls: the accepted `T=6` and `T=10` samples with the same
  `n`, `gamma`, `dt`, burn-in, block time, stream count, and sample count;
  seeds `2026090106` and `2026090110`.
- Simulator source: `flux/NLS_entropy_ft.cpp`, SHA-256
  `9ae5835ed708c8794c8b00ba799b23761482953aaf0ed47cd0b4ba3966d4eaf2`.

If the offloaded driven CSV cannot be materialized safely, it may be restored
by rerunning the exact source, seed, and command into lossless Zstandard.  The
restored decompressed CSV must match the expected SHA-256 above; otherwise the
analysis stops and no substituted sample is used.

## Density coordinates and measure

The saved reduced state is

`x=(I1,I2,I3,theta1,theta3)`,

where actions are positive and angles lie in `[-pi,pi)`.  The KDE is fitted in

`z=(log(I1),log(I2),log(I3),theta1,theta3)`.

The first three kernels are ordinary Gaussian kernels.  Each angular kernel is
a wrapped Gaussian, implemented by periodic convolution on `[-pi,pi)`.  The
density with respect to the physical reduced measure is

`log rho_x = log rho_z - log(I1) - log(I2) - log(I3)`.

The Jacobian term is mandatory in every system-entropy increment.

## Frozen KDE rule

The estimator is a product, binned Gaussian KDE on a
`48 x 48 x 48 x 32 x 32` grid.  The three action axes have constant boundary
conditions; the two angle axes wrap periodically.  Linear interpolation is
used only inside the action grid, with periodic interpolation on angle axes.
There is no density floor and no extrapolation.

For a training set with effective sample size `N_eff`, Scott's
five-dimensional factor is

`s = N_eff^(-1/9)`.

`N_eff=N/g`, where `g` is the largest streamwise statistical inefficiency
estimated by Geyer's initial-positive-sequence rule over
`log(I1),log(I2),log(I3),cos(theta1),sin(theta1),cos(theta3),sin(theta3)`.
Autocovariances are computed within streams and then pooled; streams are never
concatenated.  Using the largest `g` gives one conservative bandwidth factor
for all five coordinates.  The raw endpoint count and every estimated `g` are
reported.

Action bandwidths are `s` times the sample standard deviations of `log(Ij)`.
Angular bandwidths are `s` times the circular standard deviations
`sqrt(-2 log R)`, capped above by the uniform-angle standard deviation
`pi/sqrt(3)`.  No multiplier is tuned.  Histogram-bin variance is subtracted
from the Gaussian smoothing variance (`Delta^2/12`) before convolution; the
analysis fails if a target bandwidth is too small to resolve on the frozen
grid.

Action-grid bounds use the global finite minimum and maximum transformed
coordinates in the applicable data set, padded by four bandwidths computed
from all endpoint samples.  Bounds use no entropy values or FT result.  Angle
bounds are exactly `[-pi,pi)`.

## Independent-stream cross-fitting

No endpoint is evaluated by a KDE trained on its stream.  Fold A trains on
odd-numbered stream endpoints and evaluates both endpoints of even-numbered
streams.  Fold B trains on even streams and evaluates odd streams.  Only block
endpoints train the density, as requested.  The pooled result contains one
held-out density pair for every block.  Fold-specific results are retained to
show estimator variation.

As a driven-data variance control, streams are also divided by `stream_id mod
3`.  Two independent KDEs are trained on groups 0 and 1 and both are evaluated
on group 2.  Their additive-constant-centered endpoint log densities and their
system-entropy increments are compared on exactly the same held-out blocks.
This control passes only if the centered endpoint disagreement RMSE is at most
`0.15`, the increment disagreement RMSE is at most `0.10`, the lowest-density
one-percent increment-disagreement RMSE is at most `0.25`, and its 99th
absolute-error percentile is at most `0.50`.  It measures finite-sample KDE
variation on the driven distribution; it cannot prove absence of common bias.

For every block,

`Delta s_sys = -log rho_x(x_end) + log rho_x(x_start)`

and

`Delta s_tot = Sigma_m + Delta s_sys`.

## Equilibrium accuracy gate (reported before any driven FT result)

At equilibrium the exact reduced Gibbs log density is, up to one additive
constant,

`log rho_exact(x) = -E(x)/T`,

with the same `E=H/2` convention as the simulator.  On held-out endpoints the
unknown constant is removed separately in each cross-fit fold by subtracting
the mean KDE-minus-exact log-density difference.  Report endpoint RMSE, MAE,
median absolute error, and 90/95/99-percentile absolute errors.

The normalization-free block quantity has the exact answer

`Delta s_sys_exact = (E_end-E_start)/T`.

Report the KDE-minus-exact increment RMSE and absolute-error quantiles.  Tail
accuracy is reported in four predeclared bins of the larger endpoint energy:
`[0,80]`, `(80,95]`, `(95,99]`, and `(99,100]` empirical percentiles.

The journal-level KDE gate requires, at both `T=6` and `T=10`:

1. held-out endpoint centered log-density RMSE at most `0.15` overall;
2. system-entropy-increment RMSE at most `0.10` overall;
3. increment RMSE at most `0.25` in the lowest-density one-percent tail; and
4. 99th-percentile absolute increment error at most `0.50`.

These thresholds target errors below ten percent in the central log
probability ratio while explicitly allowing a larger, reported tail error.
They are fixed before the driven result is evaluated.  Failure of any item
means the driven KDE total-entropy FT is descriptive only, not validated.

## Frozen direct-FT analysis

Medium-only and total-entropy diagnostics are produced side by side.

Negative support reports raw negative, positive, and zero counts and the
negative probability with a `2000`-replicate whole-stream bootstrap interval.
Analysis seed: `2026090291`.

Symmetric bins are centered at zero.  The unrounded Freedman--Diaconis width
from all values of the relevant observable is used.  Every bin pair with at
least `20` raw samples on each side enters the WLS fit; no contiguous subwindow
is selected.  The model is

`log[N(+a)/N(-a)] = intercept + slope*a`.

Weights are `1/(1/N_plus+1/N_minus)`.  Report every raw bin count, slope,
intercept, weighted `R^2`, and `2000` whole-stream bootstrap percentile
intervals.  A two-sided detailed-FT fit is resolved only if there are at least
`1000` negative blocks and at least `8` qualifying symmetric bin pairs.  If
the gate fails, a slope may be listed as an unsupported descriptive statistic
only when at least two pairs exist; it is never called an FT validation.

For the integral diagnostic report

`log mean exp(-Delta s)`,

its whole-stream bootstrap interval, exponential-weight ESS and ESS fraction,
and the largest single weight share.  It is numerically resolved only if ESS
is at least `1000`, the largest share is at most `0.01`, and the bootstrap
interval is stable under deleting each stream in turn (maximum leave-one-stream
change in the log mean at most `0.10`).  No plus-four points,
tail extrapolation, KDE extrapolation, or artificial counts are allowed.

The final finite-time FT diagnostic passes only if all of the following hold:
the equilibrium KDE gate passes; the density-variance control on the driven
data passes; the raw two-sided-support gate passes; the
detailed-fit 95% interval contains the reference slope `1`; the integral-fit
95% interval contains `0`; and the integral diagnostic passes its ESS,
single-weight, and leave-one-stream stability gates.  A numerical point
estimate close to the reference is not sufficient.

## Integrity and provenance gates

- Exactly `1,000,064` ordered rows and `7813` rows in every one of `128`
  streams for every input.
- All saved values finite, actions positive, and angles in `[-pi,pi)`.
- Consecutive endpoint continuity and the previously frozen heat/entropy and
  first-law identities must pass their roundoff gates.
- Every raw archive passes `zstd -t`; decompressed SHA-256 values are recorded.
- Because the host had less than 200 MiB free before production restoration,
  the per-block derived total-entropy table is not duplicated on disk.  It is
  deterministically reconstructed from the immutable raw archives and the
  committed analysis code; all raw symmetric-bin counts, fit inputs, KDE
  errors, bandwidths, and gate statistics are retained as small text tables.
- Report the source commit, source/binary/input hashes, seeds, exact production
  or restoration command, analysis command, KDE grid, bandwidths, and all
  failed as well as passed gates.
