# Frozen protocol: multi-affinity timestep and estimator-floor audit

Freeze date: 2026-09-11 (America/New_York)

## Questions

This prospective audit extends the completed `(7,5)` timestep study to
`(6.5,5.5)` with `Delta beta = 0.02797202797202797` and `(8,4)` with
`Delta beta = 0.125`. It also uses the already completed `(6,6)` runs to
measure the nonzero floor of the matched-bin slope estimator when the exact
answer is zero.

The existing `dt=5e-4` raw trajectories are immutable baselines. Four new
cases are run: both affinities at `dt=2.5e-4` and `dt=1.25e-4`. This optional
middle timestep is included because four additional raw files require about
2.5 GB while more than 20 GB was free at freeze time, and the 10-core host can
run four two-thread cases concurrently.

## Immutable dynamics and production size

- Production source commit: `1905cf4e606a4a7f4dd8930caa64bd4cc861e9d4`.
- Source SHA-256:
  `98e7f8f5f915c8ce02bd8aa10722025c09fd739184b981961692869c9356c0d3`.
- Binary SHA-256:
  `4c4880d721733897d200f2601690da873f5b226aaa1013df23df2351c1dfe7d1`.
- `n=10`, `gamma=0.1`, burn-in `500`, base block duration `5`.
- 8 batches x 16 lanes = 128 independent streams.
- 31,252 base blocks per stream = 4,000,256 rows per case.
- Stationary time per stream `156260`; total stationary time per case
  `20001280`; two worker threads; middle bond 5.

Fresh production seeds are frozen in `CASES.tsv`: `2026091101` through
`2026091104`. No source, integrator, heat accumulator, compiler option, or
output format is changed.

The immutable `dt=5e-4` rows reuse their original pooled-analysis bootstrap
seeds (`2026092108` for `(6.5,5.5)` and `2026094126` for `(8,4)`) so their
reported point estimates, bootstrap variances, and confidence intervals
exactly reproduce the previously reported baselines.  New cases use distinct
prospectively frozen bootstrap seeds.

## Frozen heat-FR estimator

For each window, `Q=(Q_left-Q_right)/2`. Base blocks are summed within each
stream to `tau=5,10,20,25,40,80,160,320,640`.

At each duration:

- `dQ=sd(Q)/20`;
- raw histogram bins are centred at integer multiples of `dQ`;
- each sign requires at least 10 raw counts;
- select the largest contiguous reliable bin block straddling zero;
- fit `log[p(+Q)/p(-Q)] = a Q + b` by the existing weighted fit;
- at least three symmetric pairs are required for a resolved per-duration fit.

A duration is admissible only when the full-sample fit resolves, the raw
negative count is at least 500, and at least 10 symmetric pairs survive. The
two predeclared selections are all admissible durations and the same set
restricted to `tau>=10`. At least three durations are required.

The primary finite-time model is `a(tau)=a_inf+c/tau`. The independently
reported diagnostic is `a_inf+c/tau+d/tau^2` when at least four durations are
available. The primary model is adequate only if reduced chi-square is at most
2 and its maximum absolute standardised residual is at most 3.

Uncertainty uses exactly 1,000 whole-stream bootstrap replicates. One
multinomial resample of the 128 streams is reused across all durations within
one case. Distinct cases use independent frozen bootstrap seeds. At least 800
jointly valid replicates are required. No plus-four point, fitted density, or
synthetic tail value enters any fit.

For each affinity and selection, the three values of `a_inf/Delta beta` are
also fitted to `r0+m*dt` with inverse-bootstrap-variance weights. The exact
duration set at every timestep is reported. This extrapolation is diagnostic;
the direct finest-step CI and the finite-time adequacy gate take priority.

## Estimator zero-slope baseline

The existing `(6,6)` runs at `dt=5e-4` and `dt=2.5e-4` are analysed with the
same nine durations, matched bins, and joint stream bootstrap. For each
timestep, define fixed inverse-variance weights

`w_tau = 1 / SE_boot[a(tau)]^2`

from the full-sample duration fits, and

`a_eq = sum_tau w_tau a(tau) / sum_tau w_tau`.

Each joint bootstrap replicate uses the same fixed weights and its nine
jointly resampled slopes. Report the percentile CI, bootstrap SE, z-score,
sign pattern, and `a_eq/0.05714285714285716`. The two timestep baselines are
compared by subtracting independent bootstrap replicates. Their baselines are
called statistically different only if that difference CI excludes zero.

The nominal estimator floor is `max_dt |a_eq(dt)|`. It is called the same
order as the residual `(7,5)` finest-step deficit only when their magnitude
ratio lies between one third and three; the exact ratio is always reported.

## Interpretation gate frozen before new output

- An affinity is numerically consistent with the FR at the finest timestep
  only when the primary `1/tau` fit is adequate and the finest-step CI for
  `a_inf/Delta beta` contains one.
- If `(8,4)` still fails the primary adequacy gate, its asymptotic slope is
  unresolved and its intercept is not used as evidence for or against the FR.
- A three-affinity verification statement is allowed only if the completed
  `(7,5)` result and both new affinities satisfy the preceding finest-step
  criterion. The equilibrium estimator floor is quoted as the numerical
  resolution limit.
- If a finest-step CI excludes one, report whether the absolute exclusion in
  slope units exceeds the nominal estimator floor.

No estimator, support threshold, duration, selection, fit form, bootstrap
rule, or gate may be changed after a new production output exists.
