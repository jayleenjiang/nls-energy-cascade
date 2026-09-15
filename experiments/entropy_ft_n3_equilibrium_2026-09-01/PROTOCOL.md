# Frozen protocol: n=3 equilibrium bath-heat known-answer audit

Frozen before production.  This experiment does not alter
`flux/NLS_entropy_ft.cpp` or its heat/entropy definitions.

## Scientific question

At equal bath temperatures, test the Cartesian bath-heat accumulator against
the exact equilibrium requirements: zero mean currents and medium entropy,
and a distribution of finite-time medium entropy symmetric about zero.

## Production parameters

- source commit: `3d89659432fac0e512a1cc86fea2b63f8f849762`
- source SHA-256: `9ae5835ed708c8794c8b00ba799b23761482953aaf0ed47cd0b4ba3966d4eaf2`
- mode: `sample_n3`
- chain length: `n=3`
- temperatures: `(T_L,T_R)=(6,6)` and `(10,10)`
- gamma: `0.1` (source constant)
- timestep: `dt=0.0005`
- burn-in: `500`
- block duration: `20`
- batches: `8`, SIMD lanes per batch: `16`, independent streams: `128`
- blocks per stream: `7813`
- blocks per temperature: `1,000,064`
- measured bond: `1` in the code's zero-based right-endpoint convention
- threads: `8`
- seeds: `2026090106` for `T=6`; `2026090110` for `T=10`

The complete blocks CSV is written by the unchanged sampler into a FIFO and
compressed losslessly with Zstandard.  Compression changes storage only, not
the simulation or rows.

## Frozen estimands

For each temperature and block,

`qL_rate = Q_L/20`, `qR_rate = Q_R/20`,
`entropy_rate = Sigma_m/20`, and
`J_E = (Q_L-Q_R)/(2*20)`.

Means are averages of 128 independent stream means.  Standard errors are the
sample standard deviation of stream means divided by `sqrt(128)`.  The 95%
confidence intervals use 5,000 nonparametric bootstrap resamples of whole
streams with analysis seed `2026090191`.  Sigma-from-zero is estimate divided
by that stream-level standard error.

Sign counts report positive, negative, and exact-zero medium-entropy blocks.
Both the naive counting z-score and the stream-level sign-imbalance z-score
are reported.

For the detailed symmetry diagnostic, let `a=Sigma_m/20`.  The plotted value
is `y(a)=log[p(a)/p(-a)]/20`.  Symmetric bins span zero to the empirical 99th
percentile of `|a|`; the bin width follows the Freedman--Diaconis rule, with
the number of bins clamped to `[20,80]`.  The WLS fit uses every bin pair with
at least 200 observations on each side.  Weights are the inverse Poisson
variance `1 / ((1/20^2)*(1/n_plus+1/n_minus))`.  The fit includes an
intercept.  Its 95% CI is a 5,000-resample whole-stream bootstrap using the
same fixed edges and support mask.  No bin or fit window may be selected by
closeness to zero.

## Frozen pass/fail boundary

The heat/entropy definition fails this known-answer audit if, at either
temperature, the mean medium-entropy rate excludes zero at 95% confidence or
is more than 3 stream-level standard errors from zero, or if the equilibrium
symmetry-slope 95% CI excludes zero.  All other observables and raw diagnostics
are reported without changing this rule.

The source-level temperature-factor check is analytic: the boundary update
must have drift `-gamma grad(E)` and noise amplitude `sqrt(2 gamma T)`, so the
continuous bath generator annihilates `exp(-E/T)`.  The finite-step sampler is
still judged by the equilibrium numerical controls above.
