# Validation report

## Final density

The estimator is an arithmetic ensemble of three separately normalized
Gibbs-anchored density ratios,

`rho_m(x) = exp(-E(x)/5) exp(f_m(x)) / (Z_eq Z_ratio,m)`.

The equilibrium partition is `Z_eq = 514.3433631454566` with Monte Carlo SE
`1.1432018569e-5`.  Per-model ratio normalizers and bootstrap intervals are in
`final_analysis/normalizer_summary.csv`.

## Independent fine test, dt=2.5e-4

All three driven models passed.  ESS fractions are 0.9160--0.9178; relative
Fokker--Planck residual medians are 0.08875--0.09501; maximum absolute moment
z-scores are 1.220--1.443; maximum marginal TVs are 0.02425--0.02440; and the
maximum pairwise centered log-ratio RMS is 0.01273.  The exact equilibrium
control also passed.

## Independent coarse test, dt=1e-3

The original two-moment coarse validation failed before test access because
one right-current moment reached `|z|=3.962779` against the 3.748287 gate.
That failure is preserved.  R10 added the right-current moment as a fixed third
calibration observable, retained the inherited 0.625 shrink factor, passed
validation, and then opened the V4 coarse test once.  On test, ESS fractions
are 0.9087--0.9114; FP medians are 0.08943--0.09310; maximum moment z-scores
are 2.652--2.830; maximum TVs are 0.02405--0.02419; and pairwise RMS is
0.01383.  The exact equilibrium test passed.

## Timestep replication

Fine/coarse ensemble centered log-density RMS is 0.024257 on the fine test
support and 0.024272 on the coarse test support, passing the frozen 0.10 gate.

## Claim boundary

The numerical NESS is validated on sampled stationary support using held-out
moments, marginals, a relative stationary Fokker--Planck residual, independent
seeds, an exact equilibrium control, and timestep replication.  No validation
is available for extreme unsampled tails, and the result is not a closed-form
or computer-assisted mathematical proof of the global stationary density.
