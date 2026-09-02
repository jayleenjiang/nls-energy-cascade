# Validation report

## KDE accuracy first

| case | unsupported / total | endpoint log-density RMSE | increment RMSE | lowest-density 1% increment RMSE | increment q99 abs. error | gate |
|---|---:|---:|---:|---:|---:|:---:|
| equilibrium T=6 | 52 / 1,000,064 | 0.273456 | 0.387016 | 1.70194 | 1.38458 | FAIL |
| equilibrium T=10 | 62 / 1,000,064 | 0.270944 | 0.383225 | 1.72094 | 1.38283 | FAIL |

Frozen maxima were 0.15, 0.10, 0.25, and 0.50, respectively.  Metrics above
are calculated only where the frozen KDE is nonzero; unsupported points are
reported and force failure.

On 328,146 held-out driven blocks, two independently trained KDEs have 26
blocks outside their common support.  On the 328,120 supported blocks, their
centered endpoint log densities disagree with RMSE 0.220699 and their system
entropy increments disagree with RMSE 0.311976; the latter rises to 1.46711 in
the lowest-density one percent.  All exceed the frozen maxima.

## Driven direct-sampling result

- KDE-supported endpoint pairs: 1,000,022 / 1,000,064.
- KDE-unsupported endpoint pairs: 42.
- Total-entropy negative count, DFT slope, and IFT: not computed under the
  no-extrapolation rule.
- Medium entropy: 41 negatives, probability
  4.0997376e-5, 95% whole-stream bootstrap interval
  [2.8998144e-5, 5.3996544e-5].
- Medium detailed fit: unavailable; zero symmetric bin pairs have at least 20
  observations per side.
- Medium log IFT estimate: -5.7735605, interval [-8.0873157,-4.8624448],
  exponential-weight ESS 2.2880, maximum weight share 0.61551, maximum
  leave-one-stream change 0.94809.  This diagnostic is unresolved.

## Integrity

The restored driven decompressed CSV SHA-256 exactly matches the accepted
original sample: `4f728b3d0e007d704d90734b0888c00ec05b60f09385c6cfd079f3417d7a088f`.
All three archives contain 1,000,064 ordered finite rows, endpoint continuity
errors are zero at saved precision, and the driven entropy identity error is
2.84217e-14.  No simulation was repeated after the matching restoration.
