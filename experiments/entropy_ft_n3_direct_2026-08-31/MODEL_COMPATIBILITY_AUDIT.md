# Section-4 stationary-density compatibility audit

This audit was completed before the `n=3`, `(10,2)` direct production run.

## Saved models

The Section-4 notebook is
`KDE/4:15_NN/FKE_5d_NLS.ipynb`.  Its executed parameter cell sets
`gamma=0.1`, `T1=2.0`, and `T3=8.0`.  The notebook saves its trained
nonequilibrium model as:

```
KDE/4:15_NN/h5_files/final.keras
SHA-256 906f5e606cf4631c9c3234ee5e09a6b8d19a5d18ed0b7897638d0e95bc07d536
```

The separately loaded equilibrium model is:

```
KDE/4:15_NN/h5_files_eq/final.keras
SHA-256 d9b2294b789d3cea2e753820c8164eaf1191a070fbaf16d43345cc7321ced163
```

The manuscript also states that the nonequilibrium calculation uses
`T1=2,T3=8`.  Neither file is a stationary density for `T1=10,T3=2`.

## Recorded validation accuracy

For the `2,8` nonequilibrium fit, the executed notebook reports:

- scaled log-density RMSE: `0.409876`;
- log-density scaling standard deviation: `3.7679`, so the corresponding
  unscaled log-density RMSE is approximately `1.544`;
- mean absolute normalized FP residual `|L^dagger rho|/rho`: `0.3733`;
- median: `0.2303`; 90th percentile: `0.8784`.

For the separate equilibrium model, the published checks are only fixed-action
angular slices.  Mean relative errors are `1.47%`, `2.62%`, and `6.40%` at
`I=0.5,1,2`.  These slices do not establish the accuracy of endpoint
log-density differences over the full five-dimensional nonequilibrium
support.

## Frozen decision

Using the `2,8` network on `10,2` endpoints would combine a parameter mismatch
with an uncontrolled endpoint log-density error.  It is therefore excluded
from the total-entropy calculation.  This is a pre-data compatibility failure,
not a judgment based on whether an FT slope looks favorable.
