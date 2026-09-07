# Frozen protocol: EDMD spectrum of the controlled three-mode dynamics

This protocol is frozen before any EDMD matrix or spectrum is computed.  It
reuses the completed controlled-integrator trajectories in
`../spectral_gap_controlled_2026-09-07/raw`; no new simulation is authorized or
needed because the stored columns reconstruct the full reduced state.

## Question and inputs

The question is whether the observable-dependent stationary-autocorrelation
rates are mixtures of several Koopman/generator modes.  The four input cases
are driven `(T1,T3)=(2,8)` and equilibrium `(5,5)`, each at integration
timesteps `1e-3` and `2.5e-4`.  Every case has 64 independent streams, burn-in
500, stationary duration 1000, and snapshot spacing 0.01.

The saved columns reconstruct the reduced state as

```
I1 = (I1_plus_I3 + I1_minus_I3)/2
I2 = I2
I3 = (I1_plus_I3 - I1_minus_I3)/2
theta1 = atan2(sin_theta1, cos_theta1)
theta3 = atan2(sin_theta3, cos_theta3).
```

All state values used by EDMD are therefore measured snapshots, not imputed
states.  The simulator used float64 internally and wrote these observables as
float32.

## Frozen dictionary

Actions use `z_j=(log(max(I_j,1e-12))-mean_j)/sd_j`.  Means and standard
deviations are computed once per case from all saved stationary snapshots and
then held fixed for every dictionary, lag, and bootstrap replicate.

Action functions are monomials `z1^a z2^b z3^c` with total degree at most
`d`.  Angle functions are the real Fourier basis on the two-torus with
`|m1|,|m3|<=M`: a constant plus cosine/sine pairs for one representative of
each `m`/`-m` pair.  The EDMD dictionary is the tensor product of the action
and angle functions.  It always contains exactly one constant.

Three nested dictionaries are fixed:

| label | `(d,M)` | action terms | angle terms | `K` |
|---|---:|---:|---:|---:|
| D1 | `(2,1)` | 10 | 9 | 90 |
| D2 | `(2,2)` | 10 | 25 | 250 |
| D3 | `(3,2)` | 20 | 25 | 500 |

No larger dictionary will be introduced after spectra are viewed.

## Pair quadrature and lags

For computationally reproducible trajectory averaging, each stream contributes
2048 evenly spaced time origins over the common admissible interval.  The same
origins are used at all lags and in all dictionaries.  This gives 131072
origin pairs per case; inferential resampling remains at the 64-stream level.

The frozen physical lags are `tau = 0.10, 0.25, 0.50, 1.00`, corresponding to
10, 25, 50, and 100 saved-snapshot steps.  For each stream and lag,

```
A = mean psi(X_s)^T psi(X_s)
B = mean psi(X_s)^T psi(X_{s+tau})
K_tau = A^+ B.
```

The computation is performed in the retained whitened Gram basis.  The primary
relative eigenvalue cutoff for `A` is fixed at `1e-10` times its largest
eigenvalue.  Sensitivity is reported at `1e-8` and `1e-12`; the primary cutoff
is never selected from spectral symmetry or agreement.  Rates use the principal
complex logarithm `lambda=log(mu)/tau`.

## Frozen checks and gates

The constant mode must satisfy `|mu-1| <= 1e-8` in every full-sample fit.

A nontrivial mode is tracked by minimum distance in the complex-rate plane,
with conjugate pairs represented by the member having nonnegative imaginary
part.  A candidate is dictionary-stable only if it is present in D1, D2, and
D3 at the same lag and the range of real parts is at most
`max(0.15,0.20*|median(real)|)` and the range of absolute imaginary parts is at
most `max(0.50,0.20*median(|imag|))`.  It is lag-stable only if the analogous
bounds hold over all four lags in D3.  It is timestep-stable only if the same
bounds hold between the two integration timesteps for the same bath case.

Statistical intervals use 500 bootstrap replicates that resample all 64
trajectories with replacement.  The dictionary standardisation is not
re-estimated inside a replicate.  Bootstrap modes are matched to the fixed
full-sample D3, `tau=0.5` reference spectrum by minimum complex-rate distance.
A rate is reportable only if it passes dictionary, lag, and timestep gates,
has at least 90% successful bootstrap matches, and its real-part 95% interval
does not overlap that of either neighbouring reportable candidate in real-rate
order.  Failed checks remain in the raw tables.

The known oscillatory check is predeclared as a stable complex mode with
`|Im(lambda)-5.3| <= 0.8`, in addition to the general convergence gates.

## Observable decomposition and eigenfunction slices

At D3 and `tau=0.5`, each requested observable is least-squares projected into
the dictionary.  Modal autocorrelation amplitudes are computed from the EDMD
left/right eigenvectors; the table reports paired-conjugate contributions
normalised by the reconstructed zero-lag covariance.  These are diagnostic
mixture weights and may be signed for a non-normal operator.

For every mode that passes all reportability gates, the eigenfunction is saved
on the `(theta1,theta3)` grid at `I1=I2=I3` equal to 1, 2, and 4.  If no mode
passes, the leading full-sample D3 modes are saved and explicitly labelled
diagnostic/unvalidated.

## Interpretation gate

- If two or more reportable slow modes have materially different observable
  weights, EDMD supports the mode-mixture explanation.
- If one reportable mode dominates all four observables, EDMD does not support
  that explanation.
- If no nontrivial mode passes all convergence gates, the numerical spectrum
  remains unresolved.  Dictionary size, lag, cutoff, and matching tolerances
  will not be changed to rescue the result.

