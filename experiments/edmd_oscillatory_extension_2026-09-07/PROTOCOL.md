# Frozen protocol: targeted EDMD completeness test for the n=3 oscillatory mode

This protocol is frozen before any extended-dictionary EDMD matrix or spectrum
is computed.  It is analysis-only and reuses the four controlled trajectory
ensembles in `../spectral_gap_controlled_2026-09-07/raw`.  No simulation or
integrator parameter is changed.

## Scientific question

The first EDMD audit resolved one leading real rate near `-0.95` but did not
recover the independently measured oscillation with imaginary rate near 5.3.
This test asks whether a single, predeclared targeted extension that preserves
raw action differences recovers that complex mode without moving the leading
real rate.

## Inputs and state reconstruction

The four inputs are driven `(T1,T3)=(2,8)` and equilibrium `(5,5)`, each at
integration timesteps `1e-3` and `2.5e-4`.  Every case contains 64 independent
streams, burn-in 500, stationary duration 1000, and 100001 snapshots at spacing
0.01.  The controlled simulator used float64 state, fixed timestep, standard
trigonometry, no floor, and no projection; saved observables are float32.

The reduced state is reconstructed exactly from saved columns:

```
I1 = (I1_plus_I3 + I1_minus_I3)/2
I2 = I2
I3 = (I1_plus_I3 - I1_minus_I3)/2
theta1 = atan2(sin_theta1, cos_theta1)
theta3 = atan2(sin_theta3, cos_theta3).
```

## Frozen targeted dictionaries

The original basis is retained: action monomials in standardized
`log(max(Ij,1e-12))`, tensor-producted with the real Fourier basis in
`(theta1,theta3)`.  The three nested original parts are D1 `(d,M)=(2,1)`, D2
`(2,2)`, and D3 `(3,2)`, with sizes 90, 250, and 500.

Each extended dictionary appends the same 19 raw, unstandardized columns in the
following fixed order:

1. `I1`, `I2`, `I3`;
2. `I1^2`, `I1*I2`, `I1*I3`, `I2^2`, `I2*I3`, `I3^2`;
3. with `dI=I1-I3`,
   `dI`, `dI*cos(theta1)`, `dI*sin(theta1)`, `dI*cos(theta3)`,
   `dI*sin(theta3)`, `dI*cos(theta1-theta3)`,
   `dI*sin(theta1-theta3)`, `dI^2`, `dI*(I1+I3)`, `dI*I2`.

The nominal extended sizes are E1=109, E2=269, and E3=519.  Four targeted
columns are exact linear combinations of the direct raw degree-one/two columns;
they are deliberately retained as explicit requested features.  The already
frozen truncated Gram inverse must remove redundant directions.  Both nominal
and retained dimensions are reported.  No feature will be added or removed
after spectra are viewed.

## Pair quadrature, lags, and aliasing

Each stream contributes the same 2048 evenly spaced origins as the first EDMD
audit, over the common admissible interval.  Thus every case uses 131072 pairs,
with inferential resampling at the 64-stream level.

The short-lag set is `tau={0.02,0.05,0.10}`.  The long-lag set is
`tau={0.25,0.50,1.00}`.  The principal-log Nyquist bounds are reported as
`pi/tau` for every lag.  A complex candidate satisfying
`|Im(lambda)| >= 0.8*pi/tau` is marked aliased and cannot count toward lag
convergence.  The oscillatory completeness gate uses only the unaliased short
set.  Long-lag behavior is reported separately and cannot rescue or reject the
short-lag result.

For comparison with the first audit, the leading-real-rate lag gate remains
the original set `{0.10,0.25,0.50,1.00}`; its primary estimate remains E3 at
`tau=0.50`.

## EDMD estimator and Gram truncation

For every stream and lag,

```
A = mean psi(X_s)^T psi(X_s)
B = mean psi(X_s)^T psi(X_{s+tau})
K_tau = A^+ B.
```

The computation uses the retained whitened Gram basis.  The primary relative
Gram-eigenvalue cutoff remains `1e-10`; fixed sensitivity checks are `1e-8`
and `1e-12`.  Rates use the principal complex logarithm
`lambda=log(mu)/tau`.  The constant mode must satisfy `|mu-1|<=1e-8`.

## Frozen candidate selection and matching

### Oscillatory candidate

In each case independently, the reference is selected from E3 at `tau=0.05`
as the candidate with `4.5 <= Im(lambda) <= 6.0`,
`-20 <= Re(lambda) <= -0.05`, and largest real part.  This is a deterministic
band-and-order rule, not nearest-to-5.3 tuning.  Other fits are matched by the
nearest complex rate, accepted only within
`max(0.50,0.35*|lambda_reference|)`.

The oscillatory candidate passes only if:

- it is present and unaliased in E1, E2, E3 at `tau=0.05`;
- its real-part range is at most `max(0.15,0.20*|median(real)|)` and its
  imaginary-part range is at most `max(0.50,0.20*median(|imag|))` across those
  dictionaries;
- the same bounds hold across all three unaliased short lags in E3;
- the same bounds hold between integration timesteps for the same bath case;
- at least 450/500 trajectory-bootstrap fits match; and
- the bootstrap 95% interval for the imaginary part excludes zero.

### Leading real and neighbouring candidates

The E3, `tau=0.5` reference candidates use the same rule as the first audit:
the first 12 modes in decreasing real part with
`-3<=Re(lambda)<=-0.05`, `0<=Im(lambda)<=12`, and `0<|mu|<=1.05`.
Dictionary, original-lag, timestep, cutoff, 500-bootstrap, reciprocal
timestep matching, and neighbouring-CI separation gates are unchanged from
the first audit.  The previous leading-rate CIs used for the declared movement
test are:

- driven `dt=1e-3`: `[-0.954628,-0.920824]`;
- driven `dt=2.5e-4`: `[-0.980007,-0.941942]`;
- equilibrium `dt=1e-3`: `[-0.972388,-0.932803]`;
- equilibrium `dt=2.5e-4`: `[-0.954830,-0.922965]`.

A new leading point estimate counts as moved only if its new 95% interval has
empty intersection with the corresponding previous interval.

## Bootstrap and observable decomposition

Use 500 replicates, resampling all 64 trajectories with replacement.  Time
origins are never resampled as independent observations.  Standardization is
fixed from the full case and is not re-estimated inside replicates.  Bootstrap
spectra are computed separately at `tau=0.05` for the oscillatory candidate and
at `tau=0.50` for real candidates.

At E3 and `tau=0.5`, project `I2`, `I1+I3`, `cos(theta3)`, and
`cos(theta1)+cos(theta3)` into the dictionary and report signed modal covariance
weights.  Compare the leading-mode weights with the first-audit values
`0.490,0.283,0.123,0.224` and the autocorrelation plateaus
`-0.935,-1.003,-1.092,-1.122`.

If an oscillatory mode passes, save its eigenfunction and the leading real
eigenfunction on the `(theta1,theta3)` grid at `I1=I2=I3` equal to 1, 2, and 4.
If it fails, save its diagnostic slices but label them unvalidated.

## Frozen interpretation gate

- If the oscillation passes and the leading real interval overlaps its previous
  interval, the slow spectral estimate is declared dictionary-converged within
  this targeted test, and the leading rate near `-0.95` may be reported.
- If the oscillation passes but the leading interval does not overlap its
  previous interval, the earlier leading estimate is dictionary-limited; the
  new value must replace it.
- If the oscillation fails, the EDMD spectrum remains incomplete despite this
  targeted extension.  No further feature or lag tuning is permitted in this
  task, and `-0.95` remains only the slowest mode visible to the tested
  dictionaries.

