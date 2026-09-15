# Observable and controlled-cloning derivation

This note fixes the sign, temperature, and change-of-measure conventions used
by the long-chain experiment.  Heat `Q_r` is positive when bath `r` transfers
energy into the chain.

## Energy and entropy normalization

The Cartesian program advances the physical energy `E=H/2`.  Its equilibrium
density is therefore

```text
rho_eq proportional exp(-E/T) = exp[-H/(2T)],
```

which is the manuscript convention.  The medium entropy delivered to the
reservoirs is

```text
Sigma_m = -beta_L Q_L - beta_R Q_R.
```

No factor of two is missing: if heat is instead quoted in manuscript
Hamiltonian units, `Q_H=2 Q_E`, then the corresponding inverse-temperature
difference is halved and the dimensionless entropy is unchanged.

## Boundary-current gauge

For any constant `g`, define the additive observable

```text
A_g = (-beta_L+g) Q_L + (-beta_R+g) Q_R
    = Sigma_m + g (Q_L+Q_R).
```

The split integrator records

```text
Q_L + Q_R = Delta E + first-law residual.
```

Choosing `g=beta_L=0.1` therefore gives

```text
A_g = -(beta_R-beta_L) Q_R = -0.4 Q_R = Sigma_R,
Sigma_R - Sigma_m = beta_L Delta E + O(first-law residual).
```

The two observables have the same infinite-time SCGF when the required
endpoint-energy exponential moments exist.  At finite time they need not
agree, so their measured difference is retained and tested for decay rather
than set to zero by assumption.

## Exact finite-step controlled proposal

For one Cartesian boundary site, the original Euler bath kernel is

```text
delta = -gamma F(z) dt + sqrt(2 gamma T dt) xi.
```

The controlled sampler proposes the same Gaussian covariance with drift
coefficient `a`:

```text
delta = a F(z) dt + sqrt(2 gamma T dt) xi.
```

For every realized step it accumulates the exact transition-density ratio

```text
log[p_original(delta|z) / p_proposal(delta|z)]
  = (|delta-a F dt|^2 - |delta+gamma F dt|^2)/(4 gamma T dt).
```

If `Delta A_g` is the observed heat functional over a selection interval, the
incremental Feynman--Kac weight is exactly

```text
exp[-k Delta A_g] * p_original(path) / p_proposal(path).
```

Consequently the proposal control changes variance and genealogy, not the
target SCGF.  The frozen production proposal uses `control_scale=0.5`; its
choice was based only on weight and ancestry support diagnostics.  Zero-control
identity and Gaussian-model self-tests are required before production.

## Symmetry target and claim boundary

The estimated object is

```text
psi_n(k) = lim_(t->infinity) t^(-1) log E_NESS exp[-k Sigma_R(t)].
```

The Gallavotti--Cohen target is `psi_n(k)=psi_n(1-k)`.  A finite simulation can
only establish numerical consistency on resolved tilt pairs, horizons,
populations, timesteps, and finite chain lengths.  It cannot prove the full
analytic symmetry or its domain of finiteness.
