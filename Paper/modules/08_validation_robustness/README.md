# Module 08: numerical validation and robustness

## Manuscript coverage

Local validation appendix supporting the numerical methods.  This module is
kept separate because most of its experiments were not part of the synchronized
advisor-facing draft.

## Experiments and purposes

### Canonical pre-production gates

- Warning-clean and sanitizer smoke builds.
- Thread-count determinism under fixed seeds.
- Equal-temperature zero current.
- Bath-swap current reversal.
- SDE equilibrium profile versus independent Gibbs MCMC reference.
- Purpose: verify implementation, stochastic normalization, and sign
  conventions before interpreting production currents.

### Timestep and stationarity

- Compare dt=1e-3, 5e-4, and selected 2.5e-4 runs.
- Inspect first-half/second-half current blocks.
- Purpose: bound discretization bias and visible finite-time drift.

### Larger chains

- n=50 and n=60 at matched/currently available resolutions.
- Purpose: test whether the four-point n=10--40 exponent immediately drifts
  toward Fourier scaling.

### Bath-temperature robustness

- Production upgrade at (T1,Tn)=(8,4); pilot at (5,1).
- Purpose: test dependence on the primary (10,2) contrast.

### Thermostat-coupling robustness

- gamma=0.05 and gamma=0.2.
- Purpose: test sensitivity to endpoint coupling strength.

## Status

These directories contain real local code and outputs, but they are
VALIDATION_ONLY until discussed with the advisor.  Do not present them as
already approved, and do not merge the full appendix into the synchronized
draft without a section-level decision.
