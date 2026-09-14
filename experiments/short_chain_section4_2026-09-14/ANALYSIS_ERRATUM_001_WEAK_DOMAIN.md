# Analysis erratum 001: weak-generator domain

Date frozen: 2026-09-14, during Stage 1 and before any density-model fit.

## Problem

The first weak-identity calibration included six observables whose generator
values have non-integrable or infinite-variance endpoint singularities in the
action-angle coordinates. For example,

`L sin(theta1)` contains `-4 gamma T1 sin(theta1)/I1`.

The equilibrium density is nonzero as `I1 -> 0`, so a whole-stream
studentized mean of this quantity does not have a stable finite-variance
standard error. The observed 95% max-|z| calibration threshold of 16.805 is
therefore not an acceptable resolution calibration.

The same issue affects:

- `sin_theta1`, `cos_theta1`, `sin_theta3`, `cos_theta3`;
- `I2_sin_theta1`, `I2_sin_theta3`.

This is an analytic generator-domain error, not a data-selected exclusion.

## Frozen correction

The six singular observables remain in raw output as diagnostic-only rows, but
they do not define the family-wise threshold and cannot pass or fail a weak
identity gate. The simultaneous gate is refrozen on the remaining 18
observables, whose action prefactors cancel all angular `1/I` singularities:

`I1,I2,I3,I1_sq,I2_sq,I3_sq,I1_I2,I2_I3,I1_I3,`
`I1_sin_theta1,I3_sin_theta3,I1_I2_sin_theta1,`
`I1_I2_cos_theta1,I3_I2_sin_theta3,I3_I2_cos_theta3,M,M_sq,E`.

The bootstrap count, seed, stream splits, timestep floor, max-|z| construction
and blind-test rule remain unchanged. The original invalid calibration output
is preserved with the suffix `_initial_invalid_domain`.

Any later density-integral weak test uses exactly this same 18-function gate.

