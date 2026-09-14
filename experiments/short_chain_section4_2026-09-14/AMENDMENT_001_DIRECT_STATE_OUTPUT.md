# Amendment 001: lossless reduced-state output

Date frozen: 2026-09-14, before any density-model fit.

## Trigger

The Stage-1 input audit rejected `equilibrium_dt2p5e-4`: two of 6,400,064
snapshots reconstruct one endpoint action as exactly zero.  The controlled
integrator itself reported strictly positive double-precision actions at those
streams.  The zeros arise because the historical output stores only
single-precision `I1_plus_I3` and `I1_minus_I3`; subtracting those rounded
numbers cancels a small action.

The affected historical rows are fixed audit facts, not exclusions:

- stream 2, row 40629: reconstructed `I3=0`;
- stream 34, row 46751: reconstructed `I1=0`.

No model has been fitted and no validation/test density statistic has been
inspected.

## Frozen repair

Replay all four controlled trajectory cases with the original physical
parameters, integration code, base seeds, stream IDs, burn-in, measurement
time and snapshot interval.  The only change is the writer: it saves

`I1,I2,I3,theta1,theta3`

directly as little-endian float64 values.  The old files remain immutable.
The trajectory-level train/validation/test assignment remains exactly the
assignment frozen in `STREAM_SPLITS.csv`.

## Acceptance gate

Before any model fit, all of the following must hold:

1. 64 streams and 100,001 rows per stream exist in every case;
2. all saved values are finite and all three actions are strictly positive;
3. projection, floor, zero-radius, midpoint-failure and non-finite counts are
   zero;
4. temperatures, timestep and base seed match the original run;
5. observables reconstructed from the direct state agree with every original
   float32 observable to an absolute tolerance of `5e-7`;
6. the two historical cancellation rows are strictly positive in the direct
   output.

Failure stops the density and mode branches.  Passing this gate replaces only
the input representation; it does not change any model, split, hyperparameter,
weak identity, statistical threshold or verdict boundary in `PROTOCOL.md`.

