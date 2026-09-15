# Post-hoc scope audit of the current-balance gate

This note was written after the frozen analysis had produced its first
verdict.  It does not replace or alter `PROTOCOL.md`, its raw masked-current
table, or `analysis/verdict_frozen_v1.json`.

The frozen protocol applied

`E[J12 + J23] = 0`

to expectations conditioned on the Section 4.1 action-support mask.  That is
not an exact stationarity identity: for an invariant density `rho`, the
unconditioned expectation of the middle-mode generator identity is zero, but
conditioning on a strict subset of state space introduces a boundary/selection
term and need not preserve the zero expectation.

The audit therefore makes two additional distinctions without changing the
data or density:

1. the exact zero-current identity is checked on the complete held-out target
   trajectories, with no support conditioning; and
2. on the support mask, the density and trajectory estimates are compared to
   each other, but neither is required individually to equal zero.

The density is not evaluated outside the validated support mask.  Therefore
the first item is a direct-trajectory known-answer check only; it is not a
global density-current estimate.  The second item is the applicable density
reconstruction check.  The original frozen gate remains reported side by
side for auditability.
