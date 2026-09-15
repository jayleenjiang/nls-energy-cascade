# R10 coarse-timestep recovery protocol

Frozen after the original coarse replication failed validation and before the
V4 coarse test split was opened.

## Motivation and audit boundary

The original two-moment coarse model failed its predeclared validation gate in
one of three seeds: the right-bond transport moment
`I2*I3*sin(theta3)` had `|z| = 3.962779`, above the frozen
`3.748286557303567` threshold.  All original outputs remain immutable under
`formal_v7_coarse/`; that protocol is recorded as a failure.

The V4 coarse validation split is now development data.  The V4 coarse test
split remains sealed and may be opened once only after this recovery passes
every unchanged validation gate.

## Fixed recovery

Keep the three independently trained R4 neural models unchanged.  Replace the
two-dimensional exponential-family correction by a three-dimensional one with
the observables, in this fixed order,

1. `I1`;
2. `I1*I2*sin(theta1)`;
3. `I2*I3*sin(theta3)`.

Fit the three coefficients on the original coarse train+validation target and
equilibrium-reference streams, then multiply all three coefficients by the
already selected fine-timestep shrink factor 0.625.  There is no factor search,
neural continuation, candidate selection, or test-set access.

## Unchanged gates

On V4 coarse validation require, for every seed: finite weights without
clipping; ESS fraction at least 0.10; relative Fokker--Planck residual absolute
median/p90 at most 0.10/0.50; all 12 moment `|z|` values at most
3.748286557303567; every marginal TV at most 0.05; and maximum pairwise
centered log-ratio RMS at most 0.10.  The exact-zero equilibrium ratio must pass
the same known-answer control.

Only a complete driven and equilibrium validation PASS authorizes the one-shot
V4 coarse test.  A validation failure is reported without another recovery in
this phase.
