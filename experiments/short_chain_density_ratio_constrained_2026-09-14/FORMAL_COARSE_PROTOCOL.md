# Coarse-timestep replication protocol

Frozen on 2026-09-14 before any constrained coarse model was trained or any V4
coarse physical statistic was read.

The purpose is a numerical-timestep replication of the final fine model at
`dt=1e-3`, not further method selection.  Use the original independent coarse
trajectories for fitting and the still-sealed V4 coarse validation/test streams
for formal evaluation.

## Neural fit

Train three models from zero with base seeds 8201, 8202, 8203 using exactly the
ten stages in `run_formal_v3_training.sh coarse`, followed by the exact R4
100-epoch transport stage: architecture regular128; regular, hybrid and
transport dictionaries in their frozen order; original schedule, batches,
physics batches, stratification, learning rates, moment boosts and terminal
checkpoint rules.  The coarse fit uses the combined original train and
validation streams and never reads an original or V4 test stream.

Frozen source hashes:

- `run_formal_v3_training.sh`: `a8a5f8fed18425dc29c754b77f82949ebb9fbea86eedaedbf054e35be75e38fa`;
- `train_formal_v3.py`: `9253638e5401dc77380919a9ba2d0089a41125574eeaab6318254e6633e720cc`;
- `train_formal_v3_hybrid.py`: `4c57ba28324af3ec603a47d6646fe0a2ca681b04a66fb2e4eca6138e9630b800`;
- `train_formal_v3_transport.py`: `0468ea29fb0c09abee52868b687c354fe2e0bf530fb05db10a473ba6f1c54e84`;
- `train_recovery_r4.py`: `c5badc2f0f4dbd28fab8164930e378f034a450b8488cada3c83a3d9ac91a0abd`.

## Frozen calibration transfer

After R4, fit the same two exponential-family moments (`I1` and
`I1*I2*sin(theta1)`) on the combined original coarse train+validation streams,
then multiply both fitted coefficients by the fine-selected factor 0.625.
There is no coarse factor search and no neural continuation after calibration.

## Evaluation

Evaluate first on the V4 coarse validation streams with the exact original
gates: ESS >= 0.10; FP median/p90 <= 0.10/0.50; all 12 moment |z| <=
3.748286557303567; all marginal TVs <= 0.05; and pairwise centered log-ratio
RMS <= 0.10.  The exact-zero equilibrium ratio must pass the same known-answer
control.  Only if both complete validation verdicts pass may V4 coarse test be
opened once.  A failure is reported without coarse tuning.
