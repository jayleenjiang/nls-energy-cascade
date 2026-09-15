# Authorization for the one-shot R10 coarse blind test

Frozen after the R10 driven validation and unchanged equilibrium validation
passed, and before any V4 coarse test statistic was read.

The original two-moment coarse validation remains a recorded failure in
`formal_v7_coarse/driven_validation`.  R10 used that validation failure as
development information and added the right-bond transport moment under
`RECOVERY_PROTOCOL_R10_COARSE.md`.  The V4 coarse test split has not been used
for fitting, calibration, selection, or evaluation.

## Validation authorization

- driven R10 verdict SHA-256:
  `988f7d08acf8c47ecab2d99d8e493a5971fcd3fd2a9256ff4544cb4c6d750fbb`;
- driven R10 moment table SHA-256:
  `1ec97b180e002183d3593aecd69f10e5ca6940d2575d664e5b28f10fe9d15bf9`;
- unchanged equilibrium verdict SHA-256:
  `f10aab180e489558f47d7d23b0548150d4f5a7c0c2a13eabf4b8c11d19aeeb72`.

The driven complete verdict is PASS.  Per-model maximum moment `|z|` values
are 2.872554, 2.493684, and 2.592332; relative Fokker--Planck residual medians
are 0.093617, 0.093916, and 0.094729; ESS fractions are 0.911531, 0.912212,
and 0.913943.  Maximum pairwise centered log-ratio RMS is 0.014309.  The
exact-zero equilibrium validation also has a complete PASS.

The following model weight SHA-256 values are frozen:

- seed 8201: `ec10a5a4712e61c45232cac16268fb71969c22f708a47e337f242bfeec548ffe`;
- seed 8202: `e3b753ac849f2d669c86aea729cb32696785676dd09ded2b6ccd50af1e2612cf`;
- seed 8203: `7b9628d2ba5fce8fc01b84d266232da499a7b0f6025ac506140c186f3d699f32`.

The V4 coarse test may now be evaluated exactly once with the unchanged gates
and validation-frozen marginal edges.  No test failure will be followed by
another R10 model change or test-set reuse.
