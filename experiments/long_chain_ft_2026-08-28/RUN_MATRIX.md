# Long-chain GC run matrix

This table separates accepted evidence from diagnostics.  A row is accepted
only after independent-run, support, paired-GC, time, population, timestep, and
selection-interval controls have been evaluated.  Raw failed runs are retained.

| chain | tilt pair | populations | horizon | current interpretation |
|---:|:---:|:---:|---:|:---|
| 10 | (0.3, 0.7) | 1024, 2048 | 60 | paired/time/population gates and supported `N_c=1024` timestep/selection controls pass |
| 20 | (0.3, 0.7) | 512, 1024 | 60 | strong-pair diagnostic unresolved at the high tilt |
| 20 | (0.4, 0.6) | 512, 1024 | 60 | paired/time/population core gates pass; not promoted to a fully controlled cross-chain claim because the endpoint selection audit fails at n=40 |
| 30 | (0.4, 0.6) | 512, 1024, 2048 | 60 | unresolved: final pair/time/support gates pass, but the frozen 1024-to-2048 `k=0.4` member population gate fails |
| 40 | (0.4, 0.6) | 512, 1024 | 120 | core remedial series and timestep control pass; supported selection `2 -> 1` control fails, so full status is unresolved |

Common replacement-production settings are `T_L=10`, `T_R=2`,
`burnin=500`, `dt=0.0005`, `selection_time=2`, `gauge_shift=0.1`, and
`control_scale=0.5`, with four independent seeds per member.  The gauge makes
the additive observable

```text
Sigma_R = -(1/T_R-1/T_L) Q_R = -0.4 Q_R.
```

It differs from medium entropy only by the endpoint term
`(1/T_L) Delta E` up to the measured first-law residual, and hence targets the
same long-time SCGF under the stated endpoint-moment condition.

The table is a live index, not a claim by itself.  Final numerical status is
read from the generated audit CSV/JSON files and the compiled report.

## Numerical-control matrix

Timestep and selection-interval controls are evaluated at both ends of the
tested chain-length range.  `n=10` uses its resolved pair `(0.3,0.7)`, while
`n=40` uses its resolved pair `(0.4,0.6)`.  The `n=40` baseline is the existing
`N_c=512`, `dt=5e-4`, `selection_time=2` remedial `t=120` production series
and is not rerun.

- `run_numerical_controls.sh`: retained failed `n=10`, `N_c=512` diagnostic;
- `run_numerical_controls_n10_remedial.sh`: frozen supported replacement at
  `n=10`, `N_c=1024`, four seeds per member;
- `run_numerical_controls_n40.sh`: `n=40`, `N_c=512`, `t=120`, four seeds
  per member; its timestep series is the retained control and its
  selection-interval-4 series is retained as an unsupported diagnostic;
- `run_selection_control_n40_remedial.sh`: frozen supported replacement using
  selection interval 1 at `n=40`, `N_c=512`, `t=120`.
