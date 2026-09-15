# Validation report

The analysis protocol was frozen in `PROTOCOL.md` before diagnostic values
were computed.  All new-density quantities use only the exchange-symmetric
Section 4.1 support mask.  Its direct probability is 0.994199 for the driven
case and 0.994034 for equilibrium.

Raw results are in `analysis/`:

- `asymmetry_pointwise.csv` and `asymmetry_angular_tv.csv`;
- `phase_locking.csv`, `phase_histograms.csv`, and
  `density_trajectory_comparisons.csv`;
- `current_balance.csv`;
- `current_balance_full_trajectory_audit.csv` and
  `current_balance_masked_comparison_audit.csv`;
- compact report tables named `report_*.csv`;
- `verdict_frozen_v1.json` and `current_scope_audit_verdict.json`.

The first frozen current gate incorrectly applied the unconditioned identity
E[J12+J23]=0 after support conditioning.  `POST_HOC_CURRENT_SCOPE_AUDIT.md`
documents the issue.  The original failure remains preserved.  The corrected
scope audit checks the exact identity only on full direct trajectories and
checks density-versus-trajectory agreement on the shared mask.  The latter
still fails at 2.306 standard errors, so the correction does not turn the
overall result into a pass.

Reproducibility commands:

```text
./run_analysis.sh
./run_current_scope_audit.sh
./run_report_summary.sh
```

Core analysis SHA-256:
`31c4288612631fdd092dd6688f6267d9cf4fd9a5617a8ec57ac5f0269cdda192`.

Current-scope audit SHA-256:
`d7e48ff24dcfec856e2743199bf24bd05c9a20b7eb4fc26d99b5a8405028c988`.

Training seeds are 8201, 8202, and 8203; the whole-stream bootstrap seed is
20260915042.  Exact model and input hashes are in
`analysis/input_hashes.csv`.
