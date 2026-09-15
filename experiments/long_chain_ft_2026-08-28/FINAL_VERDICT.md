# Final verdict: long-chain Gallavotti--Cohen diagnostic

## Bottom line

The simulation does **not** verify the fluctuation theorem uniformly for
`n=10,20,30,40`.  It provides a fully controlled numerical-consistency result
at `n=10` and positive but incomplete evidence at the longer chains.  This is
not evidence that the theorem is false at the unresolved lengths; it means the
finite rare-event estimator has not passed every frozen convergence and
numerical-sensitivity gate there.

| n | tested pair | final horizon/population | core paired-SCGF audit | extra numerical audit | admissible status |
|---:|:---:|:---:|:---|:---|:---|
| 10 | `(0.3,0.7)` | `t=60`, `N_c=2048` | pass | supported timestep and selection controls pass | numerically consistent with GC over the tested pair |
| 20 | `(0.4,0.6)` | `t=60`, `N_c=1024` | pass | no chain-specific control; the predeclared endpoint selection suite does not pass globally | positive core evidence, not a fully controlled cross-chain result |
| 30 | `(0.4,0.6)` | `t=60`, `N_c=2048` | paired/time/support pass | `k=0.4` population-member gate fails | unresolved population convergence |
| 40 | `(0.4,0.6)` | `t=120`, `N_c=1024` | pass after independent remedial series | timestep passes; supported selection control fails | unresolved selection sensitivity |

The final paired residuals are stored in `final_summary/final_gc_summary.csv`.
The core audit passes 3 of 4 rows; the endpoint numerical audit passes 3 of 4
control rows.  These counts are mechanical summaries of frozen gates, not
formal family-wise hypothesis tests.

## What can be claimed

> For the boundary-driven Cartesian NLS chain at `(T_L,T_R)=(10,2)`, the
> complete frozen rare-event and numerical-control suite is consistent with
> the Gallavotti--Cohen SCGF symmetry at `n=10` over the resolved pair
> `(k,1-k)=(0.3,0.7)`.  The resolved core SCGF pairs at `n=20` and `n=40` are
> also consistent with the symmetry, but the complete long-chain control suite
> is not passed: population convergence remains unresolved at `n=30`, and
> selection-interval sensitivity remains unresolved at `n=40`.

## What cannot be claimed

- a proof of the fluctuation theorem;
- symmetry over every `k` in `[0,1]`;
- a verified statement for all four tested chain lengths;
- an infinite-chain limit;
- failure of the physical FT at `n=30` or `n=40` (the failed gates diagnose
  estimator convergence/sensitivity, not a statistically resolved symmetry
  violation).

## Traceability

- frozen core rows: `FINAL_RESULT_SPEC.csv`;
- fixed gates and amendments: `PREDECLARED_GATES.md`,
  `REMEDIAL_AMENDMENT_2026-08-29.md`,
  `NUMERICAL_CONTROL_AMENDMENT_2026-08-29.md`, and
  `NUMERICAL_CONTROL_AMENDMENT_N40_2026-08-29.md`;
- full claim history: `CLAIM_LEDGER.md`;
- generated core audit: `final_summary/audit.json`;
- generated numerical audit: `final_controls/audit.json`;
- report source: `report/long_chain_ft_report.tex`.
