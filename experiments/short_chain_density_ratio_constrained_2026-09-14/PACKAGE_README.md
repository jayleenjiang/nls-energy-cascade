# Five-dimensional NESS density package

This package contains the final normalized density model, frozen test results,
publication figures, the standalone report, and the replacement paper section.

## Main entry points

- `report/five_dimensional_ness_density_report.pdf`: detailed seven-page report;
- `paper/section4_1_ness_density.tex`: manuscript-ready Section 4.1;
- `scripts/final_ness_density.py`: callable absolute-density evaluator;
- `final_analysis/normalizer_summary.csv`: ratio and absolute normalizers;
- `final_analysis/final_gate_summary.csv`: all final fine/coarse/equilibrium gates;
- `final_timestep_comparison/verdict.json`: final timestep gate;
- `FINAL_VERDICT.md` and `VALIDATION_REPORT.md`: claim and audit boundary.

## Density evaluation

Input CSV columns are `I1,I2,I3,theta1,theta3`.  From the experiment root:

```bash
/Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11 \
  scripts/final_ness_density.py input.csv output.csv --timestep fine
```

The output adds `log_rho_ss` and `rho_ss`.  The fine estimator is the primary
result; the coarse estimator is the timestep replication.

## Raw data

The raw trajectory arrays are deliberately not duplicated in the ZIP.  They
remain at:

- `/Users/jayleenjiang/NLS_section4_density_cache/short_chain_density_v6_holdout_2026-09-14` (final fine holdout);
- `/Users/jayleenjiang/NLS_section4_density_cache/short_chain_density_v4_holdout_2026-09-14` (coarse holdout);
- `/Users/jayleenjiang/NLS_section4_density_cache/short_chain_section4_2026-09-14` (fitting data).

`V4_HOLDOUT_SPLITS.csv`, `V6_HOLDOUT_SPLITS.csv`, and the corresponding audit
JSON files record paths, roles, and per-stream SHA-256 hashes.

## Claim boundary

The complete reduced five-dimensional NESS density is numerically resolved on
sampled stationary support with held-out and timestep replication.  Accuracy
in extreme unsampled tails is not established, and this is not an analytic
global solution.
