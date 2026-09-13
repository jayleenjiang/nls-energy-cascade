# Validation report

- New simulation: none.
- Fit window fixed before analysis: `[0.05,1.00]` for every early damped fit.
- Independent units: 64 saved trajectory correlations per case.
- Bootstrap: 2,000 trajectory-level resamples per fit; all 8,000 fits accepted.
- Damped-fit gate: frequency CI excludes zero, `Delta AIC <= -10`, and at least
  half a cycle in the fixed window (`lambda_I >= 3.30694`).
- Even-candidate plateau gate: unchanged from the original frozen experiment.
- Input SHA-256 values: recorded in `PROTOCOL.md`.
- Output tables:
  - `results/sign_changing_early_damped_fits.csv`;
  - `results/sign_changing_fit_curves.csv`;
  - `results/even_candidate_common_plateaus.json`;
  - `results/sector_summary.json`.
- Visual report: `sector_followup_report.pdf`, inspected page by page.
- Claim limitation: the saved observable set does not span exact endpoint-
  exchange sectors because `sin(theta3)` was not recorded.
