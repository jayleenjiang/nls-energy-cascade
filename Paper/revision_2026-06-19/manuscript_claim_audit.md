# Manuscript claim audit

Generated: `2026-06-20T19:29:46.208529+00:00`

Scope note: Core numerical/data claims only; author declarations and external plagiarism checks remain outside code-verifiable scope.

Summary: **18 / 18** claims verified; **0** failed.

| ID | Section | Verdict | Evidence |
|---|---|---:|---|
| `flux_scaling_main` | abstract / thermal conductivity | VERIFIED | `Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/flux_primary_scaling.json`<br>`Paper/revision_2026-06-19/experiments/flux_validation/validation_report.md` |
| `conductivity_scaling` | abstract / thermal conductivity | VERIFIED | `Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/flux_primary_scaling.json` |
| `flux_table_values` | thermal conductivity | VERIFIED | `Paper/revision_2026-06-19/experiments/flux_validation/validation_report.md`<br>`Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/flux_primary_scaling.json` |
| `flux_diagnostics` | thermal conductivity | VERIFIED | `Paper/revision_2026-06-19/experiments/flux_validation/validation_report.md` |
| `larger_n_current_robustness` | thermal conductivity | VERIFIED | `Paper/revision_2026-06-19/experiments/flux_validation/larger_n_pilot_2026-06-20/n10_50_b64_scaling_scaling.json`<br>`Paper/revision_2026-06-19/experiments/flux_validation/larger_n_pilot_2026-06-20/n50_b64_summary.csv`<br>`Paper/revision_2026-06-19/experiments/flux_validation/larger_n_pilot_2026-06-20/n50_b16_dt2p5e-4_summary.csv`<br>`Paper/revision_2026-06-19/experiments/flux_validation/larger_n_pilot_2026-06-20/n50_b16_burn10000_summary.csv`<br>`Paper/revision_2026-06-19/experiments/flux_validation/larger_n_pilot_2026-06-20/README.md` |
| `n60_current_robustness` | thermal conductivity | VERIFIED | `Paper/revision_2026-06-19/experiments/flux_validation/larger_n60_pilot_2026-06-20/n10_60_b64_scaling_scaling.json`<br>`Paper/revision_2026-06-19/experiments/flux_validation/larger_n60_pilot_2026-06-20/n60_b64_summary.csv`<br>`Paper/revision_2026-06-19/experiments/flux_validation/larger_n60_pilot_2026-06-20/README.md` |
| `bath_parameter_current_robustness` | thermal conductivity | VERIFIED | `Paper/revision_2026-06-19/experiments/flux_validation/parameter_robustness_2026-06-20/moderate_contrast_T8_T4_prod/b64_scaling_scaling.json`<br>`Paper/revision_2026-06-19/experiments/flux_validation/parameter_robustness_2026-06-20/moderate_contrast_T8_T4_prod/n10_b64_summary.csv`<br>`Paper/revision_2026-06-19/experiments/flux_validation/parameter_robustness_2026-06-20/moderate_contrast_T8_T4_prod/n20_b64_summary.csv`<br>`Paper/revision_2026-06-19/experiments/flux_validation/parameter_robustness_2026-06-20/moderate_contrast_T8_T4_prod/n30_b64_summary.csv`<br>`Paper/revision_2026-06-19/experiments/flux_validation/parameter_robustness_2026-06-20/moderate_contrast_T8_T4_prod/n40_b64_summary.csv`<br>`Paper/revision_2026-06-19/experiments/flux_validation/parameter_robustness_2026-06-20/README.md`<br>`Paper/revision_2026-06-19/experiments/flux_validation/parameter_robustness_2026-06-20/production_summary.csv` |
| `flux_scaling_fit_sensitivity` | thermal conductivity | VERIFIED | `Paper/revision_2026-06-19/experiments/flux_validation/larger_n60_pilot_2026-06-20/flux_scaling_sensitivity_n10_60.json`<br>`Paper/revision_2026-06-19/experiments/flux_validation/larger_n_pilot_2026-06-20/n50_b64_summary.csv`<br>`Paper/revision_2026-06-19/experiments/flux_validation/larger_n60_pilot_2026-06-20/n60_b64_summary.csv`<br>`Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/flux_primary_scaling.json` |
| `finite_window_current_statistics` | finite-time current fluctuations | VERIFIED | `Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/current_windows_window_statistics.csv` |
| `long_chain_action_profiles` | nonequilibrium steady state | VERIFIED | `Paper/revision_2026-06-19/manuscript_figure_metrics.json` |
| `lte_table_values` | local thermodynamic equilibrium | VERIFIED | `Paper/revision_2026-06-19/source_trace_metrics.json` |
| `lte_control_values` | local thermodynamic equilibrium | VERIFIED | `Paper/revision_2026-06-19/source_trace_metrics.json` |
| `lte_equilibrium_convention` | local thermodynamic equilibrium | VERIFIED | `Paper/revision_2026-06-19/draft.tex` |
| `short_chain_equilibrium_validation` | short-chain Fokker--Planck | VERIFIED | `Paper/revision_2026-06-19/short_chain_nn_rerun_metrics.json`<br>`Paper/revision_2026-06-19/source_trace_metrics.json` |
| `short_chain_symmetry_scope` | stabilization | VERIFIED | `Paper/revision_2026-06-19/short_chain_nn_rerun_metrics.json` |
| `eigen_relaxation_diagnostic` | eigenfunction | VERIFIED | `Paper/revision_2026-06-19/eigen_fit_sensitivity.json`<br>`Paper/revision_2026-06-19/short_chain_nn_rerun_metrics.json` |
| `reproducibility_summary_table` | numerical reproducibility summary | VERIFIED | `Paper/revision_2026-06-19/draft.tex`<br>`Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/flux_primary_scaling.json`<br>`Paper/revision_2026-06-19/experiments/flux_validation/larger_n60_pilot_2026-06-20/flux_scaling_sensitivity_n10_60.json`<br>`Paper/revision_2026-06-19/experiments/flux_validation/parameter_robustness_2026-06-20/moderate_contrast_T8_T4_prod/b64_scaling_scaling.json`<br>`Paper/revision_2026-06-19/experiments/flux_validation/validation_report.md`<br>`Paper/revision_2026-06-19/manuscript_figure_metrics.json`<br>`Paper/revision_2026-06-19/source_trace_metrics.json`<br>`Paper/revision_2026-06-19/short_chain_nn_rerun_metrics.json`<br>`Paper/revision_2026-06-19/eigen_fit_sensitivity.json` |
| `data_availability_artifacts` | data and code availability | VERIFIED | `Paper/revision_2026-06-19/draft.tex`<br>`Paper/revision_2026-06-19/availability_path_audit.json`<br>`Paper/revision_2026-06-19/availability_path_audit.md` |

## Failed text checks

None.

## Notes

- This report is a local claim/data audit. It does not replace author confirmation of funding, contributions, or competing interests.
- It also does not replace a professional plagiarism/self-plagiarism service check before formal journal submission.
