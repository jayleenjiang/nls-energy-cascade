# Section 3.4: entropy production and the fluctuation relation (Fig. 8, Table `tab:heat-fr-summary`)

Program: `sec3_4_heat.cpp` (fixed frame, energy E = H/2; integrates the bath heat and the bond action current over blocks of length 5 and writes one row per block).

## Runs
- n = 10 affinity sweep (Table, Fig. 8): `jobs_as_run/affinity_sweep_CASES.tsv` (temperature pairs, dt, seeds). One-line commands in `../PAPER_RUNS.md`.
- Equilibrium baseline at n = 10 (finest dt): `jobs_as_run/equilibrium_baseline_CASES.tsv`, `equilibrium_baseline_COMMANDS.tsv`.
- Longer chains (correlations in the text): `jobs_as_run/long_chain_production.sh`.

## Where each number comes from
| item | file |
|---|---|
| Table, (6.5,5.5): a∞/Δβ = 0.99352 [0.98358, 1.00121] | `results/affinity_sweep_n10/extrapolation_summary.csv`, row `weak_dt1p25e4` |
| Table, (8,4): 0.95988 [0.95229, 0.96548] | same file, row `moderate_dt1p25e4` |
| Table, (7,5): 0.99573 [0.98854, 1.00088] | see "Not included" below |
| per-window slopes | `results/affinity_sweep_n10/per_tau_results.csv` |
| equilibrium slope (zero within error) | `results/affinity_sweep_n10/equilibrium_weighted_baseline.csv`, `results/equilibrium_baseline_n10/` |
| correlation of heat and action current, 0.81 (n = 10, τ = 20) and 0.10 (n = 40) | `results/long_chain_production/coupling_summary.csv` |

Scripts: `analysis/analyze_affinity_sweep.py` (with `analyze_pooled_extrapolation.py`), `analysis/analyze_equilibrium_baseline.py`, `analysis/analyze_long_chain_production.py`. They are the scripts as run; they read the raw block files from folders given by relative paths, which are not part of the repository.

## Figure
`figures/make_figures.py [OUT]` draws `heat_ft_grid.pdf` from `figures/fig8_data_extracted_from_pdf.json`. The per-bin histogram ratios behind Fig. 8 were computed from the raw driven block files, which are no longer stored. The points, error bars and fit lines were therefore read back from the vector graphics of the original figure (`figures/heat_ft_paper_grid_original.pdf`) and redrawn in the paper style. The numbers in the Table do not depend on this.

## Not included
- The raw driven block files of the n = 10 affinity sweep (only the equilibrium raw blocks were kept).
- The (7,5) row of the Table comes from the timestep study `analysis/extrapolation_summary.csv` (sha256 91f4413c…), which is not among the stored files; the value is recorded in the report `FT_RESULTS_REPORT.tex`, section on the primary verdict at the finest timestep.
