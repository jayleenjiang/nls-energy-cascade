# Section 4.2: mechanism of stabilization (Figs. 13–16)

Program: `sec4_2_stabilization.cpp` (RK4 in action–angle variables for the Hamiltonian part, Cartesian Euler–Maruyama baths; modes can be held at zero with `--hold`; initial states with `--M0`, `--start worst|anti|random`, `--init`, `--phases`).

## Runs
| runs | script | used for |
|---|---|---|
| local runs: held modes, Λ sweep, drain from a full middle mode, NESS | `jobs_as_run/local_runs/run_rk4a_all.sh` | Fig. 13, Λ ≈ 9, Λ > 0 for all 25 pairs, drain times for n = 3 |
| Unity batch 06: n = 4 slow decrease (E0 scan), NESS references | `jobs_as_run/unity_batch06/` | n = 4 full-middle-mode results |
| Unity batch 07: three initial phases, n = 3–6, M0 = 50, 100, 200 | `jobs_as_run/unity_batch07_worstcase/` | Fig. 14, last time with M ≥ 30 |
| small γ, end actions, t ≤ 0.3 | `jobs_as_run/local_smallgamma/run_smallg.sh` | Fig. 15 |
| small γ, mean of M | `jobs_as_run/local_Mrise/run_Mrise.sh` | Fig. 16 |
| Unity batch 08: small γ, run until M has decayed | `jobs_as_run/unity_batch08_smallgamma/` | not used in the current text |

The scripts call the program under its original names (`rk4a_chain`, `rk4a_chain_v2`, `v2`); `sec4_2_stabilization` takes the same options. One-line commands of the Unity and Mrise runs: `../PAPER_RUNS.md`.

## Analysis
- `analysis/report13/`: the analysis of the local and batch 06 runs (`tables_n3_n4.py`, `rk4a_angles_sweep.py`, `drain_analyze_rk4a.py`, `trap_scan_analyze_rk4a.py`, …); results in `results/report13/`. `README_report13_original.md` maps each number to its script and result file. Some of these results (NESS tails, n = 5, events) are not used in the paper.
- `analysis/fokker_planck/`: stationary Fokker–Planck equation of the two end angles at fixed (I_1, I_3), for the canonical bath and the bath of HLNS; table `results/fokker_planck/frozen_map_table.md` (the ratios −0.15 and close to −√3/2 in the text).
- `analysis/worstcase/worst_case_analyze.py`: last passage times and opposite-phase fractions for batch 07 → `results/worstcase/worst_case.json`, `worst_case_hires.json`; `last_passage_mean.py`, `last_above_M0.py` → `results/worstcase/last_above_M0.json`.
- `results/Mrise/Mrise.json`: largest rise of E M(t) − M(0) and its time (also printed by the figure script).
- `analysis/validation/schemes4_corot.cpp`: the comparison of four schemes from which the RK4 action–angle scheme was taken.

## Figures
- `analysis/reduce_figure_data.py ANGLES WORST SMALLG MRISE OUT` reduces the dumps to `results/figure_data/`.
- `figures/make_figures.py [OUT]` draws `fig43_boundary_angles.pdf`, `fig43_worstcase_M.pdf`, `fig43_smallgamma_endpoints.pdf`, `fig43_smallgamma_M.pdf`.

## Not included (large)
All trajectory dumps (batch 06/07 outputs, local runs; the small-γ dumps are 0.2 GB and the Mrise dumps 2 GB). They are reproducible from the scripts and seeds.
