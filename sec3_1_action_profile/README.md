# Section 3.1: action profiles and finite-size scaling (Figs. 2–4)

Program: `sec3_1_action_profile.cpp` (open or closed chain, Cartesian implicit midpoint in the co-rotating frame, Euler–Maruyama baths). Output per run: `<label>_summary.txt` with the mean action of every site, standard errors over the 256 trajectories, the bond currents and checkpoints of the profile.

## Runs
- `jobs_as_run/cases_16cpu.csv`, `cases_4cpu.csv` (143 cases) with `run_s31_16cpu.sbatch`, `run_s31_4cpu.sbatch`: groups `atlas` (BC1, BC2, BC3, BC3b × (T1,Tn) = (10,5), (4,6), (6,6) × n = 25, 50, 100, two replicates), `closed`, `dense` (n = 20 … 150 at (4,6)) and `dtcheck` (n = 100 at dt = 5e-4, compared with the dt = 1e-3 runs of `atlas`).
- `jobs_as_run/cases_rerun.csv`, `run_rerun.sbatch`: two BC3b cases in which some trajectories ran away, run again with step halving (`--safe-step 1`) and the runaway monitor (`--runaway-M 500`). The summaries in `results/` for these two cases are the reruns; the runaway trajectories (2 of 256 and 1 of 64) are left out of the means.
- One-line commands of all runs: `../PAPER_RUNS.md`.

## Analysis and figures
- `results/{atlas,closed,dense,dtcheck}/`: the 143 run summaries.
- `analysis/load31.py` reads them and merges replicates (stream-weighted means). `analysis/figs31.py` writes the tables in `analysis/data/` (scaling fits, Gibbs comparison, dt check). `analysis/gibbs_chain.py` computes the Gibbs reference profiles by Metropolis (`analysis/data/gibbs_*.csv`).
- `figures/make_figures.py [OUT]` draws `open_action_profile_atlas.pdf`, `open_action_profile_collapse.pdf`, `closed_action_profile.pdf`.
