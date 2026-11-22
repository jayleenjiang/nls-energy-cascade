# Section 3.3: action-current conductivity (Figs. 6–7, Table `tab:current-exponents`)

Program: `sec3_3_action_current.cpp` (the integrator of Sec. 3.1; it also writes the time-averaged bond currents of every trajectory, `<label>_current_streams.txt`, used for the standard errors).

## Runs
- `jobs_as_run/make_tasks.py` writes `tasks_A.csv` (linear response: T1, Tn = 5 ∓ ΔT/2, n = 15–45, and chain length n = 20–150 at (4,6)) and `tasks_B.csv` (chain length n = 160–250); `run_A.sbatch`, `run_B.sbatch`; `build_check.sbatch` built the program and compared it with the previous version before the runs.

## Analysis, table and figures
- `results/runs/{resp,len}/`: summary and per-trajectory currents of every run.
- `analysis/fig33_reproduce.py results/runs OUT` computes the centre-bond current with standard errors over trajectories (runaway trajectories excluded and counted), the fits, and writes `current_response.csv`, `current_length.csv`, `fits_length.txt`. Its output is identical to `results/tables/`.
- Table `tab:current-exponents`: `results/tables/fits_length.txt`.
- `figures/make_figures.py [OUT]` draws `fig33_current_response.pdf` and `fig33_current_length.pdf` from `results/tables/`.
