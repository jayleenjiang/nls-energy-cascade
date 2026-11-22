# Section 4.3: eigenfunctions of the generator (Figs. 17–19, Table `tab:spectrum`)

Programs:
- `sec4_3_spectrum.cpp`: n = 3 chain; accumulates the gEDMD matrices A = E[φ φᵀ], B = E[φ (Lφ)ᵀ] for the dictionary d = 5, M = 3 (2744 functions) while it simulates, without storing samples.
- `sec4_3_return.cpp`: E_x[g(X_t)] from a fixed initial state x, averaged over 4096 × 16 trajectories in 32 groups.

## Spectrum (Table `tab:spectrum`, Figs. 17–18)
- Runs: `jobs_as_run/run_step2.sbatch` (NESS (2,8) and equilibrium (5,5), dt = 1e-3 and 2.5e-4). Original notes: `jobs_as_run/README_step2_original.md`.
- `analysis/analyze_step2.py` (with `gedmd_n3.py`, `ged_io.py`) solves B c = λ A c for the nested dictionaries D1–D6 and validates every mode on two halves of the data. Results: `results/spectrum/*_analysis.json` (Table), `*.json` (metadata and dictionary), logs.
- `analysis/eigfun.py` evaluates and normalises the eigenfunctions.

## Return probabilities (Fig. 19)
- Runs: `jobs_as_run/run_ex.sh` (γ = 0.1, 0.01, 0.001), `local_runs_g005_g001.sh` (γ = 0.05), `runs_g002.sh` (γ = 0.02). They call the program under its original name `local_test_v2`, with positional arguments; the equivalent one-line commands are in `../PAPER_RUNS.md`.
- `analysis/ham_orbit.py`: return times of the deterministic orbit (dotted lines in Fig. 19). `analysis/an_stat.py`: analysis of the stationary runs.

## Figures
- `analysis/reduce_figure_data.py MODES SAMPLES EX OUT` writes `results/figure_data/`: the eigenvalues and normalised coefficients of ψ0–ψ3 (`eigenfunction_coefficients.npz`), the arrays of Figs. 17–18, and the group means and standard errors of every observable of the ten return runs (`return_probabilities.npz`, run parameters in `return_probabilities_runs.json`).
- `figures/make_figures.py [OUT]` draws `fig_psi_MH.pdf`, `fig_psi1_theta.pdf`, `fig_loc_backward_gscan.pdf`.

## Not included (large)
The gEDMD accumulator files (`*.bin`, about 480 MB per case), the mode files `*_modes.npz` (27 MB per case), the n = 3 samples of Sec. 4.1 used to normalise the eigenfunctions, and the raw return-run JSON files (2.3 MB each).
