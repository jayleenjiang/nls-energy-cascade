# Section 3.2: local equilibrium (Fig. 5)

Program: `sec3_2_local_equilibrium.cpp` (the integrator of Sec. 3.1 plus histograms of (I_a, I_{a+1}, θ_{a+1}) at three bonds, 80³ bins, one file per batch and bond).

## Runs
- `jobs_as_run/cases_lte.csv`, `run_lte.sbatch`: NESS (10,5) and equilibrium T = 7.5, n = 25, 50, 100, two replicates; histogram bonds at 24%, 48%, 72% of the chain.
- `jobs_as_run/analyze_lte.sbatch` runs the analysis on the cluster.

## Analysis and figure
- `analysis/analyze_lte_v2.py`: the fixed analysis (fit of log q against log p of the NESS and rescaled equilibrium histograms). Its outputs are in `results/`: `lte_summary_v2.csv` (slopes quoted in the text), `aggregated_hist_v2.npz` (histograms summed over batches, input of the figure), `analysis_manifest_v2.json`, and the 12 run summaries in `results/summaries/`.
- The per-batch histogram files (`*.histbin`) are not included.
- `figures/make_figures.py [OUT]` draws `lte_residual_mesh_hot_side.pdf` and checks the fitted slopes against `lte_summary_v2.csv`.
