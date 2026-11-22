# Section 4.1: NESS density of the three- and four-mode chains (Figs. 9–12, Table `tab:ness-density`)

Program: `sec4_1_density.cpp` (the Cartesian integrator; it writes state dumps, or a streamed kernel density estimate at given query states).

## Steps
1. **Labels** (`jobs_as_run/prep_labels.sbatch`): Gibbs references at T = 5 by Metropolis (`labels/gibbs_reference_n.py`), independent trajectories (10⁶ samples each), query states (`labels/make_queries*.py`) and KDE labels from 10⁹ samples for n = 3, 4, equilibrium (5,5) and NESS (2,8). `labels/eval_eq*.py` checks the equilibrium labels against the exact density. Label run summaries: `results/label_summaries/`.
2. **Long validation trajectories** (`jobs_as_run/val_long.sbatch`): 256 streams, burn-in 500, 10⁴ time units, two independent sets A and B per case. Summaries: `results/validation_summaries/`.
3. **Training** (`jobs_as_run/train_v5.sbatch`): `training/train_kde_nn_v5.py` (n = 4), `training/train_kde_nn_n3_v5.py` (n = 3); network and reduced Hamiltonian in `training/n{3,4}_{model,common}.py`. The trained networks are in `models/` (equilibrium and two seeds of the NESS for each n).
4. **Evaluation** (`jobs_as_run/eval_v5.sbatch`, `training/eval_check_long.py`, `training/diag_sin_n4.py`): `results/eval/`.
5. **Table** `tab:ness-density`: `analysis/an_sec4_v5.py MODELS DENSITY VAL OUT`, run from `training/` (it imports the model files from the working folder); output `results/an_sec4_v5.json`, log `results/an_sec4_v5_run.log`.

SLURM logs of steps 2–4: `jobs_as_run/slurm_logs/`. Original notes: `jobs_as_run/README_v5_original.md`.

## Figures
- `analysis/reduce_figure_data.py GIBBS VAL OUT` evaluates the networks on the Gibbs samples and on angle grids and writes the small arrays in `results/figure_data/` (needs TensorFlow and the large inputs below).
- `figures/make_figures.py [OUT]` draws `fig_marginals_n3.pdf`, `fig_marginals_n4.pdf`, `fig_theta_slices_n3.pdf`, `fig_theta_slices_n4.pdf` from these arrays (no TensorFlow needed).

## Not included (large)
Gibbs reference samples `gibbs_n{3,4}_T5.bin`, the trajectory dumps, the query and label files, and the validation dumps `n{3,4}_neq_val{A,B}.bin`. The validation dumps used for the figures were regenerated with the commands of `val_long.sbatch`; their summaries agree with `results/validation_summaries/` except for the run time.
