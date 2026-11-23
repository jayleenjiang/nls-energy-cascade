# Numerical study of the boundary-driven cascade chain of cubic NLS

| paper | folder | program | figures / tables |
|---|---|---|---|
| Sec. 3 (intro), Fig. 1 | `sec1_model_geometry/` | `geometry/main.cpp` (exact resonant families) | `cascade_embedding` |
| 3.1 Action profiles | `sec3_1_action_profile/` | `sec3_1_action_profile.cpp` | Figs. 2–4 |
| 3.2 Local equilibrium | `sec3_2_local_equilibrium/` | `sec3_2_local_equilibrium.cpp` | Fig. 5 |
| 3.3 Action-current conductivity | `sec3_3_action_current/` | `sec3_3_action_current.cpp` | Figs. 6–7, Table `tab:current-exponents` |
| 3.4 Fluctuation relation | `sec3_4_heat_fluctuation/` | `sec3_4_heat.cpp` | Fig. 8, Table `tab:heat-fr-summary` |
| 4.1 NESS density (n = 3, 4) | `sec4_1_stationary_density/` | `sec4_1_density.cpp` + TensorFlow solver | Figs. 9–12, Table `tab:ness-density` |
| 4.2 Mechanism of stabilization | `sec4_2_stabilization/` | `sec4_2_stabilization.cpp` | Figs. 13–16 |
| 4.3 Eigenfunctions of the generator | `sec4_3_eigenfunctions/` | `sec4_3_spectrum.cpp`, `sec4_3_return.cpp` | Figs. 17–19, Table `tab:spectrum` |

Each section folder has its own `README.md`, which explains where every number and figure comes from and lists what is not in the repository.

## Layout of a section folder

```
secX_Y_name/
  <program>.cpp      the simulation program; its default options are the paper settings
  jobs_as_run/       case tables, SLURM scripts and shell scripts exactly as they were run
  analysis/          scripts that turn the raw output into the numbers of the paper
  results/           small outputs: run summaries, tables, and figure_data/ (arrays reduced from the raw data)
  figures/           make_figures.py: redraws the paper figures from results/
```

The raw trajectory data (several GB per section) are not included. Every raw file can be regenerated with the commands in `PAPER_RUNS.md` and the scripts in `jobs_as_run/`.

## Requirements

- C++17 compiler with OpenMP (tested with g++ 13.3 and clang++ 18.1). Eigen 3.4 headers are included in `eigen/` (MPL2, see `eigen/COPYRIGHT_Eigen_3.4.0.txt`).
- Boost (header only) for `sec1_model_geometry/geometry/main.cpp`.
- Python 3 with numpy, scipy, matplotlib; a LaTeX installation with `lmodern` (figure text uses `text.usetex`).
- TensorFlow 2.21 only for the density solver of Sec. 4.1 and for `sec4_1_stationary_density/analysis/reduce_figure_data.py`.

## Build

From the repository root, for example:

```
g++ -O3 -march=native -fopenmp -std=c++17 -I eigen sec3_1_action_profile/sec3_1_action_profile.cpp -o sec3_1_action_profile
```

and in the same way for the other seven programs. The options and their defaults are listed in the USAGE block at the top of each program. With `--dry-run 1` a program prints the parameters of the run (in the option format of the program originally used) and stops.


