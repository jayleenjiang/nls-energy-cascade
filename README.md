# Numerical study of the boundary-driven cascade chain of cubic NLS

Code and the small data behind the figures and tables of the paper. There is one folder per section of the paper, and each folder holds the one program used for that section, together with the job files as they were run, the analysis scripts, the reduced results and a script that redraws the figures.

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

## Paper runs

`PAPER_RUNS.md` gives a one-line command for each of the 552 runs of the paper, for example

```
./sec3_1_action_profile --n 50 --T1 10 --Tn 5 --seed 202610110002
```

Only the options that differ from the paper settings are written; `--threads` and `--out` are left to the user (the values used are in `jobs_as_run/`).

## Figures

```
python3 secX_Y_name/figures/make_figures.py OUTPUT_FOLDER
```

redraws the figures of that section from `results/`. Figs. 2–19 were checked against the PDF files of the paper: rendered at 100 dpi, the redrawn figures and the paper figures differ in 0 pixels. Fig. 1 has the same content and a slightly different size. Fig. 8 is drawn from data read back from the original figure file, not from computed data; see `sec3_4_heat_fluctuation/README.md`.

## Verification of the programs

The eight programs are rewritten, commented versions of the programs used for the paper. `verify/` contains the original programs and two checks:

- `verify/verify_against_original.sh` runs 20 short cases with both the new and the original program and compares every output file. All 20 cases give identical output with g++ 13.3 and with clang++ 18.1 (`verify/results_gcc13.md`, `verify/results_clang18.md`).
- `verify/check_paper_runs.py` checks, for every line of `PAPER_RUNS.md`, that the program expands it to exactly the parameters of the original run (552 runs, 0 mismatches) and writes `PAPER_RUNS.md`.

The model, the integrators and the bath are described in Section 2 of the paper; the header of each program summarizes what it computes, the scheme and the output files.
