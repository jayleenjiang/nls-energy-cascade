# Section 4.4 — bath strength scan gamma = 0.1, 0.05, 0.02, 0.01 (2026-10-08)

4.4_report_gscan.tex/.pdf: your 4.4_report.tex unchanged (Sections up to "Localized test"), plus one new
section "Bath strength scan" with four tables and the same-format figure set for gamma = 0.05, 0.02, 0.01.
The figure PNGs are in the same folder, so the .tex compiles as is.

Inputs: Unity runs gamma001_run_2026-10-07, gamma005_run_2026-10-07 (and step2_run, gamma002_run).
Local supplementary runs (localized test, stationary samples for eigenfunctions): data/.

code/
- scan4.py      slow eigenvalues per dictionary, gamma-scaled out-of-sample rule 0.05|lam|+0.5 gamma -> data/scan4.json
- band.py       complex band at K = 875 (Re > -200 gamma) -> data/band875.json, band.out
- eig_gen.py    eigenfunctions at K = 1715: (M,H) R^2, exchange symmetry -> data/eig_g0.0*_D4.json
- figs_gamma.py the five same-format figures for gamma = 0.01 / 0.05
- an_stat.py    stationary indicator autocorrelations (local runs)
- local_runs_g005_g001.sh  commands of the local runs
