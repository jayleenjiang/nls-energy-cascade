# report_13: all code and raw data (2026-10-27)

Every number in `report/report_13.tex` can be recomputed from this package:
```
bash code/run_all_analysis.sh
```
- Needs python3 with numpy, scipy and matplotlib.
- Runs in about 20–30 min.
- Writes to a new folder `recomputed/`; nothing in the package is overwritten.
- Tested: `tables_n3_n4.py`, the sweep range, n5, trap, drain and trap_scan reproduce `results/` exactly.

## Integrator (all simulations)
- `code/simulation/rk4a_chain.cpp`, sha256 eaa74c03…. It is the same file in `batch_06_unity/`.
  - Hamiltonian part: RK4 in action–angle variables.
  - Bath: Cartesian Euler–Maruyama at sites 1 and n.
  - Functions copied verbatim from `code/reference/schemes4_corot.cpp` (the four-scheme test, 10/14), scheme RK4_A.
- Build: `g++ -O3 -fopenmp -std=c++17 rk4a_chain.cpp -o rk4a_chain`.
- Local runs (in the cloud workspace, 2 cores): `code/simulation/run_rk4a_all.sh`. Its lock-angle runs are not used in report_13 and not included here.
- Unity runs: `code/simulation/batch_06_unity/` (tasks.tsv, build/run sbatch, prompt).
  - Unity jobs 66192915 (build) and 66192916 (array, 35 tasks); every task finished with `nonfinite=0`.
  - The Unity slurm logs and per-task `.log` / `.command.txt` files are not in this package. They are still in `Unity_batch06_download_2026-10-27_2Yrwa9/`.

## Report item → code → data
| report_13 section | numbers | code (`code/…`) | raw data (`data/…`) | stored result (`results/…`) |
|---|---|---|---|---|
| Numerics | integrator | simulation/rk4a_chain.cpp | local_rk4a/val (Gibbs and dt checks) | rk4a_compare.json → validation |
| n=3 case 1 | E, Λ table | analysis/tables_n3_n4.py | local_rk4a/case1/n3_hold2_* | rk4a_compare.json → case1 |
| n=3 case 1 | sweep "Λ = 2.7–12.2" | analysis/rk4a_angles_sweep.py | local_rk4a/sweep (25 runs) | sweep_rk4a.json (Λ 2.732–12.119) |
| n=3 case 1 | angle peak at 4π/3 | analysis/rk4a_angles_sweep.py | local_rk4a/angles | fig43_boundary_angles_rk4a.png |
| n=3 case 2 | F(I1, I3) table (canonical (2,8), (5,5), HLNS) | fokker_planck/frozen_phase_fp.py, frozen_phase_map.py, fig_frozen_map.py; checks: frozen_phase_checks.py, frozen_phase_sde.py | none: Fokker–Planck solve, no simulation | fokker_planck/frozen_map_part_*.json, frozen_map_table.md, checks .txt |
| n=3 case 3 | drain / angle / return table, ⟨M⟩_NESS = 4.70 | analysis/tables_n3_n4.py | local_rk4a/drain (n3_E10…200), local_rk4a/ness | rk4a_compare.json → case3 |
| n=4 cases 1, 2 | E, Λ₂, Λ₃ tables | analysis/tables_n3_n4.py | local_rk4a/case1/n4_hold2_*, n4_hold23_* | rk4a_compare.json → case1 |
| n=4 case 3 | formed fraction, τ, leak, prediction, M at formation and break-up | analysis/trap_scan_analyze_rk4a.py (+ theory_trap.py for theory_trap.json) | unity_batch06/e0scan, gscan, tscan; ness/n4_T2_8 (M_ref) | trap_scan.json |
| n=5 | leak by pair, lifetimes, cluster fractions | analysis/n5_analyze.py | unity_batch06/n5 | n5.json |
| NESS | mean M, 99.9 % quantile, max, Gibbs quantiles | analysis/tail_analyze.py | unity_batch06/ness | tail.json |
| NESS | κ, dips, rise slope, angles in dips | analysis/events_analyze.py | unity_batch06/events, ness | events.json |
| NESS | pair fractions with M > 15 | analysis/pair_gibbs.py | unity_batch06/ness/n4_T2_8 | pair_gibbs.json |
| NESS | "long-lived pair needs E0 ≳ 50" | analysis/trap_analyze_rk4a.py | unity_batch06/trap | trap.json |
| (figure only, not in report) | drain figure n=3, 4 | analysis/drain_analyze_rk4a.py | local_rk4a/drain + unity_batch06/drain | drain.json, fig_drain.png |

`analysis/rk4a_compare.py` is the script that first produced the local-run numbers. It also compares with the older aa_chain and Cartesian runs. The aa_chain runs are not in this package. The Cartesian reference JSONs it reads are in `results/cartesian_reference_for_rk4a_compare/`. `tables_n3_n4.py` reproduces the same numbers from the rk4a runs alone.

## Data format
Each run has the following files.
- `NAME_summary.txt`:
  - parameters;
  - `meanI_j`, `rate_j` (time average of the Hamiltonian d log I_j/dt), `Sleft_j`, `Sright_j`, each given as mean and SE over the streams;
  - `nonfinite_steps`.
- `NAME_dump.f32` (if sampled):
  - float32, layout (streams/16, snapshots, 16, 2n−1);
  - columns I_1..I_n, then θ_2..θ_n with θ_j = 2(φ_j − φ_{j−1}) in (−π, π];
  - the notes' θ_1 is −θ_2.
- The sampling interval is in `code/simulation/run_rk4a_all.sh` (local runs) or `data/unity_batch06/tasks.tsv` (Unity runs).

## Not included
- Raw data of the older Cartesian and aa_chain runs: not used in report_13.
- The lock-angle runs: that table was removed from the report.
- The −0.48 bound: it comes from Yao Li's slides (UofA.pdf, p. 34), not from code.
