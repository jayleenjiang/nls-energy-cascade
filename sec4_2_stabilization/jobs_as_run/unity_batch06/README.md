# batch_06 (2026-10-26): remaining §4.3 runs with rk4a_chain (schemes4 RK4_A + Cartesian bath)

These are the parts of the §4.3 report that still came from the Cartesian code:
- n = 4 case 3 (first drain, trap, E0 / γ / T scans);
- n = 5;
- NESS (ness, events).

They are redone with the same integrator as `04_03_stabilization/report_rk4a_2026-10-26/`:
- `rk4a_chain.cpp`, identical to that folder's `code/rk4a_chain.cpp` (sha256 eaa74c03…);
- same cases, parameters and seeds as `run_hi_events.sh`, `run_trap.sh` and `run_ext.sh`.

Not repeated, because they were already done with rk4a in report_rk4a: the n = 3 drain, the lock runs, and the n = 3 trap runs (the return time comes from the drain runs).

| file | use |
|---|---|
| `tasks.tsv` | 35 cases: folder, name, arguments |
| `build.sbatch` | compiles (without -march=native) and runs a short Gibbs check |
| `run_array.sbatch` | array 0-34, 8 cores each; output goes to `out/<folder>/` and is never overwritten |
| `UNITY_PROMPT.md` | prompt for Unity |

- Cost: about 14 core-hours in total. The largest task, gscan n4_E200_g0.02, takes about 20–30 min on 8 cores.
- Dumps: about 1.6 GB in total.
- Analysis is done afterwards on the Mac, not on Unity.
  - The existing scripts read this dump layout unchanged.
  - Two changes are needed. `trap_scan_analyze.py` reads `interval` from the summary, and rk4a summaries do not have that line, so the value is taken from `tasks.tsv`. The hard-coded Cartesian NESS means are replaced by the rk4a NESS means.
  - Both changes go in new copies of the scripts.
