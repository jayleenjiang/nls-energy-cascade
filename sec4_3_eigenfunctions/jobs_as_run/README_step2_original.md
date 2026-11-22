# Section 4.4, Step 2 — streaming gEDMD production (n = 3)

Generator Galerkin projection (gEDMD) of the n = 3 chain, gamma = 0.1, for NESS (T1,T3) = (2,8) and
equilibrium (5,5), each at dt = 1e-3 and 2.5e-4.

| file | role |
|---|---|
| `NLS3_gedmd_stream.cpp` | simulation + streaming accumulation. Dynamics copied unchanged from `spectral_gap_controlled_2026-09-07/source/NLS_stationary_autocorr_controlled.cpp` (Cartesian, double, projection-free). Accumulates A = E[phi phi^T], B = E[phi (L phi)^T] for the largest dictionary (d=5, M=3, K=2744) per stream group, plus lagged products of 6 observables (autocorrelation). Samples are never stored. |
| `gedmd_n3.py` | reference Python implementation (dictionary, generator, solver) |
| `ged_io.py` | reader for the C++ output |
| `analyze_step2.py` | eigenvalues for nested dictionaries D1..D6, out-of-sample validation of every mode (half A vs half B), modal weights, predicted vs measured autocorrelation |
| `run_step2.sbatch` | SLURM array, 4 tasks (one per case) |

Per case: 256 streams x 1e4 time units, snapshot 0.1 -> 2.56e7 samples (8x Step 1 in total time).

Verified before delivery (cloud container):
- generator: L* rho_eq = 0 symbolically; Hamiltonian part matches paper (1.10); analytic L phi vs finite differences 7e-9;
- C++ `--from` on the Step-1 sample file reproduces the Python A, B to 2e-15 (relative), including the nested sub-block extraction;
- simulation mode runs with K = 2744 (smoke test).

Outputs per case: `CASE.bin` (~480 MB, raw accumulators), `CASE.json` (metadata + dictionary),
`CASE.log`, `CASE_analysis.json`, `CASE_modes.npz`, `CASE_analysis.log`.
Only the `*_analysis.json`, `*_modes.npz`, `*.json`, `*.log` files are needed back on the Mac.
