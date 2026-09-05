# Numerical study of an energy cascade model derived from a dispersive equation

Code for a computational study of the stochastic NLS energy-cascade model: a
closed chain of `n` modes in action–angle coordinates `(I_j, φ_j)`, driven out
of equilibrium by stochastic heat baths at the two ends and simulated with
Euler–Maruyama. The project accompanies Hani–Li–Nahmod–Staffilani
(arXiv:2505.16018), extending their three-mode analysis numerically to long
chains and to a 5D Fokker–Planck / eigenfunction study of the three-mode system.

## What's here

| Path | Contents |
|------|----------|
| `cpp/simulation/` | Long-chain SDE integrators and on-the-fly accumulators |
| `cpp/fp5d/` | Three-mode (5D) SDE simulator + density estimation (NC-KDE) |
| `cpp/backward/` | Backward / adjoint solvers for the eigenfunction analysis |
| `notebooks/` | Neural-network Fokker–Planck and eigenfunction solvers (TF/Keras) |
| `python/analysis/` | Histogram processing, LTE fits, phase-locking peak finding |
| `python/plotting/` | Figure generation |
| `python/data_gen/` | Training-data generation for the backward/NN solvers |
| `scripts/` | CMake build file and run scripts |
| `paper/` | Working draft and the LTE results report (LaTeX) |

### C++ — `cpp/simulation/`
- `main.cpp` — baseline Euler–Maruyama integrator for the `n`-mode closed chain.
- `main_fixed.cpp` — same, with two corrections: boundary `φ`-drift sign and the
  `√2` factor in the boundary noise.
- `lte_dump.cpp` — dumps decorrelated samples for offline LTE fitting.
- `lte_histogram.cpp` — accumulates 3D histograms `(I_j, I_{j+1}, θ_j)` on the fly
  (Eigen-based).
- `lte_histogram_simd.cpp` — SIMD/structure-of-arrays rewrite: vectorizes across
  trajectories, polynomial `fast_sin`/`fast_log`, stores 1D marginals in C++, and
  reaches ~1e9 samples by default.
- `flux_V1.cpp`, `flux_V2.cpp` — energy-flux measurement for the thermal-conductivity
  scaling. (Note: `flux_V2` carried a `φ`-sign bug; re-confirm scaling results
  against the corrected integrator.)

### C++ — `cpp/fp5d/`
Three-mode system in coordinates `(I_1, I_2, I_3, θ_1, θ_3)`, simulated with SIMD,
with several density-estimation and validation variants:
`NLS5D_SIMD_NCKDE5D[_v2].cpp` (main NC-KDE density estimator),
`NLS5D_SIMD_marginal[_2].cpp`, `NLS5D_SIMD_nckde[_slice].cpp`,
`NLS5D_SIMD_check2/3.cpp` (validation against the Gibbs marginal),
`NLS5D_slice3.cpp` (conditional slices).

### C++ — `cpp/backward/`
- `NLS_backward.cpp` — backward Monte-Carlo / generator solver for the spectral
  (eigenfunction) analysis.
- `Kuramoto_SIMD_backward.cpp` — Kuramoto analogue used to validate the method.

### Notebooks
- `FKE_5d_NLS.ipynb` — neural-network solver for the 5D Fokker–Planck operator,
  validated against the exact Gibbs measure.
- `FKE_4d.ipynb` — reduced 4D system.
- `FKE_eigen.ipynb` — neural-network eigenfunction / spectral-gap estimation.

## Building the C++ code

Linux (GCC):
```
g++ -O3 -march=native -ffast-math -fopenmp -I/usr/include/eigen3 \
    cpp/simulation/lte_histogram_simd.cpp -o lte_histogram_simd
```

macOS (Apple Silicon, clang + libomp):
```
clang++ -O3 -mcpu=native -ffast-math -std=c++17 \
    -I/opt/homebrew/include/eigen3 \
    -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include \
    -L/opt/homebrew/opt/libomp/lib -lomp \
    cpp/simulation/lte_histogram.cpp -o lte_histogram
```

Dependencies: a C++17 compiler, OpenMP, and (for the Eigen-based files) Eigen 3.
The SIMD rewrite `lte_histogram_simd.cpp` is Eigen-free.

Example run (NESS, n=25, three bulk sites):
```
./lte_histogram_simd 10 2 25 histo_N25 6 12 18
```

## Python / notebooks

Python analysis needs `numpy`, `scipy`, `matplotlib`, `pandas`; the notebooks
additionally need `tensorflow`/`keras`. The analysis scripts read the `.hist`
and `.marg` files produced by the C++ accumulators.

## Notes

- Hamiltonian convention: the code uses `H_code = 2 · H_paper`, so the invariant
  measure is `exp(−H_code / 2T) = exp(−H_paper / T)`.
- Data files, compiled binaries, and figures are intentionally **not** tracked
  (see `.gitignore`); regenerate them from the source here.
- The `paper/` folder is included for convenience; remove it if you want the
  repository to be code-only.
