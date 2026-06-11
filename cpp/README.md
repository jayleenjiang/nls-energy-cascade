# C++ code — what each file is for

All long-chain integrators share the same physics (Case-0 closed chain, Gibbs/FD
bath, the two boundary bug-fixes: `φ`-drift sign `+=` and the `√2` boundary
noise). The variants below differ only in step size, sampling, conditioning, or
instrumentation.

## `simulation/` — long-chain SDE integrators

| File | What it does | Feeds |
|------|--------------|-------|
| `main.cpp` | Baseline Euler–Maruyama integrator (pre-fix). | — (historical) |
| `main_fixed.cpp` | Corrected integrator (boundary `φ`-sign, `√2` noise). Physics base for everything below. | all later runs |
| `lte_dump.cpp` | Early approach: dumps decorrelated samples to CSV for offline fitting. | superseded by the accumulators |
| `lte_histogram.cpp` | First on-the-fly 3D histogram of `(I_j, I_{j+1}, θ_j)` (Eigen). | §3.1③ LTE |
| **`lte_histogram_simd.cpp`** | **Canonical LTE accumulator.** SIMD over trajectories, `fast_sin`, 1D marginals in C++, stationary-mean-field init, `n²`-scaled burn-in. Writes `.hist` + `.marg`. | §3.1③ LTE — `.hist/.marg` → `python/analysis/hist_analysis_bootstrap.py` + `matlab/fit.m` |
| `lte_simd_conditioned.cpp` | Same, but **fixed-S conditioning**: only accumulates when \|Σⱼ Iⱼ − S\| < δ (S = 0.5n). | LTE local-Gibbs-under-conditioning test |
| `lte_dttest.cpp` | Convergence run at `dt = 5e-4`. | large-`n` Euler-artifact diagnosis |
| `lte_simd_dt25.cpp` | Convergence run at `dt = 2.5e-4`. | same — with `dttest`, shows overflow rate ∝ `dt²` (n=100 fat tails) → §3.1 n=100 caveat / box-convention fit |
| `lte_drift.cpp` | Instrumented dt-stability probe: tracks running max drift / relative step, prints to stderr. | same Euler-artifact diagnosis |
| `lte_adaptive.cpp` | Adaptive shared step (`dt` chosen from current drift; `dt_max=1e-3`, `dt_min=1e-7`). | large-`n` stability alternative |
| `flux_V1.cpp`, `flux_V2.cpp` | Energy-flux measurement `J_j = 4 I_{j-1}I_j sinθ_j` for the conductivity scaling. | §3.2 thermal conductivity (note: `flux_V2` had a `φ`-sign bug — re-confirm) |

## `fp5d/` — three-mode (5D) system `(I₁,I₂,I₃,θ₁,θ₃)`

| File | What it does | Feeds |
|------|--------------|-------|
| `NLS5D_SIMD_NCKDE5D.cpp`, `_v2.cpp` | 5D SDE + noisy-channel KDE density estimate. | §4.1 5D density |
| `NLS5D_SIMD_marginal.cpp`, `_2.cpp` | Marginal density extraction. | §4.1 |
| `NLS5D_SIMD_nckde.cpp`, `_nckde_slice.cpp` | NC-KDE density / conditional slices. | §4.1 |
| `NLS5D_SIMD_check2.cpp`, `_check3.cpp` | Validation against the exact Gibbs marginal. | §4.1 equilibrium validation |
| `NLS5D_slice3.cpp` | Conditional-slice extraction. | §4.1 |

## `backward/` — generator / eigenfunction

| File | What it does | Feeds |
|------|--------------|-------|
| `NLS_backward.cpp` | Backward Monte-Carlo / generator solver. | §4.3 eigenfunction, spectral gap λ_R |
| `Kuramoto_SIMD_backward.cpp` | Kuramoto analogue — method validation. | §4.3 method check |
