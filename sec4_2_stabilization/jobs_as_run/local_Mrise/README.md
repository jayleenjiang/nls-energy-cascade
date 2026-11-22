# Does M first rise, then fall? (2026-10-29)

Same start as the small-γ test: I_1 = I_n = 0.01, middle modes share M0 = 100, (T1, Tn) = (2, 8).
- worst: all angles 4π/3; random: random phases.
- n = 3–6; γ = 1e-3 (16384 trajectories) and γ = 1e-4 (65536 trajectories).
- dt 1e-4, t ≤ 0.1, sampled every 1e-3, seed 901.
- Integrator rk4a_chain_v2, the same file as batch_07; build it as `../code/v2`.
- Run with `run_Mrise.sh` (needs `make_tasks.py` in `../code/`), then `fig_Mrise.py`.
- Cloud container, not Unity. All 16 runs nonfinite = 0.
- Dumps are not copied (several GB); the runs are reproducible from the script.

## What is plotted
- Solid line: mean of M(t) − M(0) over all trajectories, divided by γ; band ±2 standard errors.
- Dashed line: integral of the mean bath drift, computed from the recorded I and angles:

  dM/dt (mean) = γ Σ_e [4T_e − 2(2M − I_e)I_e − 4 I_e I_neighbour cos θ_e],  e = endpoints 1 and n.

## Results (`Mrise.json`)
z = largest rise / its standard error.

| n | start | γ | largest rise / γ | z | at t | drift-integral prediction | drift turns negative at t |
|---|---|---|---|---|---|---|---|
| 3 | worst | 1e-3 | 0.74 | 22 | 0.021 | 0.79 | 0.023 |
| 3 | random | 1e-3 | 0.19 | 6 | 0.007 | 0.18 | 0.008 |
| 4 | worst | 1e-3 | 0.67 | 23 | 0.021 | 0.71 | 0.023 |
| 4 | random | 1e-3 | 0.32 | 7 | 0.014 | 0.28 | 0.012 |
| 5 | worst | 1e-3 | 0.76 | 23 | 0.025 | 0.79 | 0.025 |
| 5 | random | 1e-3 | 0.39 | 9 | 0.017 | 0.39 | 0.017 |
| 6 | worst | 1e-3 | 0.93 | 27 | 0.029 | 0.94 | 0.030 |
| 6 | random | 1e-3 | 0.53 | 10 | 0.023 | 0.50 | 0.022 |
| 3 | worst | 1e-4 | 0.96 | 17 | 0.028 | 1.05 | 0.030 |
| 3 | random | 1e-4 | 0.30 | 4 | 0.009 | 0.18 | 0.008 |
| 4 | worst | 1e-4 | 0.72 | 13 | 0.023 | 0.72 | 0.023 |
| 4 | random | 1e-4 | 0.37 | 5 | 0.014 | 0.28 | 0.012 |
| 5 | worst | 1e-4 | 0.80 | 14 | 0.026 | 0.80 | 0.026 |
| 5 | random | 1e-4 | 0.42 | 6 | 0.016 | 0.39 | 0.017 |
| 6 | worst | 1e-4 | 0.95 | 18 | 0.029 | 0.95 | 0.030 |
| 6 | random | 1e-4 | 0.52 | 7 | 0.020 | 0.50 | 0.022 |
