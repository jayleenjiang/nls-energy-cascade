# batch_07 (2026-10-29): §4.3 worst-case start, n = 3, 4, 5, 6

Mentor's plan: start with tiny endpoints, large total action M in the middle, and the worst phases θ = 4π/3. That angle is
the unstable zero of I_2(1 + 2cos θ_1) and of I_{n-1}(1 + 2cos θ_n). Measure (1) the decay of M(t) and (2) the last t with
M(t) ≥ L.

## Choices (made without the mentor; see make_tasks.py)
- Initial actions: I_1 = I_n = 0.01, and the middle modes share M0 equally.
- Initial angles, in the report convention:
  - endpoint bonds measured outward: θ_1 = 2(φ_1 − φ_2), θ_n = 2(φ_n − φ_{n−1});
  - interior bonds measured left to right: θ_j = 2(φ_j − φ_{j−1}).
- Three phase types:
  - worst: all angles 4π/3;
  - anti: endpoint angles 4π/3, interior angles π (n ≥ 4);
  - random: random phases.
- All three phase types use the same seed. The noise of trajectory s is then identical across them, so they are compared on the same noise.
- (T1, Tn) = (2, 8), γ = 0.1, 128 trajectories, sampled every 0.1.

| M0 | dt | run length |
|---|---|---|
| 50 | 1e-4 | 300 |
| 100 | 1e-4 | 800 |
| 200 | 5e-5 | 1500 |

- Thresholds L = 30, 50 and M0/2 are applied in the analysis. The largest M seen in the NESS is 25.2 (n=3), 24.7 (n=4) and 20.1 (n=5), so L = 30 is never reached by stationary fluctuations.
- hires/: M0 = 100, first time unit, sampled every 1e-3. Used to see how fast the endpoint angles leave 4π/3.
- ness/: n = 6 NESS run, which gives <M>_NESS for n = 6.

## Code
- `rk4a_chain_v2.cpp` is `rk4a_chain.cpp` of batch_06 with one added option, `--phases`.
  - Without the option the output is byte-identical; this is checked in `build.sbatch`.
  - `rk4a_chain_v1_reference.cpp` is the unchanged batch_06 file (sha256 eaa74c03…).
- `analysis/worst_case_analyze.py` runs afterwards on the Mac, not on Unity.

Estimated cost: 42 tasks, about 38 core-hours. The largest task takes about 25 min on 8 cores. Dumps take about 1.4 GB.
