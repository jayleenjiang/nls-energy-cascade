# batch_08 (2026-10-29): §4.3 worst-case start at small γ, run until M has decayed

Follow-up of batch_07 (γ = 0.1). The mentor asked for very small γ. Same start as batch_07:
- I_1 = I_n = 0.01, the middle modes share M0 = 100 equally;
- (T1, Tn) = (2, 8);
- phases worst / anti / random, as in batch_07.

## Parts
| part | γ | run length | sampled every | tasks |
|---|---|---|---|---|
| A (`tasks_A.tsv`, `run_A.sbatch`) | 1e-3 | 3·10⁴ | 1 | 11 |
| B (`tasks_B.tsv`, `run_B.sbatch`) | 1e-4 | 3·10⁵ | 10 | 11 |

- Tasks per part: n = 3, 4, 5, 6 × {worst, anti (n ≥ 4), random}.
- Common settings: 32 trajectories, dt 1e-4, seed 801.
- Output: `out/g1e-3/` and `out/g1e-4/`.

## Why these run lengths
Local test in the cloud container: n = 4, worst and random, 16 trajectories.

Last t with M ≥ 30:

| γ | worst | random |
|---|---|---|
| 1e-2 | median ≈ 320, max 523 | either ≤ 4, or 260–890 |
| 1e-3 | 34–78 for 5 trajectories; 1000–5900 for the other 11 | ≤ 300 for 10 trajectories; 5500–7500 for 3; not finished by t = 8000 for 3 |

- The slow trajectories sit in the state "two neighbouring middle modes hold ≥ 90% of M with cos θ < −0.9". At γ = 1e-3 this was 62% of worst and 38% of random trajectories at t = 2000.
- Their time grows like 1/γ. The run length is therefore 30/γ: about 5× the slowest finished n = 4 trajectory at γ = 1e-3.
- γ = 0 check (no bath, n = 4 worst start, t ≤ 1000): M is conserved to 4e-8.
- In the same γ = 0 run only 2.6% of samples are in that state. So the state is reached because of the bath, not by the Hamiltonian alone.

γ = 1e-5 and 1e-6 are not included: they would need t ~ 3·10⁶ to 3·10⁷, which is too expensive. The early times at those γ are in the local runs (`worstcase_smallgamma_2026-10-29/`).

## Cost
Estimate from batch_07 on Unity: 3.1e-6 core-seconds per step per trajectory for n = 6; less for smaller n.

| part | core-hours | wall time per task (16 cores) | dumps |
|---|---|---|---|
| A | ≈ 50–80 | ≈ 0.5 h | ≈ 0.3 GB |
| B | ≈ 500–800 | ≈ 4–6 h (limit set to 2 days) | ≈ 0.3 GB |

## Code
- `rk4a_chain_v2.cpp`: the same file as batch_07 (sha256 171c60e5…). Not modified.
- `build.sbatch` checks:
  - CHECK1: Gibbs <M> for n = 4, T = 5;
  - CHECK2: M drift with γ = 0.
- Every task line was dry-run in the cloud container with a 0.5 time-unit window: all 22 ran, nonfinite = 0.
- The analysis runs afterwards on the Mac.
