# Burn-in and finite-time flux-tail diagnostic, 2026-07-10

This folder contains a new diagnostic experiment for two questions:

1. How should the burn-in time be chosen from no-burn-in relaxation data?
2. How does `log P[flux(tau) > A]` depend on the averaging time `tau` and threshold `A`?

The simulator is kept separate from the existing production current-scaling workflow:

- source: `flux/NLS_flux_relaxation_tau.cpp`
- binary used here: `flux/flux_relax_tau`
- analysis script: `flux/analyze_burnin_ld.py`

The measured quantity is the conserved action current

```text
J_j = 4 I_{j-1} I_j sin(2(phi_j - phi_{j-1})).
```

It is called "flux" in the diagnostic plots for consistency with the meeting
question, but it is the action current used in the manuscript.

## 1. No-burn-in transient pilot

Parameters:

- `(T1,Tn,gamma) = (10,2,0.1)`
- `dt = 5e-4`
- `n = 10,20,40,80`
- `1024` trajectories per `n`
- no burn-in
- checkpoints `T = 50,100,...,500`

For each checkpoint, the simulator records:

- terminal mean action profile `E[I_j(T)]`
- cumulative profile average over `[0,T]`
- cumulative current average over `[0,T]`
- last-interval current average over `(T-50,T]`

Main outputs:

- `analysis/burnin_terminal_profile.png`
- `analysis/burnin_cumulative_profile.png`
- `analysis/burnin_flux_timeseries.png`

Preliminary interpretation:

- For `n=10`, the current is essentially stable by `T=100`.
- For `n=20`, the current is close to stable by `T=100--200`.
- For `n=40`, the last-interval current is near the production value by about `T=200--300`; the cumulative no-burn-in current continues to drift because early transient data are included.
- For `n=80`, the profile at `T=50` is visibly underdeveloped, while by about `T=150--200` the terminal profile shape is much closer to the later profiles. The current is small and noisy, so this is only a qualitative relaxation diagnostic.
- For the finite-time tail pilot below, we used a conservative common burn-in `B=500` for `n <= 40`.

## 2. Finite-time flux-tail pilot

Parameters:

- `(T1,Tn,gamma) = (10,2,0.1)`
- `dt = 5e-4`
- `burnin = 500`
- `measure = 200`
- base block `tau_block = 20`
- `tau = 20,40,...,200`
- `n = 10,20,30,40`
- `1024` trajectories per `n`

For each trajectory, the simulator stores ten consecutive block-current values
of length `20`.  The analysis uses prefix averages to construct one
`J_tau` sample per trajectory for each `tau`.

Main outputs:

- `analysis/tau_tail_survival_n10.png`
- `analysis/tau_tail_survival_n20.png`
- `analysis/tau_tail_survival_n30.png`
- `analysis/tau_tail_survival_n40.png`
- `analysis/tau_logprob_vs_tau_n10.png`
- `analysis/tau_logprob_vs_tau_n20.png`
- `analysis/tau_logprob_vs_tau_n30.png`
- `analysis/tau_logprob_vs_tau_n40.png`
- `analysis/tau_survival_surface.csv`
- `analysis/tau_survival_slope_fits.csv`

Preliminary interpretation:

- The histogram center narrows as `tau` grows, as expected for a time-averaged current.
- For thresholds below the mean current, `P[J_tau > A]` can increase with `tau` because the distribution concentrates around a positive mean.
- For thresholds above the mean current, `P[J_tau > A]` decreases with `tau`.
- At fixed `tau`, the upper-middle survival tail is locally close to linear in `A` on a log plot, but the fitted window is sample-limited for `1024` trajectories.
- The rate proxy `-tau^{-1} log P[J_tau > A]` does not cleanly collapse across `tau` in this pilot. Therefore this should still be described as a finite-time tail diagnostic, not as a demonstrated large-deviation principle.

At `tau=200`, the pilot gives:

| n | mean J | std | descriptive lambda | R^2 | fitted A window |
|---:|---:|---:|---:|---:|---:|
| 10 | 0.39595 | 0.06071 | 31.4 | 0.9918 | [0.45,0.55] |
| 20 | 0.11753 | 0.03026 | 62.0 | 0.9909 | [0.15,0.19] |
| 30 | 0.05387 | 0.02037 | 91.0 | 0.9955 | [0.08,0.10] |
| 40 | 0.02980 | 0.01632 | 126.1 | 1.0000 | [0.05,0.07] |

These values are consistent in scale with the earlier canonical finite-time
distribution table, but the burn-in and seed are from this new diagnostic
pilot.

## Full 100000-trajectory run

The simulator uses `LANES=16`, so `100000` trajectories corresponds to
`6250` batches.  With the current local 8-thread timing and `burnin=500`, the
four-length full run is roughly an overnight job:

- `n=10`: about 1.5 hours
- `n=20`: about 2.3 hours
- `n=30`: about 3.3 hours
- `n=40`: about 4.5 hours

Example commands:

```bash
for n in 10 20 30 40; do
  prefix="Paper/revision/experiments/flux_validation/burnin_ld_2026-07-10/full100k_tau_n${n}"
  flux/flux_relax_tau tau 10 2 "$n" 6250 500 200 0.0005 20 20260710 8 "$prefix"
done

python3 flux/analyze_burnin_ld.py \
  --experiment-dir Paper/revision/experiments/flux_validation/burnin_ld_2026-07-10 \
  --threshold-max 0.7
```

If we decide to use the older manuscript burn-ins instead of the transient-based
`B=500`, the full run becomes substantially longer, especially for `n=30,40`.
