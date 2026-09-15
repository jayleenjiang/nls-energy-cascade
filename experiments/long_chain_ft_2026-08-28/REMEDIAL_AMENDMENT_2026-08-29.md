# Frozen remedial amendment: n=30 population convergence

Date frozen: 2026-08-29, before generating any `N_c=2048` data for this
amendment.

The first frozen final-summary audit exposed a failed population-member gate
for the `n=30`, `(k,1-k)=(0.4,0.6)`, `t=60`, `N_c=512 -> 1024`
comparison.  The paired GC residual, support, final-time symmetry, and
`N_c=1024` time-convergence gates passed, but both individual SCGF members did
not satisfy the predeclared absolute-and-two-combined-SE population criterion.
The failed comparison remains part of the record and no threshold is changed.

The only remedial data permitted by this amendment are four new independent
runs per member at `N_c=2048`, using the unchanged settings

```text
T_L=10, T_R=2, n=30, burnin=500, horizon=60,
selection_time=2, dt=0.0005, gauge_shift=0.1,
control_scale=0.5, resample_threshold=1,
k in {0.4,0.6}, seeds 89801,...,89808.
```

The frozen final candidate is the `N_c=2048`, `t=60` pair.  It is admissible
only if all existing gates pass without modification:

1. support for both pair members;
2. paired GC gate at `t=60`;
3. time-convergence gate from `t=40` to `t=60` at `N_c=2048`;
4. individual-member population gates from `N_c=1024` to `N_c=2048`;
5. paired-residual population gate from `N_c=1024` to `N_c=2048`.

If any gate fails, `n=30` remains unresolved.  No further population size,
tilt pair, or tolerance will be selected from these data.
