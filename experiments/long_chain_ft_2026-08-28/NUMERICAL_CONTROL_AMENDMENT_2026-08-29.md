# Frozen remedial amendment: n=10 numerical controls

Date frozen: 2026-08-29, before generating any `N_c=1024` numerical-control
data described below.

The first `n=10`, `N_c=512` endpoint-control audit is retained as a failed
diagnostic.  It cannot validate numerical invariance for two reasons:

1. the baseline `k=0.7` member has minimum weight ESS `45.91`, below the
   frozen `0.1 N_c=51.2` support floor, so the timestep comparison lacks full
   support even though its observed SCGF changes are small;
2. the selection-interval-4 series has still lower high-tilt ESS, and
   `t/2=30` is not on its four-unit output grid, so it cannot use the same
   late-half estimator at the frozen `t=60` horizon.

No tolerance is changed and neither failed series is discarded.  The
replacement controls use `N_c=1024`, because the existing independent
production baseline at this population has full support for both members at
`t=60`.  The exact replacement matrix is frozen as follows:

```text
baseline (existing): selection_time=2, dt=0.0005,
                     selection2_n10/N1024_t80, evaluated at t=60
timestep control:    selection_time=2, dt=0.00025, seeds 90101,...,90108
selection control:   selection_time=1, dt=0.0005, seeds 90201,...,90208
common settings:     T_L=10, T_R=2, n=10, N_c=1024, burnin=500,
                     horizon=60, gauge_shift=0.1, control_scale=0.5,
                     resample_threshold=1, k in {0.3,0.7}
```

The shorter selection interval is chosen because it keeps `t/2=30` on the
selection grid and avoids the unsupported long interval.  It is not chosen
from a GC residual.  Each replacement control is admissible only if both
members pass the original support gate and all member, full-window residual,
and late-half residual changes pass the unchanged absolute-and-two-combined-SE
criteria.  Failure leaves the corresponding numerical sensitivity unresolved;
no additional interval, timestep, population, or tolerance will be selected
from these data.
