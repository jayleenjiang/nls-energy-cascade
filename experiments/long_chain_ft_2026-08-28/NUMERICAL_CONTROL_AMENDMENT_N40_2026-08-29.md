# Frozen remedial amendment: n=40 selection-interval control

Date frozen: 2026-08-29, after the selection-interval-4 support failure became
observable and before generating any selection-interval-1 data below.

The pre-specified `n=40`, `N_c=512`, `selection_time=4` control is retained as
an unsupported diagnostic.  Its first `k=0.6` run already has minimum
selection-weight ESS far below the unchanged `0.1 N_c=51.2` floor.  The full
four-seed series remains in the experiment record, but a comparison without
full support cannot validate estimator invariance.

The replacement control changes only the selection interval from the supported
production baseline value 2 to 1:

```text
T_L=10, T_R=2, n=40, N_c=512, burnin=500, horizon=120,
selection_time=1, dt=0.0005, gauge_shift=0.1,
control_scale=0.5, resample_threshold=1,
k in {0.4,0.6}, seeds 90301,...,90308.
```

The interval is selected because shorter blocks improve per-selection weight
support; no observed GC residual is used.  The replacement is admissible only
if both members pass the original support gate and all member, full-window
residual, and late-half residual changes pass the unchanged absolute
`0.01`-and-two-combined-SE criteria.  Failure leaves the selection sensitivity
unresolved.  No additional interval, population, control scale, or tolerance
will be selected from these data.
