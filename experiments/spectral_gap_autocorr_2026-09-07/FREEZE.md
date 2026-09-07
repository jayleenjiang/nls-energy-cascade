# Pre-production freeze record

The protocol, source, analysis, and exact commands were frozen before either
production case was launched.

- Historical autocorrelation target: none.  The files
  `KDE/NLS_backward_Y_train.txt` and `Paper/revision/eigen_fit_sensitivity.json`
  contain conditional expectations, not stationary trajectories.
- Historical curve SHA-256:
  `e2b0474d566a6837e660da30e5fbcb7f36441a7600895d33259e8589f8f58141`
- Historical source SHA-256 (`cpp/backward/NLS_backward.cpp`):
  `789926a6abc8b839dfad818a82e4c4fd4d247bc87bb715edf5212d6ed2d1dd24`
- New sampler source SHA-256:
  `a1d558c267a4eca00c3b8f445e1236c1d097be50e72845f18ed4a4e12d7c9df1`
- New production binary SHA-256:
  `447ba1141ba5f126bccd6f39661c939eee4683ef96fcd12a1595420e7e532316`
- Frozen analysis SHA-256:
  `35d93b7d721972bdd0f847db76ac9656a71bf69d17021edc83a5e9f0eaf4b064`

Smoke test: four independent driven streams, burn `0.1`, measurement `0.2`,
snapshot interval `0.01`; all four files had the exact expected shape, finite
values, and zero non-finite flags.  A synthetic exponential with rate `-1.5`
was recovered by the local-rate differencing function to machine precision.
