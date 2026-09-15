# Instrumentation audit

The original and instrumented sources were compiled with the same optimization, Eigen, and OpenMP flags. Both were run with:

```text
T1=6 Tn=6 n=3 batches=1 burnin=0.01 measure=0.01 dt=0.0005 seed=2026090601 threads=1
```

The complete current-sample CSVs are byte-identical, with common SHA-256:

```text
a636ec29ea69cf46bdee8a7338a1d184004082fb77c6e58b8250d5d4c5e9cc41
```

The maximum absolute difference between the three original and instrumented mean-action profile values is `2.220446049250313e-16`. Thus the added profile and near-floor accumulators did not alter the sampled dynamics in this deterministic smoke comparison.

Original-audit binary SHA-256:

```text
0056a4cfc9805dad04ac7b2870d2a39c50660364f8aa15707358f31cebfd93db
```

