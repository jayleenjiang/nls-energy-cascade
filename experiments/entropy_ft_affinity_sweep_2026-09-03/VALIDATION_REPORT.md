# Validation report: n=10 heat-flux affinity sweep

Final reanalysis date: 2026-09-04

## Change made

The production data were not regenerated or modified. The histogram spacing
(standard deviation divided by 20), minimum raw count (10 on both sides),
minimum number of symmetric pairs (3), and contiguous straddle-zero rule are
unchanged.

The original gate required all full-sample-resolved times to remain resolved
simultaneously in a stream-bootstrap replicate. It rejected the known
equilibrium control and reported 0/1000 for every case. Analysis amendment
002 adds:

1. an independent stream bootstrap for each time; and
2. an intercept bootstrap over reliably resolved times only.

The complete old analysis and report are preserved in
pre_gate_respec_2026-09-04/. The amendment specification is
AMENDMENT_002_BOOTSTRAP_GATE.md.

## Old and amended gates

For driven cases, reliable times satisfy both n_negative at least 500 and
full-sample matched-bin R-squared at least 0.98. At equilibrium the literal
R-squared rule gives an empty set because the exact reference is a flat zero
response. The equilibrium known-answer intercept therefore uses all
full-sample-resolved times with n_negative at least 500 and explicitly omits
the inapplicable nonzero-signal R-squared condition.

| Delta beta | old times | old reported | old raw valid | amended times | amended valid |
|---:|---|---:|---:|---|---:|
| 0 | 20,40,80,160,320,640 | 0/1000 | 29/1000 | 20,40,80,160,320,640 (equilibrium exception) | 1000/1000 |
| 0.027972 | 20,40,80,160,320,640 | 0/1000 | 19/1000 | 20,40,80,160,320 | 1000/1000 |
| 0.057143 | 20,40,80,160,320 | 0/1000 | 103/1000 | 20,40,80,160 | 1000/1000 |
| 0.125000 | 20,40,80 | 0/1000 | 38/1000 | 20,40,80 | 1000/1000 |
| 0.222222 | 20,40 | 0/1000 | 0/1000 | 20,40 | no intercept |
| 0.400000 | 20 | 0/1000 | 0/1000 | 20 | no intercept |

The original output displayed zero whenever raw acceptance was below the
predeclared 800 threshold. The raw-valid column above exposes the actual
sub-threshold intersection count without changing the old verdict.

## Per-time intervals

All 23 full-sample-resolved rows now have a per-time interval. Every row has
1000/1000 accepted stream-bootstrap replicates except Delta beta 0.057143 at
t=320, which has 976/1000. Exact values are in
analysis/window_summary.csv and in the PDF report.

Finite-time intervals containing their reference:

- equilibrium: t = 20,40,80,160,320,640;
- Delta beta 0.027972: t = 80,320,640;
- Delta beta 0.057143: t = 80,160,320;
- Delta beta 0.125000: t = 80;
- Delta beta 0.222222: none;
- Delta beta 0.400000: none.

The t=640 row at Delta beta 0.027972 is excluded from the amended intercept
because its R-squared is 0.9704. The t=320 row at Delta beta 0.057143 is
excluded because n_negative is 150 and R-squared is 0.6020.

## Reliable-time intercept results

| Delta beta | a_infinity | 95% CI | a_infinity / Delta beta | ratio 95% CI | accepted | verdict |
|---:|---:|---|---:|---|---:|---|
| 0 | 0.00007803 | [-0.00011784, 0.00025898] | N/A | N/A | 1000/1000 | CONSISTENT_WITH_EQUILIBRIUM |
| 0.027972 | 0.0276114 | [0.0269972, 0.0280475] | 0.98711 | [0.96515, 1.00270] | 1000/1000 | CONSISTENT_WITH_FT |
| 0.057143 | 0.0566825 | [0.0550349, 0.0578691] | 0.99194 | [0.96311, 1.01271] | 1000/1000 | CONSISTENT_WITH_FT |
| 0.125000 | 0.1277489 | [0.1145073, 0.1376372] | 1.02199 | [0.91606, 1.10110] | 1000/1000 | CONSISTENT_WITH_FT |
| 0.222222 | N/A | N/A | N/A | N/A | 0/1000 | UNRESOLVED |
| 0.400000 | N/A | N/A | N/A | N/A | 0/1000 | UNRESOLVED |

The equilibrium CI contains zero. The first three driven ratio CIs contain
one. The Delta beta 0.222222 and 0.400000 cases do not have three reliable
times and remain unresolved.

The full-sample regressions of slope on 1/t have R-squared values 0.017 at
equilibrium and 0.135, 0.855, and 0.985 for the three driven intercepts.
These are reported diagnostics, not extra post-hoc gates. In particular, the
0.027972 intercept is statistically consistent with the FT reference under the
requested gate, but the low extrapolation R-squared shows that a linear 1/t
description is noisy and should not be oversold.

## Integrity and provenance

- Each of the five new raw files contains exactly 1,000,064 rows, 128 streams,
  and 7,813 contiguous block IDs per stream.
- The accepted (10,2) production was reused unchanged.
- Production source commit:
  1905cf4e606a4a7f4dd8930caa64bd4cc861e9d4.
- Production source SHA-256:
  98e7f8f5f915c8ce02bd8aa10722025c09fd739184b981961692869c9356c0d3.
- Driven protocol commit:
  e170205c5d9a3d9b2bfca3714a4447fa965a1e40.
- Equilibrium amendment commit:
  cf5282d3ab77ce93c05578ff15fb6598026081c7.
- Analysis seed: 2026090499; bootstrap replicates: 1000.

Raw input SHA-256 values:

| case | rows | SHA-256 |
|---|---:|---|
| (6,6) | 1,000,064 | 7fd6827f12b42844de1142613af124702ca5535ba63f3c54ab71006b256aef4d |
| (6.5,5.5) | 1,000,064 | 8b44be37ea8f566fb3862b6cd796dcfcf433be0888752babb876754c84d4c421 |
| (7,5) | 1,000,064 | 382496fcd194154b0f6d486f2d4da9de80d2d0791f94cb0f3e695d8e666a29cc |
| (8,4) | 1,000,064 | 0404492717ce335c1e2eff7970d5d514435fabec9b91604fa92b7d4846cdf685 |
| (9,3) | 1,000,064 | 342f3783b40db98263cd685c8ab4e531c9d83842db6762b0bd0a06061642958d |
| (10,2) reused | 1,000,064 | a23806e82f5514a9c3375d10a6644946b6b57d0efe8452b3bbf397b6230f9929 |

## Claim boundary

Under the amended, auditable gate, the equilibrium control is consistent with
zero and the three weakest driven affinities are numerically consistent with
the heat-flux FT after reliable-time extrapolation. This is finite-sample,
finite-chain numerical evidence, not a proof. The two stronger affinities
remain unresolved because directly sampled negative-tail support does not
persist for at least three reliable times. Early-time deviations and the
noisy weakest-drive extrapolation remain visible in the report.
