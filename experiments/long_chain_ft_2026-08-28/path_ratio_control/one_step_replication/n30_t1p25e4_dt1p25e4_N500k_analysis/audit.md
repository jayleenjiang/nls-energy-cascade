# Discrete path-ratio FT audit

Overall: **PASS**

| Gate | Status |
|---|---:|
| numerical_integrity | PASS |
| forward_ift | PASS |
| reverse_ift | PASS |
| crooks_histogram | PASS |
| overall | PASS |

## Integral fluctuation relation

| ensemble | log mean exp(-Sigma) | 95% CI | ESS | max fraction |
|---|---:|---:|---:|---:|
| forward | 0.000203272 | [-0.00119645, 0.00155126] | 413615.1 | 5.393e-05 |
| reverse | 0.000446827 | [-0.000767381, 0.00172518] | 408938.4 | 0.0001277 |

## Bidirectional Crooks histogram

Accepted bins: 69; slope = 0.994052 +/- 0.00538; intercept = -0.000182324 +/- 0.00204; identity weighted RMSE = 0.0162166.
