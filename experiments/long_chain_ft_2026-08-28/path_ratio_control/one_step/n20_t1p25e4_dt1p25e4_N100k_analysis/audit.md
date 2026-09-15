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
| forward | -0.000533409 | [-0.00248277, 0.00139616] | 92561.7 | 9.51e-05 |
| reverse | 2.16008e-05 | [-0.00181859, 0.00190128] | 92493.8 | 6.95e-05 |

## Bidirectional Crooks histogram

Accepted bins: 56; slope = 1.01341 +/- 0.0187; intercept = -0.000426662 +/- 0.00453; identity weighted RMSE = 0.0297196.
