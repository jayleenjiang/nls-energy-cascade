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
| forward | -0.00281136 | [-0.0096322, 0.00464091] | 47282.5 | 0.002596 |
| reverse | 0.00157372 | [-0.00743402, 0.0125271] | 30874.6 | 0.003071 |

## Bidirectional Crooks histogram

Accepted bins: 25; slope = 0.978539 +/- 0.00554; intercept = 0.000481233 +/- 0.00337; identity weighted RMSE = 0.0195217.
