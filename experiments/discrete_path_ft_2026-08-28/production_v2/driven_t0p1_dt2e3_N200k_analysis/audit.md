# Discrete path-ratio FT audit

Overall: **BLOCKED**

| Gate | Status |
|---|---:|
| numerical_integrity | PASS |
| forward_ift | BLOCKED |
| reverse_ift | BLOCKED |
| crooks_histogram | PASS |
| overall | BLOCKED |

## Integral fluctuation relation

| ensemble | log mean exp(-Sigma) | 95% CI | ESS | max fraction |
|---|---:|---:|---:|---:|
| forward | -0.0100544 | [-0.0167994, -0.00167313] | 55785.9 | 0.001718 |
| reverse | -0.0082907 | [-0.0151235, -0.000235049] | 54534.9 | 0.002189 |

## Bidirectional Crooks histogram

Accepted bins: 25; slope = 0.997141 +/- 0.00559; intercept = -0.000620035 +/- 0.00337; identity weighted RMSE = 0.0201817.
