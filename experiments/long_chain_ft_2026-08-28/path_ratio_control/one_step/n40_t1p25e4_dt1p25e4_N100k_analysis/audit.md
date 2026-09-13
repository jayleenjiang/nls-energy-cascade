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
| forward | -0.000144004 | [-0.00403164, 0.00357616] | 68082.3 | 0.0006925 |
| reverse | -0.000325367 | [-0.00385917, 0.00312637] | 71814.4 | 0.0001705 |

## Bidirectional Crooks histogram

Accepted bins: 49; slope = 0.989127 +/- 0.00989; intercept = 0.000140684 +/- 0.00465; identity weighted RMSE = 0.0405927.
