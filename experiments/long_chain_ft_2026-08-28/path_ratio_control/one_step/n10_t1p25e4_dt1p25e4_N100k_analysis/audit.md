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
| forward | 0.000114943 | [-0.00068693, 0.000991139] | 98236.2 | 3.808e-05 |
| reverse | 0.000471102 | [-0.000303445, 0.00126488] | 98222.7 | 4.745e-05 |

## Bidirectional Crooks histogram

Accepted bins: 56; slope = 0.942641 +/- 0.0388; intercept = -0.000228409 +/- 0.0045; identity weighted RMSE = 0.029585.
