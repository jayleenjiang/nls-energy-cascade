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
| forward | -0.00353241 | [-0.00703047, 0.000244389] | 226928.4 | 0.0006929 |
| reverse | -0.00214484 | [-0.00728945, 0.0045505] | 87731.8 | 0.00228 |

## Bidirectional Crooks histogram

Accepted bins: 33; slope = 0.993352 +/- 0.00243; intercept = 0.000258796 +/- 0.00151; identity weighted RMSE = 0.0116872.
