# Discrete path-ratio FT audit

Overall: **BLOCKED**

| Gate | Status |
|---|---:|
| numerical_integrity | PASS |
| forward_ift | BLOCKED |
| reverse_ift | PASS |
| crooks_histogram | PASS |
| overall | BLOCKED |

## Integral fluctuation relation

| ensemble | log mean exp(-Sigma) | 95% CI | ESS | max fraction |
|---|---:|---:|---:|---:|
| forward | -0.00310973 | [-0.00576673, -0.000350492] | 82093.9 | 0.000296 |
| reverse | -0.000868817 | [-0.00358574, 0.00196346] | 82316.0 | 0.000241 |

## Bidirectional Crooks histogram

Accepted bins: 51; slope = 1.02343 +/- 0.0128; intercept = -0.000910846 +/- 0.00459; identity weighted RMSE = 0.0312083.
