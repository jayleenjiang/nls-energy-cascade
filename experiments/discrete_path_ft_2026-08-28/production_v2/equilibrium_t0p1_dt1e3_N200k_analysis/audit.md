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
| forward | 0.00504293 | [-0.0014766, 0.0140398] | 48438.2 | 0.003279 |
| reverse | -0.000418568 | [-0.00305684, 0.00233037] | 141057.2 | 0.0005593 |

## Bidirectional Crooks histogram

Accepted bins: 40; slope = 0.981888 +/- 0.0073; intercept = 0.000123879 +/- 0.00326; identity weighted RMSE = 0.0215388.
