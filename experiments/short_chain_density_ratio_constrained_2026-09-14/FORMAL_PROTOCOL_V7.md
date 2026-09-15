# Formal V7 fine-timestep validation freeze

Frozen on 2026-09-14 before any V6 physical statistic was read.  The V6
integrity-only audit passed.  R9 selected the common calibration shrink factor
0.625 by its predeclared cross-development safety score.

## Exact driven bundles

| base seed | reported seed | weights SHA-256 | config SHA-256 |
|---:|---:|---|---|
| 8201 | 43201 | `7c76f5cbd00b01db27a1d9f831c4b809c4db95deff47079a43a69a9f937ea322` | `aa232ad08f02453abece3f371023225bea3180819dfc319bd7de0e3a33c0e207` |
| 8202 | 43202 | `881b96ef4aa69cc94ac6f8b7ca239cbbefdcc8a4d582f0efa0802587aa221738` | `ee9f7b80d26c23d7989f37711e914551b2e0970a21220f7cd7816ab0ac8e8907` |
| 8203 | 43203 | `8fd45218b837b7026b6128d7e2578f06415ce1ca4986b0bf6faf392dc447084f` | `aa844eebf4ee9b524006ec895fe9d8f376934831c9724a02af66a138f789d887` |

R9 selection SHA-256:
`c5920463e706bc3d342aa852724e90a3968124fec63ae6e93529a0a41c1e7d3b`.
V6 split SHA-256:
`c76d5f1317905a54aefef7ccb82a73a32fbd9653825465db2cb3bbebcb00ec24`.
V6 audit SHA-256:
`35173c93c46d3de9cb892fb77fc53681b223f445183d59b7b480a8dcf8325e52`.

All original formal gates remain unchanged, including the moment cutoff
3.748286557303567 and the FP median/p90 cutoffs 0.10/0.50.  Validation creates
the marginal edges.  V6 test remains sealed unless all three driven models and
the exact-zero equilibrium known-answer control pass validation completely.
