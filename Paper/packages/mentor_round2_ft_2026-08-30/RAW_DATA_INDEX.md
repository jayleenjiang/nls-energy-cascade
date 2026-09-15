# Large raw-data index

The delivery ZIP is intentionally compact.  The following accepted raw files
remain in the local repository and are identified by immutable SHA-256 hashes.
All derived CSVs, audits, and figures needed to inspect the conclusions are
inside this package.

## Direct `t=20` blocks

| chain | bytes | SHA-256 | absolute path |
|---:|---:|:---|:---|
| 10 | 147,784,241 | `a23806e82f5514a9c3375d10a6644946b6b57d0efe8452b3bbf397b6230f9929` | `/Users/jayleenjiang/Documents/NLS/experiments/entropy_ft_2026-08-26/production/n10_blocks.csv` |
| 20 | 149,305,409 | `68617ba13ad37f65f0f2cd9b1b61b14c41396b1a9bf315b99e5a611a45777d1f` | `/Users/jayleenjiang/Documents/NLS/experiments/entropy_ft_2026-08-26/production/n20_blocks.csv` |
| 30 | 150,281,336 | `5660003a206a86c930fcf00b65f87dee78b1187aa097bf344d64451484782dbc` | `/Users/jayleenjiang/Documents/NLS/experiments/entropy_ft_2026-08-26/production/n30_blocks.csv` |
| 40 | 150,787,130 | `70c511177fa141257f711586676a4444a539f17ee74936ec7a7926a5a23d7e9b` | `/Users/jayleenjiang/Documents/NLS/experiments/entropy_ft_2026-08-26/production/n40_blocks.csv` |

Each file contains 1,000,064 finite, ordered, non-overlapping blocks.  The
package includes the independent raw and analysis audit reports.

## Two-site NESS total-entropy blocks

| condition | bytes | SHA-256 | absolute path |
|:---|---:|:---|:---|
| driven `(10,2)` | 301,923,355 | `65aa2c9a4da2c23a016d88ca9f05ead4f3029c630a70c24ed83e5e297d166612` | `/Users/jayleenjiang/Documents/NLS/experiments/entropy_ft_scgf_2026-08-27/total_entropy_n2_short/production/driven_blocks.csv` |
| equilibrium `(6,6)` | 302,064,684 | `5d8de682435bf1de52ea9cd17573634382c1be89db77b9cd647fe805e1c9e4b2` | `/Users/jayleenjiang/Documents/NLS/experiments/entropy_ft_scgf_2026-08-27/total_entropy_n2_short/production/equilibrium_blocks.csv` |

## Exact discrete path-ratio trajectories

| direction | bytes | SHA-256 | absolute path |
|:---|---:|:---|:---|
| forward | 226,988,237 | `254b689856c029c686c416b2beeacee602c25ca31fb4ed042b1d7f86f4e96aef` | `/Users/jayleenjiang/Documents/NLS/experiments/discrete_path_ft_2026-08-28/production_v2/driven_t0p1_dt1e3_N1m_forward.csv` |
| reverse | 226,985,525 | `a7307635f65f9f92b62ac5d26fe4731674116977f3a52f122d2f6de760637ee0` | `/Users/jayleenjiang/Documents/NLS/experiments/discrete_path_ft_2026-08-28/production_v2/driven_t0p1_dt1e3_N1m_reverse.csv` |

## Long-chain Phase-II raw output

The 56 accepted Phase-II summaries, 56 timeseries, and 56 logs total only a
few megabytes and are included directly under
`results/long_chain_phase2/raw/`.
