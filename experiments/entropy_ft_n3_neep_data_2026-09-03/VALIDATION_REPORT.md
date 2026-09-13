# n=3 consecutive-pair data for a future NEEP test

## Outcome

Data generation and file-integrity validation passed for both requested
conditions.  No NEEP model was trained or evaluated.

Each condition contains 128 independent streams and 39,063 consecutive
`delta_t=0.1` transitions per stream: 5,000,064 transition rows and 5,000,192
actual state snapshots after counting the unrepeated first state of each
stream.  Every transition spans 200 unchanged Cartesian integration steps.

The two NEEP-ready archives total 1,592,478,048 compressed bytes.  Including
the exact raw simulator archives, the complete retained dataset is
2,754,680,210 bytes (2.5655 GiB).

## Integrity results

All four Zstandard archives pass decompression integrity checks.  Each archive
contains 5,000,065 lines including its header.  Both cases pass all of the
following full-data checks:

- 128 stream ids and every interval id are present in frozen order;
- all values are finite, every action is positive, and all angles lie in
  `[-pi,pi)`;
- adjacent interval endpoints agree exactly at saved precision;
- the stored first-law residual equals `Q_L+Q_R-Delta E` exactly;
- the stored medium entropy agrees with `-Q_L/T_L-Q_R/T_R` to at most
  `1.77636e-15` in the driven data and `8.88178e-16` at equilibrium;
- every appended sine/cosine value agrees with its raw angle exactly under the
  same double-precision recomputation; the largest unit-circle identity error
  is `2.22045e-16`.

There were zero midpoint-solver failures in either production run.

## First-law residual per 0.1 interval

| case | mean | standard deviation | RMS | RMS per unit time | maximum absolute |
|---|---:|---:|---:|---:|---:|
| driven `(10,2)` | `-7.59971e-8` | `3.45284e-5` | `3.45285e-5` | `3.45285e-4` | `3.26722e-3` |
| equilibrium `(6,6)` | `-2.76947e-9` | `2.79177e-5` | `2.79177e-5` | `2.79177e-4` | `1.90710e-3` |

Driven residual quantiles `(0.1%,1%,50%,99%,99.9%)` are
`(-2.96998e-4,-9.72514e-5,-1.05853e-10,9.62679e-5,2.97245e-4)`.
Equilibrium values are
`(-2.34559e-4,-8.36251e-5,-5.04999e-12,8.34989e-5,2.36698e-4)`.

## Consecutive-state correlation

The angle statistic is periodic: it uses the centered complex representation
`exp(i theta1)`, not a linear correlation of wrapped angle values.

| case | observable | lag-1 rho | 1/e time | integrated time | `delta_t/(1/e time)` |
|---|---|---:|---:|---:|---:|
| driven | `I2` | 0.684385 | 0.214157 | 0.544955 | 0.466948 |
| driven | periodic `theta1` | 0.406770 | 0.113788 | 0.144758 | 0.878830 |
| equilibrium | `I2` | 0.688453 | 0.198284 | 0.531220 | 0.504327 |
| equilibrium | periodic `theta1` | 0.466315 | 0.132835 | 0.151443 | 0.752813 |

Across driven streams, the 5th--95th percentile range of the `I2` 1/e time is
`[0.197560,0.254385]`, and its integrated-time range is
`[0.481441,0.640631]`.  For periodic `theta1` the corresponding ranges are
`[0.111504,0.115807]` and `[0.133236,0.159884]`.

At equilibrium the `I2` ranges are `[0.192636,0.223398]` and
`[0.485762,0.615645]`; periodic `theta1` gives `[0.130631,0.135029]` and
`[0.138048,0.164579]`.

Therefore `delta_t=0.1` retains nonzero consecutive-state correlation, but it
is not asymptotically small compared with the fastest measured coordinate.
In particular it is 0.879 of the driven periodic-angle 1/e time and 0.753 of
the equilibrium value.  This is a temporal-resolution warning for later NEEP
training, not a failed data-integrity check.  No finer-interval data were
silently generated because the requested interval was frozen at 0.1.

## Data locations

The directly trainable archives are:

- `raw/driven_neep_transitions.csv.zst`;
- `raw/equilibrium_neep_transitions.csv.zst`.

Each row keeps all exact simulator columns and appends cosine and sine columns
for `theta1` and `theta3` at both endpoints.  The exact unchanged simulator
outputs are retained as `raw/driven_blocks.csv.zst` and
`raw/equilibrium_blocks.csv.zst`.  Hashes and byte sizes are in
`DATA_MANIFEST.csv`; detailed numerical output is in `analysis/`.
