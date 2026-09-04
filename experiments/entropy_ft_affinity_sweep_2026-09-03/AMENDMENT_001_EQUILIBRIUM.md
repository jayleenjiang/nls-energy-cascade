# Protocol amendment 001: n=10 equilibrium production

Amendment date: 2026-09-03 (America/New_York)

## Authorization and scope

After the four driven productions had started under `PROTOCOL.md`, the user
explicitly authorized adding the previously unavailable equilibrium case
`(T_L,T_R)=(6,6)` at the same requested production size.  This amendment does
not modify, stop, restart, or select among any of the four active driven runs.
It adds one independently manifested production after those runs finish.

## Frozen equilibrium case

- Chain length: `n=10`.
- Baths: `T_L=T_R=6`, hence `Delta beta=0`.
- Source commit: `1905cf4e606a4a7f4dd8930caa64bd4cc861e9d4`.
- Frozen source: `source/NLS_entropy_ft_1905cf.cpp`.
- Source SHA-256:
  `98e7f8f5f915c8ce02bd8aa10722025c09fd739184b981961692869c9356c0d3`.
- Cartesian, projection-free dynamics; no integrator change.
- `gamma=0.1`, `dt=5e-4`, burn-in `500`.
- Base block duration: `20`; measured bond: `5`.
- 8 SIMD batches x 16 lanes = 128 independent streams.
- 7,813 blocks per stream = exactly 1,000,064 base blocks.
- OpenMP threads: 2.
- Seed: `2026090405`.

The exact command is frozen as

```text
bin/entropy_ft_affinity_equilibrium sample 6 6 10 8 500 20 7813 0.0005 2026090405 2 raw/dbeta_0p000000_n10 5
```

## Scheduling and integrity

The equilibrium run is queued behind the original four-case pipeline so it
does not change their resource allocation.  It starts only on AC power.  The
runner refuses duplicate processes and refuses to overwrite an existing raw
CSV.  Completion requires the exact expected header, 1,000,064 data rows,
finite output, contiguous stream/block identifiers, and a recorded SHA-256.

## Analysis amendment

After equilibrium production is complete, the same frozen aggregation,
matched-bin, raw-count, and stream-bootstrap rules in `PROTOCOL.md` will be
applied.  At `Delta beta=0`, ratios that divide by `Delta beta` are undefined
and must be reported as `N/A`; the equilibrium direct-tail reference is slope
zero.  The final crossover table will replace the formerly unavailable row
with the measured equilibrium moments and negative-count support.  Preliminary
analysis produced before this fifth run must be preserved before the final
six-case analysis is generated.
