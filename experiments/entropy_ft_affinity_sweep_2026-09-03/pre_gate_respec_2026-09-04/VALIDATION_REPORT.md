# Validation report: n=10 heat-flux affinity sweep

Final analysis date: 2026-09-04

## Integrity and provenance

- The four new driven files and the new equilibrium file each contain exactly
  1,000,064 data rows, with 128 streams and 7,813 consecutive block IDs per
  stream.
- All five files have the expected nine-column CSV schema and finite selected
  heat fields.
- The accepted (10,2) production was reused without modification.
- Production source commit:
  1905cf4e606a4a7f4dd8930caa64bd4cc861e9d4.
- Production source SHA-256:
  98e7f8f5f915c8ce02bd8aa10722025c09fd739184b981961692869c9356c0d3.
- Driven protocol commit: e170205c5d9a3d9b2bfca3714a4447fa965a1e40.
- Equilibrium amendment commit:
  cf5282d3ab77ce93c05578ff15fb6598026081c7.
- The preliminary analysis made before equilibrium was available is preserved
  locally in pre_equilibrium_outputs/.

Raw input SHA-256 values:

| case | rows | SHA-256 |
|---|---:|---|
| (6,6) | 1,000,064 | 7fd6827f12b42844de1142613af124702ca5535ba63f3c54ab71006b256aef4d |
| (6.5,5.5) | 1,000,064 | 8b44be37ea8f566fb3862b6cd796dcfcf433be0888752babb876754c84d4c421 |
| (7,5) | 1,000,064 | 382496fcd194154b0f6d486f2d4da9de80d2d0791f94cb0f3e695d8e666a29cc |
| (8,4) | 1,000,064 | 0404492717ce335c1e2eff7970d5d514435fabec9b91604fa92b7d4846cdf685 |
| (9,3) | 1,000,064 | 342f3783b40db98263cd685c8ab4e531c9d83842db6762b0bd0a06061642958d |
| (10,2) reused | 1,000,064 | a23806e82f5514a9c3375d10a6644946b6b57d0efe8452b3bbf397b6230f9929 |

## Frozen-gate verdict

Every nonzero-affinity case is **UNRESOLVED** under the predeclared final gate.
This is not recorded as a pass or a failure.

| Delta beta | resolved full-sample times | full-sample a_inf | a_inf / Delta beta | valid joint bootstraps / 1000 | verdict |
|---:|---|---:|---:|---:|---|
| 0.027972 | 20, 40, 80, 160, 320, 640 | 0.0284613 | 1.01749 | 0 | UNRESOLVED |
| 0.057143 | 20, 40, 80, 160, 320 | 0.0580178 | 1.01531 | 0 | UNRESOLVED |
| 0.125000 | 20, 40, 80 | 0.127749 | 1.02199 | 0 | UNRESOLVED |
| 0.222222 | 20, 40 | N/A | N/A | 0 | UNRESOLVED |
| 0.400000 | 20 | N/A | N/A | 0 | UNRESOLVED |

The first three full-sample intercepts are close to the FT references. They
cannot be promoted to a numerical FT confirmation because the fixed
stream-bootstrap construction did not yield the required 800 jointly resolved
intercepts. No bins, times, or bootstrap rules were changed after seeing this
outcome.

## Equilibrium control

At Delta beta = 0, the reference asymmetry slope is zero. All six full-sample
times are directly resolved. The full-sample extrapolated intercept is
7.80343e-05, but its joint bootstrap also has 0/1000 accepted replicates, so
the formal long-time equilibrium verdict is **UNRESOLVED**.

At t=20, where the individual bootstrap passes the 800-replicate rule:

- direct slope: 1.63012e-04;
- WLS SE: 1.28209e-04;
- stream-bootstrap 95% CI: [-8.97958e-05, 4.18259e-04];
- accepted bootstrap replicates: 818/1000;
- negative windows: 499,786/1,000,064.

The interval contains zero and the signs are nearly balanced, so this
individual finite-time control is consistent with equilibrium. It does not
repair the missing joint long-time CI.

## CLT and Gaussian crossover

The Gaussian bulk diagnostic approaches the FT reference at weak drive but
departs from it as the affinity grows:

| Delta beta | a_Gauss / Delta beta at largest resolved two-tail t | gaussFT at t=640 | skew at t=160 | excess kurtosis at t=160 | n_neg at t=160 |
|---:|---:|---:|---:|---:|---:|
| 0.027972 | 0.9599 | 1.0418 | 0.0512 | 0.0388 | 20,136 |
| 0.057143 | 0.9310 | 1.0659 | 0.0951 | 0.0409 | 2,997 |
| 0.125000 | 0.7450 | 1.3175 | 0.1806 | 0.0808 | 2 |
| 0.222222 | 0.5147 | 1.9406 | 0.2323 | 0.0889 | 0 |
| 0.400000 | 0.3007 | 3.2298 | 0.2361 | 0.0836 | 0 |

Thus the directly sampled negative tail remains useful at weak affinity and
collapses rapidly at stronger drive. The Gaussian slope is a bulk/CLT
diagnostic and is not substituted for missing two-tail evidence.

## Artifact map

- analysis/window_summary.csv: all 36 case/window moment and fit rows.
- analysis/symmetric_bin_raw_counts.csv: all 1,031 raw positive/negative bin
  pairs used by the fixed fits.
- analysis/infinite_time_extrapolation.csv: the six extrapolation records.
- analysis/crossover_summary.csv: the requested six-row crossover table.
- analysis/input_hashes.csv: absolute input paths, row counts, and hashes.
- figures/: publication-format PNG and PDF figures.
- report/affinity_sweep_report.tex and .pdf: final report.
- COMMANDS.tsv and EQUILIBRIUM_COMMAND.tsv: exact production commands.

## Claim boundary

The supported statement is: the full-sample direct-tail slopes move toward the
FT reference as the drive weakens, while the Gaussian bulk ratio approaches
one and the negative tail becomes increasingly resolvable. Under the frozen
bootstrap gate, this sweep does **not** establish the long-time heat-flux FT;
all formal long-time verdicts remain UNRESOLVED.
