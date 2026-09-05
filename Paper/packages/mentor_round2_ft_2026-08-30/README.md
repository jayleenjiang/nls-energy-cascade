# Mentor round-II fluctuation-theorem package

Date: 2026-08-30

This is the compact, self-contained delivery package for the mentor's second
round of fluctuation questions.  It contains the detailed report, the source
and analysis code used for the accepted calculations, audited result tables,
figures, protocols, compact Phase-II raw outputs, and integrity manifests.

## Start here

1. `mentor_round2_ft_report.pdf` — detailed scientific report organized by the
   mentor's original questions.
2. `REPRODUCE.md` — build, analysis, and audit commands.
3. `RAW_DATA_INDEX.md` — locations, sizes, and SHA-256 hashes of the large raw
   files that remain in the main local repository.
4. `results/` — accepted compact data and audits.
5. `code/` — the exact simulation and analysis programs needed for each part.

## Scientific bottom line

- The requested direct study is complete: each of `n=10,20,30,40` has
  1,000,064 non-overlapping `t=20` blocks, aggregated to
  `t=20,40,...,200`.
- The action-current tails show resolved finite-time large-deviation scaling,
  but direct medium-entropy and heat-current histograms do not reach their
  asymptotic FT references before the negative tail becomes unresolved.
- At `n=2`, two independent total-entropy constructions pass finite-time
  detailed and integral fluctuation-relation audits.
- At `n=10,20,30,40`, controlled cloning estimates are numerically consistent
  with `psi_n(k)=psi_n(1-k)` at the one resolved complementary pair tested for
  each length.  This is not an all-tilt proof or an `n -> infinity` theorem.

## Package layout

```text
code/
  direct_sampling/       million-block sampler and two-tail analyses
  n2_total_entropy/      endpoint-density total-entropy analyses
  discrete_path/         exact finite-step path-ratio control
  long_chain_scgf/       cloning simulator and Phase-II analyses
figures/                 figures used by the detailed report
protocols/               direct-sampling protocol and production manifest
results/
  direct_sampling/       audited compact products of the 4M-block run
  n2_total_entropy/      accepted parametric endpoint analysis
  discrete_path/         accepted path-ratio analysis and timestep control
  long_chain_phase2/     protocol, all 56 compact raw runs, analyses, verdict
```

The multi-hundred-megabyte trajectory/block CSVs are not duplicated inside
this compact package.  Their immutable hashes and exact local paths are in
`RAW_DATA_INDEX.md`; the accepted Phase-II raw summaries and timeseries are
small and are included in full.
