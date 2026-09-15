# Package contents

This directory contains the complete `section4-density-ratio-v1` analysis:

- frozen protocol and deterministic stream-role manifest;
- source code for the Gibbs-anchored models, training and selectors;
- all 9 equilibrium and 6 driven model bundles, histories and raw metrics;
- frozen equilibrium validation set and one-shot blind result;
- frozen driven validation verdict and raw moment/marginal CSV files;
- exact commands, source/data hashes and the verified-cache amendment;
- Markdown verdict/validation reports and the compiled LaTeX/PDF report.

The 256 input trajectories are intentionally not duplicated here.  They are in
the previously delivered package
`nls_short_chain_section4_2026-09-14.tar.zst`; every file used here is identified
by path and SHA-256 in `FROZEN_STREAM_ROLES.csv`, and the local extraction audit
is in `provenance/CACHE_AUDIT.json`.

The driven blind split was not opened and no coarse-timestep density-ratio fit
was run because the preregistered driven validation gate failed.

