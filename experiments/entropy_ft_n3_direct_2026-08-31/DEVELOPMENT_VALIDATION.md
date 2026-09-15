# Pre-production implementation validation

The endpoint extension was made additively in `flux/NLS_entropy_ft.cpp`.
The original source remains recoverable byte-for-byte from Git commit
`7215d47140683f95173002aade5b0653e72ca0f3`; its SHA-256 is
`98e7f8f5f915c8ce02bd8aa10722025c09fd739184b981961692869c9356c0d3`.

Before committing the prospective protocol:

1. The C++ self-test passed with maximum gradient error
   `1.55486956643e-10`, Hamiltonian energy derivative
   `-9.02131476493e-17`, and maximum boundary-Laplacian error
   `5.29314434244e-4`.
2. The legacy `sample` command was compiled both from the parent-commit source
   and from the endpoint-extended source.  With the same small-run parameters
   and seed, the block, profile, and burn-in files were byte-identical; all
   summary fields except wall time were identical.
3. A 32-row `sample_n3` schema smoke test recomputed
   `Sigma_m=-Q_left/10-Q_right/2` and
   `Q_left+Q_right-Delta E` from every row, and checked exact consecutive
   reduced-endpoint continuity within each stream.  All checks passed.

The smoke data are not retained and are not part of scientific production.
The scientific production uses a distinct frozen seed.
