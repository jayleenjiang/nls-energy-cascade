# Validation and integrity report

1. The frozen protocol was committed as `c7e2a92` before production output.
2. The simulator source was unchanged from SHA-256
   `9ae5835ed708c8794c8b00ba799b23761482953aaf0ed47cd0b4ba3966d4eaf2`.
3. The binary SHA-256 was
   `93c4aa6d046f3c35cd0ea1136091fc44ef1b64baf6237e4663ca9e86b995a156`.
4. The source self-test passed: maximum gradient error `1.55487e-10`,
   Hamiltonian `dE/dt=-9.02131e-17`, and maximum boundary-Laplacian error
   `5.29314e-4`.
5. Both Zstandard archives pass `zstd -t`.
6. Each decompressed CSV has `1,000,065` lines: one header and exactly
   `1,000,064` block rows.
7. The analysis independently verifies 128 streams and 7,813 blocks per
   stream for each temperature, with all required values finite.
8. Both production summaries report zero implicit-midpoint failures.
9. No fit window was tuned.  Both cases use all 80 frozen-rule symmetric bin
   pairs; every pair exceeds 200 observations on each side.
10. The first-law block residual RMS is `2.362029704304115e-4` at `T=6` and
    `7.432934102765286e-4` at `T=10`, corresponding to rate RMS values
    `1.1810148521520336e-5` and `3.7164670513826657e-5`.
11. All four requested means at both temperatures include zero in the frozen
    stream-bootstrap 95% CI.  Both entropy means are within 3 stream-level
    standard errors, and both frozen symmetry-slope CIs include zero.

Verdict: **PASS under the predeclared gates.**
