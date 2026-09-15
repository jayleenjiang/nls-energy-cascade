# Final verdict

**PASS: the projection-free controlled integrator satisfies the BC1 equal-temperature known-answer control at n=25 and n=100.**

- n=25: max |z| = 1.875311; 0/24 pointwise failures; Bonferroni critical value 3.078088; simultaneous PASS.
- n=100: max |z| = 2.984038; 5/99 uncorrected pointwise failures; Bonferroni critical value 3.478063; simultaneous PASS.
- Old SIMD maxima were 8.678126 and 12.928525, with 19/24 and 71/99 pointwise failures.
- The action-profile porting gate passed: spatial-mean differences were 0.407% and 1.298%.
- Projection, floor, zero-radius, midpoint-failure, and nonfinite-batch counts were all zero.

The inherited SIMD equal-temperature sine signal is therefore a numerical-scheme artefact. This comparison changes several numerical ingredients together, so it does not isolate projection as the unique cause. Driven (10,2) SIMD bond-sine profiles must be re-measured with the controlled integrator before quantitative use.
