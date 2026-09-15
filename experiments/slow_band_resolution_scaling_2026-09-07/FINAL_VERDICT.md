# Final verdict

Part 1 completed all **336 frozen cells** and **168,000 raw fit replicates**.
The requested separation, `delta=0.20`, failed the frozen resolution rule in
all 48 tested cells, including both bath cases, all six stream
counts through `N=2048`, both windows, and both amplitude constraints.

The Part-2 gate is therefore **closed**: `delta=0.20 did not pass for both cases at a common tested N<=2048`. No real trajectory
production was launched. This is the protocol-mandated result, not an
interrupted or incomplete run.

There were 22 passing cells elsewhere on the grid. They occur almost
entirely at the much wider separation `delta=0.60`; the driven fixed,
non-negative rule reaches `delta=0.40` only at `N=2048`. Thus the current
inverse problem can resolve coarse spectral splitting, but not the physically
relevant approximately 0.20 slow-band spacing.

For `delta=0.20`, 7 of 8 rule/case error curves meet the frozen
saturation definition. The fitted error exponents are much shallower than the
independent-stream `N^(-1/2)` benchmark and the extrapolated stream counts are
unstable, ranging from about `3.5e4` to effectively infinite. These are
descriptive extrapolations only and cannot authorize Part 2.

The defensible claim is a numerical resolution bound: with the existing
empirical covariance, lag support, and frozen multi-exponential pipeline,
`delta=0.20` is not recoverable with up to 2048 emulated independent streams.
This does not prove that the physical slow band is absent.
