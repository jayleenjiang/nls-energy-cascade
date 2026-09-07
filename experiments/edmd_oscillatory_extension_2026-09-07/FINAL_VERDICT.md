# Final verdict

**SPECTRUM INCOMPLETE — the frozen oscillatory completeness gate failed.**

The targeted raw-action extension exposes a complex candidate near the known
frequency, and 3/4 bath/timestep cases pass every gate. The driven
`dt=1e-3` case fails dictionary convergence: at `tau=0.05`, E1/E2/E3 give
real parts `-6.4657, -5.1121, -5.5195`; range `1.3536` exceeds the frozen
tolerance `1.1039`.

The new leading visible real rates are approximately `-0.90` to `-0.92` in
the driven runs and `-0.89` to `-0.91` at equilibrium. Three of four new 95%
intervals do not overlap their old intervals near `-0.95`, so the old estimate
was dictionary-sensitive. Because the oscillatory completeness gate failed,
the new value near `-0.91` is still only the slowest rate visible to the tested
dictionaries, not a full-spectrum spectral gap.

The two driven candidates near `-1.8` and `-2.1` remain unresolved because
their trajectory-bootstrap intervals overlap at both timesteps.
