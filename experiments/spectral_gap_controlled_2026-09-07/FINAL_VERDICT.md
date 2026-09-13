## Material Passport

- Artifact type: final numerical verdict
- Status: COMPLETE
- Scope: three-mode controlled-integrator relaxation experiment

# Final verdict

**The controlled-integrator experiment does not resolve a unique three-mode
spectral gap.**

The projection-free Cartesian sampler passed all numerical and burn-in audits,
and resolved observables were stable under `dt=1e-3 -> 2.5e-4`.  However, the
unchanged common-plateau gate failed for driven and equal-temperature data at
both timesteps.  Under the frozen interpretation rule, the remaining
observable dependence is an accessible-window property of the dynamics, not a
demonstrated projection artefact and not a timestep-limited result.

The robust positive result is narrower: the damped oscillation in `I1-I3`
survives the controlled discretization with `lambda_I` approximately 5.0--5.5,
consistent with the historical value 5.35.  Other exact exchange-odd
combinations show different frequencies, so 5.35 is not universal across the
odd sector.
