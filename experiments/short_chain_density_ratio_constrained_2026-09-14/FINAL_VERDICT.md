# Final verdict: five-dimensional three-mode NESS density

PASS, with a strict sampled-support claim.

For the driven three-mode system at `(T1,T3)=(2,8)`, `gamma=0.1`, the study
constructs an absolutely normalized numerical density in the reduced variables
`(I1,I2,I3,theta1,theta3)` with respect to
`dI1 dI2 dI3 dtheta1 dtheta3`.  The final fine-timestep estimator passed all
predeclared validation and the declared final test gates on a fresh 64-stream
holdout at `dt=2.5e-4`.  Fine R11 and coarse R10 now use the identical
three-moment exponential calibration and fixed 0.625 shrink factor, with
unchanged neural weights.  Both pass validation, test, and equilibrium
known-answer controls.  Their centered ensemble log-density RMS is 0.024355 on
the fine test support and 0.024353 on the coarse test support, below the frozen
0.10 gate.

The density may therefore be used as a numerical solution of the full reduced
five-dimensional NESS on the region sampled by the stationary trajectories.
This is not an analytic global solution, and accuracy in unsampled extreme
tails is not established.

The machine-readable evaluator is `scripts/final_ness_density.py`; final
tables are in `final_analysis/`, and publication figures are in
`final_figures/`.
