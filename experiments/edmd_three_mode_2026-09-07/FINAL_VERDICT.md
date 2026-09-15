# Final verdict

The saved controlled trajectories were sufficient to reconstruct the complete
reduced state, so no new simulation was run.

The EDMD calculation passes its basic numerical controls: every full-sample
fit returns the constant eigenfunction with `|mu-1| < 1.2e-14`; all 500 Gram
directions are retained in D3; changing the fixed relative cutoff from `1e-8`
to `1e-12` changes none of the reported spectra at printed precision.

For the driven system, exactly one nontrivial rate survives the frozen
dictionary, lag, timestep, bootstrap, and neighbour-separation gates:

- `dt=1e-3`: lambda = -0.940174, stream-bootstrap 95% CI
  [-0.954628,-0.920824];
- `dt=2.5e-4`: lambda = -0.963191, 95% CI
  [-0.980007,-0.941942].

The next two driven candidates, near -1.8 and -2.1, have overlapping
trajectory-bootstrap intervals and therefore are not separately reportable.
The known oscillation near `Im(lambda)=5.3` appears at selected lags but fails
the lag-convergence gate in all four datasets.  EDMD therefore does not pass
that known-answer check.

At equilibrium, the leading slow rate is again stable (about -0.95).  A second,
faster real rate near -2.85 also passes all gates at both timesteps.  This does
not validate an analogous second driven rate.

Consequently, the requested hypothesis is **not established**: the present
EDMD dictionary does not resolve several distinct driven slow modes accounting
for the four autocorrelation plateau values.  It resolves one common leading
slow mode.  Observable weights on that mode differ strongly (0.49 for I2,
0.28 for I1+I3, 0.12 for cos(theta3), and 0.22 for the cosine sum), while the
remaining weight is distributed over faster but non-separated candidates.
This is consistent with finite-window modal contamination, but is not enough
to claim a resolved multi-mode spectrum.

