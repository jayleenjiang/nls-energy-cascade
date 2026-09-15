# Small-chain total-entropy endpoint control

## Purpose

The exact finite-time stochastic-thermodynamic relation concerns

`Delta s_total = Sigma_medium + log rho_ss(X_0) - log rho_ss(X_t)`.

The million-sample production files contain `Sigma_medium` but not the NESS
density endpoint.  This control asks whether that endpoint can be reconstructed
reliably for `n=2`.  It is intentionally separate from the long-chain action
current and from the medium-entropy cloning calculation.

## Why `n=2`

The Cartesian state is four-dimensional, but global phase is uniform.  The
stationary density can therefore be represented in the three reduced
coordinates

`u_1 = log |c_1|^2`, `u_2 = log |c_2|^2`, and
`theta = 2(arg c_2 - arg c_1) mod 2 pi`.

If `q(u_1,u_2,theta)` is the density with respect to
`du_1 du_2 dtheta`, the Cartesian density ratio is obtained from

`log rho = log q - u_1 - u_2 + constant`.

The Jacobian term is required; omitting it would bias the system entropy.

## Data and cross-fitting

The endpoint sampler stores non-overlapping trajectory blocks from independent
streams after burn-in.  It records reduced start/end states, start/end energy,
both bath heats, medium entropy, action current, and the first-law residual.

Density estimation is cross-fitted by independent stream parity:

- fit on even streams and evaluate odd-stream endpoints;
- fit on odd streams and evaluate even-stream endpoints;
- never evaluate an endpoint with a density fitted to its own stream.

The density is a smoothed periodic three-dimensional histogram.  Candidate
smoothing scales are selected on an equal-temperature data set only.

## Exact equilibrium validation

At `T_left=T_right=T`, the exact stationary density is Gibbs and

`log rho(X_0)/rho(X_t) = [E(X_t)-E(X_0)]/T`.

This supplies a non-circular validation target for the learned density ratio.
The driven `T_left=10, T_right=2` data are not used to choose the histogram
smoothing.

The predeclared acceptance gates are in `PREDECLARED_GATES.md`.  In
particular, the learned equilibrium endpoint must pass support, slope,
correlation, RMSE, integral-FT, independent-stream, and exponential-weight ESS
checks before the driven result is interpreted.

## Interpretation boundary

Passing this control would provide an `n=2` numerical consistency test of total
entropy.  It would not prove the theorem and would not justify replacing the
unknown `2n`-dimensional NESS density at `n=10,20,30,40`.  Failure would mean
that nonparametric endpoint-density error is too large, not that the FT is
violated.
