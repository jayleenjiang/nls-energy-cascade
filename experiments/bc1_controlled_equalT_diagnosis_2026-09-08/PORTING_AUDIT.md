# Porting audit

The profile sampler is derived from the verified controlled spectral sampler:

- parent source SHA-256:
  `58c871882f1180701a7aae637a8f4be113315a72fb9d7815562f5fc90faf7976`;
- parent binary SHA-256:
  `f7cd820bb8f716a8790f688db060259da5657e4f20fa09100546e66efb696b3f`.

The parent and ported executables were run at `n=3`, `T1=T3=5`,
`gamma=0.1`, `dt=5e-4`, burn-in `0.01`, measurement `0.02`, sampling interval
`0.01`, one 16-lane batch, base seed `2026090800`, and stream offset zero.
The port's final action values agree with the parent's float32 saved values to
maximum absolute error `1.1920929e-07`. After accounting for the parent's
left-bond orientation (`2(phi1-phi2)` versus the profile convention
`2(phi2-phi1)`), the final bond-sine values agree exactly at float32 precision.

Generic-chain smoke tests at `n=25` and `n=100` produced finite output with
zero projection, floor, midpoint-failure, zero-radius, and non-finite events.
These smoke runs are not production data and are excluded from analysis.
