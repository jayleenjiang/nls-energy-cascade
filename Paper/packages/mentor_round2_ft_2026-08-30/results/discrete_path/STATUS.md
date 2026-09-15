# FT verification status

Overall: **PASS within the declared finite-time n=2 scope**.

## Evidence 1: exact discrete path ratio

- `n=2`, `(T_left,T_right)=(10,2)`, `t=0.1`, `dt=1e-3`.
- 1,000,000 independent forward and 1,000,000 independent reverse paths.
- Forward log-IFT `-0.003532`, 95% CI `[-0.007030,0.000244]`.
- Reverse log-IFT `-0.002145`, 95% CI `[-0.007289,0.004550]`.
- Raw-bin Crooks slope `0.993352 +/- 0.00243`, intercept
  `0.000259 +/- 0.00151`, 33 accepted bins.
- Kernel entropy approaches heat entropy at first order under timestep
  refinement.

Authoritative audit:
`../discrete_path_ft_2026-08-28/production_v2/driven_t0p1_dt1e3_N1m_analysis/audit.md`.

## Evidence 2: driven NESS total entropy

- `n=2`, `(T_left,T_right)=(10,2)`, burn-in `500`, `t=0.1`, `dt=5e-4`.
- 64 independent streams and 1,048,576 non-overlapping blocks.
- Stationary endpoint density cross-fitted in reduced coordinates and fixed
  before driven application.
- Independent Gibbs validation: support `0.99979`, RMSE `0.00826`, slope
  `0.99941`, correlation `0.99994`.
- Driven log-IFT `0.002175`, 95% CI `[-0.000744,0.005251]`, ESS `324,604`.
- Same-process detailed-FT slope `0.99450 +/- 0.00738`, intercept `-0.02219`,
  60 accepted raw bins.
- Three action/Fourier density families give slopes from `0.9945` to `1.0035`
  and all IFT intervals include zero.

Authoritative audit:
`../entropy_ft_scgf_2026-08-27/total_entropy_n2_short/parametric_analysis/audit.md`.

## Allowed claim

> Finite-time integral and detailed fluctuation relations are numerically
> verified for total entropy in the projection-free Cartesian NLS dynamics at
> n=2, both through an exact discrete forward/reverse path-ratio construction
> and through an independently validated driven-NESS endpoint calculation.

## Claim boundary

This is not a mathematical proof.  It does not establish the asymptotic
same-process Gallavotti--Cohen symmetry for `n=10,20,30,40`.  The long-chain
medium-entropy experiment remains endpoint-term and rare-event limited.
