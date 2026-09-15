# Final verdict: one-bath open-chain profiles

## Strict stationarity verdict

Five of the six logical conditions pass every predeclared stationarity gate:
`T10_n25`, `T10_n50`, `T6_n25`, `T6_n50`, and `T6_n100`.

`T10_n100` does not pass.  Its total-action and interior-profile gates pass,
and its Q4--Q3 free-end relative change is only 3.21%, but the paired free-end
difference is `-0.0090291 +/- 0.0043472` (SE), giving `z=-2.077`; its 95% CI
therefore narrowly excludes zero.  The difference is negative, not evidence
that the free end is still growing.  The frozen gate is nevertheless failed,
so any three-length result at `T_left=10` is descriptive rather than a fully
stationary scaling claim.

## Profile result

For the fully stationary `T_left=6` set, the midpoint action scales as
`n^-0.47689`, with the three-point, one-residual-degree-of-freedom 95% interval
`[-0.48942,-0.46435]`.  This is close to the fixed `n^-1/2` prediction.

The fixed-exponent rescaling gives a strong whole-profile collapse:

- `T_left=6`: interior median/max spread falls from 64.81%/71.60% in the raw
  profiles to 3.49%/7.58% after multiplying by `sqrt(n)`.
- `T_left=10`: interior median/max spread falls from 64.06%/72.37% to
  4.39%/6.71%; because `T10_n100` misses the strict stationarity gate, this is
  descriptive evidence only.

The predeclared 50% penetration threshold is not reached in any chain: the
action at the free end remains above half of the driven-end value for all six
conditions.  Thus the accessible chains are filled throughout rather than
showing a finite penetration front.

The numerical result is that removing the right reservoir does **not** destroy
the fixed `sqrt(n)` profile collapse over the tested sizes.  A single left BC1
bath, whose drift still contains the global mass `M`, is sufficient to produce
the observed size-dependent action scale.  This statement is restricted to
`n=25,50,100`, the two tested left-bath temperatures, and this numerical
scheme; it is not an asymptotic proof.

## Claim boundaries

- The action-profile result may be used quantitatively for the five conditions
  that pass stationarity, including the complete `T_left=6` size series.
- The `T_left=10` three-size exponent and collapse may be shown only with the
  explicit `T10_n100` stationarity caveat.
- Bond-sine profiles are archival only.  The inherited corrected-noise SIMD
  integrator accumulated 181,474,703 positive-action projection events across
  the 12 runs, and the independent controlled-integrator audit already showed
  that this scheme produces a spurious bond-sine signal.
- There were zero non-finite or discarded trajectories, so instability is not
  the reason for the one failed stationarity gate.
