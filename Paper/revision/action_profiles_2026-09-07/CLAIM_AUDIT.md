# Claim audit for revised Section 3.1 action profiles

## Supported directly by the 36-profile dataset

- Reversing `(T1,Tn)` reverses the spatial action gradient.
- BC1 and BC2 give nearly coincident one-point action profiles.
- BC3 and BC3b sustain order-one midpoint action over `n=25,50,100`.
- The three-point midpoint slopes lie in `[-0.485,-0.464]` for BC1/BC2 and
  `[0.032,0.083]` for BC3/BC3b across the three bath pairs.
- All 72 replicate runs completed with no non-finite or discarded trajectory.
- The BC1 `(10,2)` endpoint-reproduction gate passed all six endpoints.

## Deliberately not claimed

- The three-point slopes are not claimed as asymptotic exponents.
- The residual standard errors of those fits are not presented as asymptotic
  confidence intervals.
- Equal-temperature SIMD data are not claimed to verify the exact Gibbs
  measure: the bond-sine control failed systematically, and the BC1/BC2 action
  profiles have a resolved spatial bias at `n=100`.
- The boundary fixed-point argument is presented as a mechanism consistent
  with the data, not as a theorem deriving the observed profiles.

