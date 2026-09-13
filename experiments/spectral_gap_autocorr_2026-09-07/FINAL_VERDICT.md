# Final verdict

The frozen stationary-autocorrelation experiment does **not** resolve a common
spectral gap for the three-mode generator.

For the driven `(T1,T3)=(2,8)` case, only two observables pass the predeclared
plateau rule: `I2` gives `lambda_R=-0.9926` with trajectory-bootstrap 95% CI
`[-1.0517,-0.9337]` on `t=[1.12,2.37]`, while `I1+I3` gives `-1.1088`
`[-1.2013,-1.0228]` on `t=[1.00,1.58]`.  The rule requires at least three
agreeing observables, so no generator-wide value is quoted.

For the equal-temperature `(5,5)` reference, `cos(theta1)`, `cos(theta3)`, and
`I2` pass separately, but their estimates span `-1.188` to `-0.907`; their 95%
confidence intervals have no common intersection and the relative spread
exceeds 20%.  The common-plateau gate therefore also fails.

All 14 stationarity checks pass by a wide margin, all raw files have the exact
expected size and finite values, and effective-count estimates at accepted
plateau endpoints exceed 50,000.  Thus the negative result is not explained by
failed burn-in or a handful of effective samples.  The correlations instead
show observable-dependent superpositions and loss of signal before a common
late-time mode is isolated.

No damped fit has a resolved nonzero frequency under the half-cycle and AIC
identifiability rule.  This is consistent with the earlier statement that the
resolved late spectrum is essentially real, while not establishing a unique
real gap.

The two driven observable-specific rates are closest to the historical
`-0.934` full-window diagnostic.  They do not support the early-window/Prony
values near `-1.60/-1.65`, and they are not consistent with treating `-0.68`
as a common late-time plateau.

