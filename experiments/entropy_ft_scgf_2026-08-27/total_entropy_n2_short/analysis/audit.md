# n=2 total-entropy endpoint audit

Overall: **BLOCKED**

| Check | Status | Detail |
|---|---:|---|
| equilibrium_numerical_integrity | BLOCKED | rows=262144/262144, nonfinite=0, midpoint failures=0, formula max=7.99e-15, balance RMS rate=2.79e-04 |
| driven_numerical_integrity | BLOCKED | rows=262144/262144, nonfinite=0, midpoint failures=0, formula max=1.78e-14, balance RMS rate=2.95e-04 |
| equal_temperature_density_ratio | BLOCKED | support=0.9980, slope=0.8535, correlation=0.9461, RMSE=0.2427 |
| learned_equilibrium_integral_ft | BLOCKED | log-IFT CI=[0.0315,0.0368], ESS=101686.5, streams=32 |
| exact_gibbs_balance_control | PASS | log mean exp(-Delta s_total)=-2.281e-09 |
| driven_endpoint_support | PASS | support=0.9980, ESS=15860.2, streams=64 |
| endpoint_stationarity | PASS | maximum |stream-level z|=0.983 over 10 checks |

Passing would validate only the exploratory n=2 endpoint estimator;
it would not establish a long-chain detailed fluctuation theorem.
