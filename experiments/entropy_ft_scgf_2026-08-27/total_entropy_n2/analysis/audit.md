# n=2 total-entropy endpoint audit

Overall: **BLOCKED**

| Check | Status | Detail |
|---|---:|---|
| equilibrium_numerical_integrity | PASS | rows=32768/32768, nonfinite=0, midpoint failures=0, formula max=4.84e-13, balance RMS rate=1.24e-05 |
| driven_numerical_integrity | PASS | rows=32768/32768, nonfinite=0, midpoint failures=0, formula max=3.08e-12, balance RMS rate=1.51e-05 |
| equal_temperature_density_ratio | PASS | support=0.9868, slope=0.9296, correlation=0.9675, RMSE=0.3603 |
| learned_equilibrium_integral_ft | BLOCKED | log-IFT CI=[0.0700,0.1025], ESS=6661.5, streams=32 |
| exact_gibbs_balance_control | PASS | log mean exp(-Delta s_total)=-3.293e-07 |
| driven_endpoint_support | BLOCKED | support=0.9868, ESS=1.6, streams=64 |
| endpoint_stationarity | PASS | maximum |stream-level z|=2.359 over 10 checks |

Passing would validate only the exploratory n=2 endpoint estimator;
it would not establish a long-chain detailed fluctuation theorem.
