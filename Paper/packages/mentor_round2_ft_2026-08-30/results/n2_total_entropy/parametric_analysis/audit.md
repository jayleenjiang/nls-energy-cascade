# Parametric n=2 total-entropy FT audit

Overall: **PASS_SMALL_CHAIN_NESS_TOTAL_ENTROPY_FT**

| Check | Status | Detail |
|---|---:|---|
| equilibrium_numerical | PASS | rows=1048576/1048576, nonfinite=0, midpoint failures=0, formula max=1.33e-14, absolute balance RMS=2.76e-05 |
| driven_numerical | PASS | rows=1048576/1048576, nonfinite=0, midpoint failures=0, formula max=1.51e-14, absolute balance RMS=3.06e-05 |
| independent_equilibrium_density_validation | PASS | support=0.99979, RMSE=0.00826, slope=0.99941, correlation=0.99994 |
| density_optimizer | PASS | equilibrium/heldout: success=1, |grad|=9.74e-07; driven/0: success=1, |grad|=7.92e-07; driven/1: success=1, |grad|=7.76e-07 |
| equilibrium_endpoint_error_budget | PASS | log-IFT=3.41299e-05, CI=[3.33951e-05,3.50088e-05], ESS=524140.3, support=0.99979 |
| driven_ift | PASS | log-IFT=0.00217456, CI=[-0.000744111,0.00525109], ESS=324604.3, support=0.99973 |
| exact_gibbs_control | PASS | log-IFT=4.899e-09 |
| driven_detailed_ft | PASS | bins=60, slope=0.994496+/-0.00738, intercept=-0.0221896 |
| density_model_sensitivity | PASS | D2K1: IFT CI=[-0.0007441,0.005251], slope=0.9945; D2K2: IFT CI=[-0.0002062,0.005535], slope=0.9977; D3K2: IFT CI=[-0.001452,0.005606], slope=1.0035 |
| endpoint_stationarity | PASS | maximum |stream-level z|=0.853 |

The accepted scope is a finite-time, n=2 NESS total-entropy numerical
verification. It is not a mathematical proof or a long-chain GC result.
