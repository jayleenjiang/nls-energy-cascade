# Entropy cloning audit

Overall: **PASS_WITH_FT_BLOCKED**

| Check | Status | Detail |
|---|---:|---|
| successful_summary_midpoint_integrity | PASS | 0 failures over 73 successful summaries |
| numerical_gate_rejections | BLOCKED | 1 attempted run(s) rejected before summary: high_tilt_convergence/n10_t20_N1024_k0p7_run1.log |
| direct_scgf_audit | PASS | independent direct-SCGF audit status |
| k0_unbiased_recovery | PASS | runs=4, entropy z=0.538, action-current z=0.541 |
| one_interval_generator_identity | PASS | k=0.1: |difference|=0.00290, k=0.5: |difference|=0.01393 |
| low_tilt_direct_crosscheck | PASS | N=512, k=0.1: z=0.80; N=512, k=0.3: z=1.09; N=1024, k=0.3: z=0.49 |
| population_convergence_k0.3 | PASS | N=512->1024: difference=0.00950, 2SE=0.02429, upper support=True |
| selection_interval_convergence_k0.3 | PASS | Delta=0.25 vs 0.5: difference=0.00869, 2SE=0.01145, support=True; Delta=1 vs 0.5: difference=0.00452, 2SE=0.01095, support=True |
| timestep_convergence_k0.3 | PASS | difference=0.00169, 2SE=0.01214, fine support=True |
| observation_time_extension_k0.3 | BLOCKED | t=40 support=False, mean=-0.17607 +/- 0.00719 |
| high_tilt_sampling_support | BLOCKED | no k>=0.5 group passes support |
| gc_pair_ready | BLOCKED | 0 standard-setting supported pairs |

`PASS_WITH_FT_BLOCKED` means that the implemented diagnostics pass their integrity checks, but no fluctuation-theorem claim is accepted because a nontrivial `k <-> 1-k` pair lacks reliable support.
