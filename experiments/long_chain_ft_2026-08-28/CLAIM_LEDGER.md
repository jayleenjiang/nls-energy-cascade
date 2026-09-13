# Claim ledger

| candidate claim | required evidence | current status |
|:---|:---|:---|
| The simulated boundary dynamics obey the intended local path-probability ratio. | Exact kernel-ratio derivation, reversibility self-test, arbitrary-`n` transient Crooks/IFT control. | Self-tests pass; supported one-step controls pass for `n=10,20,40` and the independent `n=30`, `N=5e5` replication. |
| The entropy affinity uses the same temperature convention as the manuscript. | Check the code energy against the Gibbs factor. | PASS: code heat uses `E=H/2`, so `exp(-E/T)=exp(-H/(2T))`; `-0.4 Q_E` is identical to `-0.2 Q_H`. |
| The right-bath entropy-current gauge and medium entropy have the same long-time target. | First-law/gauge identity and decay of their directly sampled finite-time SCGF difference where both have support. | Identity RMS is `6.1e-6`--`1.9e-5`; low-tilt SCGF differences decrease with time. Equality remains conditional on controlled endpoint exponential moments. |
| The NESS SCGF is GC symmetric for `n=10`. | Paired `k=0.3,0.7`, four seeds, support, time convergence, population doubling, numerical controls. | Paired/time/population gates pass through `t=60`.  The first `N_c=512` controls are retained as unsupported diagnostics; the frozen supported `N_c=1024` timestep-halving and selection-interval controls both pass. |
| The NESS SCGF is GC symmetric for `n=20`. | Same gates on a resolved nontrivial pair. | The `k=0.4,0.6` paired/time/population core gates pass through `t=60`; the stronger `0.3,0.7` pair is unresolved.  Because the predeclared endpoint selection controls do not both pass, this positive core result is not promoted to a fully controlled all-chain claim. |
| The NESS SCGF is GC symmetric for `n=30`. | Same gates on `k=0.4,0.6`. | **UNRESOLVED.**  At `N_c=2048`, the paired residual, CI, support, and time gates pass, but the pre-specified 1024-to-2048 `k=0.4` member comparison fails the unchanged two-combined-SE population gate (`-0.00347`, about `2.36` combined SE).  No further population is selected.  The stronger `0.3,0.7` pair is also unresolved. |
| The NESS SCGF is GC symmetric for `n=40`. | Same gates on `k=0.4,0.6`. | **UNRESOLVED under the full control suite.**  The original `t=60/80` series fails the unchanged late-half gate.  The independent remedial `t=120`, `N_c=512/1024` series passes paired/time/population gates and timestep halving, but the supported selection `2 -> 1` comparison fails the unchanged member/two-SE and late-half absolute gates. |
| The model satisfies FT for all tilts or all chain lengths. | Analytic ergodicity/domain proof or reliable SCGF over the full symmetry interval and an asymptotic chain-length argument. | **Not established and must not be claimed.** |

The final admissible wording is: the complete frozen suite is numerically
consistent with GC at `n=10`; the core paired-SCGF diagnostics are consistent
at `n=20` and `n=40` but do not pass the complete cross-chain control suite;
and `n=30` is unresolved by its population gate.  “Proof of FT” and an
all-long-chain verification claim are deliberately excluded.
