# Long-chain FT Phase II

This directory contains the prospective validation phase that follows the
frozen Phase-I experiment in `../long_chain_ft_2026-08-28/`.  Phase I is not
rewritten or relabelled.  Phase II addresses only the convergence gates that
remained missing or unresolved at `n=20,30,40`.

The mentor-requested direct-sampling experiment is already complete in
`../entropy_ft_2026-08-26/`: each chain length has 1,000,064 non-overlapping
base blocks of length 20, aggregated to `t=20,40,...,200`, with both-tail,
Gaussian-benchmark, time-scaling, and symmetry analyses.  Those data are
reused rather than simulated again.

## Frozen execution order

1. Build and self-test the sampler:

   ```bash
   bash experiments/long_chain_ft_phase2_2026-08-29/build_and_selftest.sh
   ```

2. Run Stage I controls and the `n=30` population extension:

   ```bash
   bash experiments/long_chain_ft_phase2_2026-08-29/run_phase2.sh stage1
   ```

3. Analyze Stage I.  The conditional `n=30` controls are run only if the
   predeclared population decision passes:

   ```bash
   bash experiments/long_chain_ft_phase2_2026-08-29/analyze_phase2.sh stage1
   bash experiments/long_chain_ft_phase2_2026-08-29/run_phase2.sh n30-controls
   bash experiments/long_chain_ft_phase2_2026-08-29/analyze_phase2.sh final
   ```

The complete protocol and stopping rules are in `PROTOCOL.md`.  Raw simulator
outputs are written under `raw/`; generated audits and figures are written
under `analysis/`.  A failed run is reported and is not automatically retried.
