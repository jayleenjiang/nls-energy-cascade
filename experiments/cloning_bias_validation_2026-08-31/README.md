# Controlled-cloning bias validation

This directory is the prospective repair of the missing cloning safeguards:
a physical-model `n=2` direct reference, four-level `1/N_c` extrapolations,
and the `n=10,t=20,k=0.3` direct-gap check.

The scientific contract is frozen in `PROTOCOL.md`.  `RUN_MATRIX.csv` records
every seed and argument set.  Raw simulator outputs go under `raw/`; generated
tables go under `analysis/`; the final audit and LaTeX report go under
`report/`.

Run order:

```bash
bash experiments/cloning_bias_validation_2026-08-31/build_and_selftest.sh
python3 experiments/cloning_bias_validation_2026-08-31/generate_run_matrix.py
bash experiments/cloning_bias_validation_2026-08-31/run_pipeline_when_ac.sh
```

The pipeline is resumable and skips a run only when both of its expected CSV
outputs are nonempty.  It does not overwrite completed output.
