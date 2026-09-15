# Final provenance

- parent repository commit before final integration:
  `995dcdc643d9ef147bca1001f727f9cd010f6ad4`;
- branch: `codex/paper-journal-revision`;
- validated simulator source SHA-256:
  `0d7000b628257311515e8edd43492a348c9bdbc3903c95c06bff5c718bf1f8a2`;
- validated simulator binary SHA-256:
  `58c22764eb864f156bc1df4a0190a0ac1851ff2505fda5fe07b9fde2e4469c82`;
- fine driven/equilibrium holdout seeds: `2026091701`, `2026091702`;
- final fine calibration base seeds: `8201`, `8202`, `8203`;
- final fine reported model seeds: `47201`, `47202`, `47203`;
- R11 uses no new simulation and no neural retraining.

Exact original commands are in `provenance/FINAL_DENSITY_COMMANDS.txt`; R11
commands are in `R11_COMMANDS.txt`.  Raw state
paths and hashes are in `V4_HOLDOUT_SPLITS.csv` and `V6_HOLDOUT_SPLITS.csv`.
The outer ZIP hash is recorded in the adjacent `.sha256` file because an
archive cannot contain its own stable digest.
