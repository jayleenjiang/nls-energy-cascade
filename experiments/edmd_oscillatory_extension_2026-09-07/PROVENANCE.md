# Provenance

## Frozen commits

- Protocol freeze: `0af6652`
- Analysis source freeze: `70cbab9`
- Controlled trajectories launch: `899a5d345a98e4c8d5605a6a5a166a08bb86adcd`
- Original EDMD result: `8ff676b`

## Source hashes

- `run_extended_edmd.py`: `a986b0cc1d3f0e190223d348d7522b2ac81ce307e0628fd00852b2329aee2c57`
- imported `run_edmd.py`: `d4c53270ed902cbcf6692dcc64d7c1f423e40189504d03703c6e4f8cd50cc6ef`
- controlled trajectory simulator source: `58c871882f1180701a7aae637a8f4be113315a72fb9d7815562f5fc90faf7976`
- controlled trajectory binary: `f7cd820bb8f716a8790f688db060259da5657e4f20fa09100546e66efb696b3f`

## Exact commands

```sh
python3 experiments/edmd_oscillatory_extension_2026-09-07/run_extended_edmd.py \
  --input experiments/spectral_gap_controlled_2026-09-07/raw \
  --output experiments/edmd_oscillatory_extension_2026-09-07/analysis

python3 experiments/edmd_oscillatory_extension_2026-09-07/summarize_extended_edmd.py \
  --analysis experiments/edmd_oscillatory_extension_2026-09-07/analysis

python3 experiments/edmd_oscillatory_extension_2026-09-07/make_figures.py
python3 experiments/edmd_oscillatory_extension_2026-09-07/audit_results.py
```

## Seeds

- Controlled trajectory base seeds: 2026090801, 2026090802, 2026090803,
  2026090804; 64 independent stream seeds per case are in the original
  controlled-run provenance.
- Oscillatory bootstrap base seed: 2026090702 plus case index.
- Real-mode bootstrap base seed: 2026090702 plus 100 plus case index.

Full hashes are in `provenance/input_stream_hashes.sha256` and
`provenance/analysis_output_hashes.sha256`.
