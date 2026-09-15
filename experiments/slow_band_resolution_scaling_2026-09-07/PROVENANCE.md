# Provenance

- Frozen scientific protocol commit:
  `8efa4c6b6c716d010c0a740f48a8f009426f2869`
- Frozen Part-1 implementation/launcher commit:
  `b9b5f42c8321b4cf3c796737e18b8230f4d302fb`
- Analysis source SHA-256:
  `b2ad9a90fc9ecfd2799c9b778dff3d7194c28faca224da70f20f5db022dce1d4`
- Previous fitting-pipeline SHA-256:
  `ada94924e17f762519a973b28712a3af57ab6f9d1987183a2459c830587d94cb`
- Driven input SHA-256:
  `9a878a637bd7475ab0f05ea0b7991276f4de1fd53f152c274530b7f99a984a47`
- Equilibrium input SHA-256:
  `f7f97e860d75414b3999d656cf585dea9d402fafe7cc1a99a22321e9963c3d45`
- Synthetic master seed: `2026091201`
- Unique deterministic task seeds: 84, listed in
  `part1_analysis/task_seeds.csv`
- Recorded wall time: 99,175.322 seconds (27.55 hours)

Exact command:

```text
python3 experiments/slow_band_resolution_scaling_2026-09-07/run_scaling.py \
  --output experiments/slow_band_resolution_scaling_2026-09-07/part1_analysis \
  --replicates 500 \
  --workers 8
```

Part 1 is analysis-only and generated no physical trajectory. The two saved
fine-timestep empirical correlation arrays are the only data inputs. The
frozen Part-2 gate returned `part2_allowed=false`, so no Part-2 simulator was
launched.

