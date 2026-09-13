# Reproducibility audit

The frozen analysis command was run again after production completion:

```text
python3 analyze_multiaffinity_timestep.py --experiment-dir \
  /Users/jayleenjiang/Documents/NLS/experiments/entropy_ft_multiaffinity_timestep_2026-09-11
```

Result: every CSV and JSON numerical output reproduced byte for byte.  Every
PNG figure also reproduced byte for byte.  Only the four Matplotlib-generated
vector PDF hashes changed; their paired PNG renders were unchanged, and the
PDF difference is attributable to regenerated creation metadata rather than
numerical or graphical content.

This verifies deterministic analysis, fixed whole-stream bootstrap seeds, and
stable numerical results from the immutable raw inputs.
