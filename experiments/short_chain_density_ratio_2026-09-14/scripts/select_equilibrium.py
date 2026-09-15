#!/usr/bin/env python3
"""Freeze equilibrium-admissible density-ratio architectures."""

from __future__ import annotations

import json
from pathlib import Path


HERE = Path(__file__).resolve().parent.parent
architectures = ("linear", "mlp32", "mlp64")
seeds = (5201, 5202, 5203)
summaries = []
admissible = []
for architecture in architectures:
    entries = []
    for seed in seeds:
        path = HERE/"equilibrium_candidates"/f"{architecture}_seed{seed}"/"validation_metrics.json"
        if not path.exists():
            raise SystemExit(f"missing {path}")
        entries.append(json.loads(path.read_text()))
    passed = all(item.get("validation_admissible", False) for item in entries)
    summaries.append({"architecture": architecture, "all_seed_admissible": passed,
                      "seeds": entries})
    if passed:
        admissible.append(architecture)
result = {
    "protocol_version": "section4-density-ratio-v1",
    "status": "FROZEN" if admissible else "NO_ADMISSIBLE_CANDIDATE",
    "admissible_architectures": admissible,
    "candidate_summaries": summaries,
}
path = HERE/"EQUILIBRIUM_ADMISSIBLE_SET.json"
path.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
print(json.dumps(result, indent=2, sort_keys=True))
raise SystemExit(0 if admissible else 2)

