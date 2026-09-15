#!/usr/bin/env python3
"""Materialize all frozen target/reference stream roles before fitting."""

from __future__ import annotations

import csv
from pathlib import Path

from density_ratio_common import target_reference_rows


HERE = Path(__file__).resolve().parent.parent
records = []
for phase in ("equilibrium", "driven"):
    for timestep in ("fine", "coarse"):
        for split in ("train", "validation", "test"):
            target, reference = target_reference_rows(phase, timestep, split)
            for role, selected in (("target", target), ("reference", reference)):
                for row in selected:
                    records.append({
                        "protocol_version": "section4-density-ratio-v1",
                        "phase": phase, "timestep": timestep, "split": split,
                        "role": role, "stream_id": row["stream_id"],
                        "physical_case": row["case"],
                        "relative_path": row["relative_path"], "sha256": row["sha256"],
                    })
path = HERE/"FROZEN_STREAM_ROLES.csv"
with path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=records[0].keys())
    writer.writeheader(); writer.writerows(records)
print(f"wrote {len(records)} rows to {path}")

