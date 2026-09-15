#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent.parent
SPLITS = HERE / "V6_HOLDOUT_SPLITS.csv"


def rows(case, role):
    with SPLITS.open(newline="") as handle:
        selected = [row for row in csv.DictReader(handle)
                    if row["case"] == case and row["role"] == role]
    return sorted(selected, key=lambda row: int(row["stream_id"]))


def evaluation_roles(phase, timestep, split):
    if timestep != "fine" or split not in ("validation", "test"):
        raise ValueError((timestep, split))
    if phase == "driven":
        return rows("driven_fine", split), rows("equilibrium_fine", split)
    if phase == "equilibrium":
        selected = sorted(
            rows("equilibrium_fine", split),
            key=lambda row: hashlib.sha256(
                f"formal-v7-equilibrium|fine|{split}|{row['stream_id']}".encode()
            ).hexdigest(),
        )
        return selected[:16], selected[16:]
    raise ValueError(phase)


def load_rows(selected):
    streams = []
    for row in selected:
        path = Path(row["path"])
        raw = np.fromfile(path, dtype="<f8")
        if raw.size != 100_001 * 5:
            raise RuntimeError(f"bad size: {path}")
        state = raw.reshape(100_001, 5)
        if not np.isfinite(state).all() or np.any(state[:, :3] <= 0):
            raise RuntimeError(f"invalid state: {path}")
        streams.append(state)
    return np.concatenate(streams), streams
