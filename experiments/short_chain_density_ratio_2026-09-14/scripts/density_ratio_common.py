#!/usr/bin/env python3
"""Frozen inputs and deterministic role assignments for density-ratio v1."""

from __future__ import annotations

import csv
import hashlib
import os
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent.parent
REPO = HERE.parents[1]
SOURCE = REPO / "experiments/short_chain_section4_2026-09-14"
CACHE = os.environ.get("NLS_SECTION4_RAW_CACHE")
SPLITS = SOURCE / "DIRECT_STATE_SPLITS.csv"
N_SNAP = 100_001
N_COL = 5
VERSION = "section4-density-ratio-v1"


def rows(case: str, split: str) -> list[dict[str, str]]:
    with SPLITS.open(newline="") as handle:
        out = [row for row in csv.DictReader(handle)
               if row["branch"] == "density" and row["case"] == case
               and row["split"] == split]
    return sorted(out, key=lambda row: int(row["stream_id"]))


def equilibrium_roles(case: str, split: str) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    ranked = sorted(
        rows(case, split),
        key=lambda row: hashlib.sha256(
            f"{VERSION}|equilibrium-role|{split}|{row['stream_id']}".encode()
        ).hexdigest(),
    )
    midpoint = len(ranked)//2
    return ranked[:midpoint], ranked[midpoint:]


def target_reference_rows(phase: str, timestep: str, split: str):
    suffix = "dt2p5e-4" if timestep == "fine" else "dt1e-3"
    equilibrium = f"equilibrium_{suffix}"
    if phase == "equilibrium":
        return equilibrium_roles(equilibrium, split)
    if phase == "driven":
        return rows(f"driven_{suffix}", split), rows(equilibrium, split)
    raise ValueError(phase)


def load_row(row: dict[str, str]) -> np.ndarray:
    path = REPO / row["relative_path"]
    if CACHE:
        relative = Path(row["relative_path"]).relative_to(
            "experiments/short_chain_section4_2026-09-14"
        )
        path = Path(CACHE)/relative
    raw = np.fromfile(path, dtype="<f8")
    expected = N_SNAP*N_COL
    if raw.size != expected:
        raise RuntimeError(f"{path}: {raw.size} values, expected {expected}")
    state = raw.reshape(N_SNAP, N_COL)
    if not np.isfinite(state).all() or np.any(state[:, :3] <= 0):
        raise RuntimeError(f"invalid state in {path}")
    return state


def load_rows(selected: list[dict[str, str]]) -> tuple[np.ndarray, list[np.ndarray]]:
    streams = [load_row(row) for row in selected]
    return np.concatenate(streams, axis=0), streams


def fixed_indices(length: int, count: int, salt: str) -> np.ndarray:
    seed = int(hashlib.sha256(f"{VERSION}|{salt}".encode()).hexdigest()[:16], 16)
    return np.random.default_rng(seed).choice(length, size=min(length, count), replace=False)
