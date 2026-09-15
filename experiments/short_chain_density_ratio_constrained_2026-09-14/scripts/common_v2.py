#!/usr/bin/env python3
"""Frozen input loading and train-only development roles for density v2."""

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
VERSION = "section4-density-ratio-constrained-development-v1"


def rows(case: str, split: str) -> list[dict[str, str]]:
    with SPLITS.open(newline="") as handle:
        selected = [
            row for row in csv.DictReader(handle)
            if row["branch"] == "density" and row["case"] == case
            and row["split"] == split
        ]
    return sorted(selected, key=lambda row: int(row["stream_id"]))


def development_roles(case: str) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    ranked = sorted(
        rows(case, "train"),
        key=lambda row: hashlib.sha256(
            f"{VERSION}|fit-development|{case}|{row['stream_id']}".encode()
        ).hexdigest(),
    )
    if len(ranked) != 32:
        raise RuntimeError(f"expected 32 training streams for {case}, got {len(ranked)}")
    return ranked[:24], ranked[24:]


def driven_reference_roles(role: str):
    driven = "driven_dt2p5e-4"
    equilibrium = "equilibrium_dt2p5e-4"
    index = 0 if role == "fit" else 1
    if role not in ("fit", "development"):
        raise ValueError(role)
    return development_roles(driven)[index], development_roles(equilibrium)[index]


def formal_training_roles(phase: str, timestep: str):
    suffix = "dt2p5e-4" if timestep == "fine" else "dt1e-3"
    equilibrium = f"equilibrium_{suffix}"
    if phase == "driven":
        return rows(f"driven_{suffix}", "train"), rows(equilibrium, "train")
    if phase == "equilibrium":
        ranked = sorted(
            rows(equilibrium, "train"),
            key=lambda row: hashlib.sha256(
                f"{VERSION}|formal-equilibrium|{timestep}|{row['stream_id']}".encode()
            ).hexdigest(),
        )
        return ranked[:16], ranked[16:]
    raise ValueError(phase)


def formal_evaluation_roles(phase: str, timestep: str, split: str):
    if split not in ("validation", "test"):
        raise ValueError(split)
    suffix = "dt2p5e-4" if timestep == "fine" else "dt1e-3"
    equilibrium = f"equilibrium_{suffix}"
    if phase == "driven":
        return rows(f"driven_{suffix}", split), rows(equilibrium, split)
    if phase == "equilibrium":
        ranked = sorted(
            rows(equilibrium, split),
            key=lambda row: hashlib.sha256(
                f"{VERSION}|formal-equilibrium|{timestep}|{split}|{row['stream_id']}".encode()
            ).hexdigest(),
        )
        midpoint = len(ranked)//2
        return ranked[:midpoint], ranked[midpoint:]
    raise ValueError(phase)


def load_row(row: dict[str, str]) -> np.ndarray:
    path = REPO / row["relative_path"]
    if CACHE:
        relative = Path(row["relative_path"]).relative_to(
            "experiments/short_chain_section4_2026-09-14"
        )
        path = Path(CACHE) / relative
    raw = np.fromfile(path, dtype="<f8")
    if raw.size != N_SNAP * N_COL:
        raise RuntimeError(f"{path}: {raw.size} values, expected {N_SNAP*N_COL}")
    state = raw.reshape(N_SNAP, N_COL)
    if not np.isfinite(state).all() or np.any(state[:, :3] <= 0):
        raise RuntimeError(f"invalid state in {path}")
    return state


def load_rows(selected: list[dict[str, str]]) -> tuple[np.ndarray, list[np.ndarray]]:
    streams = [load_row(row) for row in selected]
    return np.concatenate(streams), streams


def fixed_indices(length: int, count: int, salt: str) -> np.ndarray:
    seed = int(hashlib.sha256(f"{VERSION}|{salt}".encode()).hexdigest()[:16], 16)
    return np.random.default_rng(seed).choice(length, min(length, count), replace=False)
