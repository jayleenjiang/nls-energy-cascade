#!/usr/bin/env python3
"""Shared frozen Section-4.3 EDMD utilities for direct float64 states."""

from __future__ import annotations

import csv
import importlib.util
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent.parent
REPO = ROOT.parents[1]
EXTENDED_SOURCE = (
    REPO / "experiments" / "edmd_oscillatory_extension_2026-09-07"
    / "run_extended_edmd.py"
)
SPEC = importlib.util.spec_from_file_location("section4_existing_extended_edmd", EXTENDED_SOURCE)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot import {EXTENDED_SOURCE}")
extended = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = extended
SPEC.loader.exec_module(extended)
base = extended.base

NROWS = 100_001
NCOLS = 5
SNAPSHOT = 0.01
N_ORIGINS = 2048
LAGS = (0.02, 0.05, 0.10, 0.25, 0.50, 1.00)
LAG_STEPS = tuple(int(round(t/SNAPSHOT)) for t in LAGS)
DICTIONARIES = ("E1", "E2", "E3")
CUTOFFS = (1.0e-8, 1.0e-10, 1.0e-12)
CASES = (
    "driven_dt1e-3", "driven_dt2p5e-4",
    "equilibrium_dt1e-3", "equilibrium_dt2p5e-4",
)
SPLITS = ROOT/"DIRECT_STATE_SPLITS.csv"


@dataclass
class Aggregate:
    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    pairs: int
    streams: int


def split_paths(case: str, split: str) -> list[Path]:
    with SPLITS.open(newline="") as handle:
        rows = [r for r in csv.DictReader(handle)
                if r["branch"] == "modes" and r["case"] == case and r["split"] == split]
    rows.sort(key=lambda r: int(r["stream_id"]))
    return [REPO/r["relative_path"] for r in rows]


def load_state(path: Path) -> np.ndarray:
    return np.fromfile(path, dtype="<f8").reshape(NROWS, NCOLS)


def state_to_raw12(state: np.ndarray) -> np.ndarray:
    i1, i2, i3, theta1, theta3 = state.T
    result = np.zeros((len(state), 12), dtype=np.float64)
    result[:, 0] = np.cos(theta1)
    result[:, 1] = np.sin(theta1)
    result[:, 2] = np.cos(theta3)
    result[:, 3] = np.sin(theta3)
    result[:, 4] = np.cos(theta1-theta3)
    result[:, 5] = i2
    result[:, 6] = i1+i3
    result[:, 7] = i1-i3
    result[:, 8] = result[:, 0]+result[:, 2]
    result[:, 9] = result[:, 0]-result[:, 2]
    result[:, 10] = result[:, 1]+result[:, 3]
    result[:, 11] = result[:, 1]-result[:, 3]
    return result


def feature_matrix(state: np.ndarray, log_mean: np.ndarray,
                   log_sd: np.ndarray) -> np.ndarray:
    return extended.extended_features(state_to_raw12(state), log_mean, log_sd)


def train_log_stats(paths: list[Path]) -> tuple[np.ndarray, np.ndarray]:
    count = 0
    total = np.zeros(3, dtype=np.longdouble)
    squares = np.zeros(3, dtype=np.longdouble)
    for path in paths:
        state = load_state(path)
        if not np.isfinite(state).all() or np.any(state[:, :3] <= 0):
            raise RuntimeError(f"invalid direct state {path}")
        values = np.log(state[:, :3]).astype(np.longdouble)
        total += values.sum(axis=0)
        squares += (values*values).sum(axis=0)
        count += len(values)
    mean = total/count
    variance = (squares-count*mean*mean)/(count-1)
    return np.asarray(mean, dtype=np.float64), np.sqrt(np.asarray(variance, dtype=np.float64))


def origin_indices() -> np.ndarray:
    maximum = max(LAG_STEPS)
    values = np.rint(np.linspace(0, NROWS-1-maximum, N_ORIGINS)).astype(np.int64)
    if len(np.unique(values)) != N_ORIGINS:
        raise RuntimeError("origin quadrature contains duplicates")
    return values


def aggregate(paths: list[Path], log_mean: np.ndarray, log_sd: np.ndarray,
              label: str) -> Aggregate:
    origins = origin_indices()
    k = len(extended.DICT_INDICES["E3"])
    A = np.zeros((k, k), dtype=np.float64)
    B = np.zeros((len(LAGS), k, k), dtype=np.float64)
    C = np.zeros((len(LAGS), k, k), dtype=np.float64)
    for number, path in enumerate(paths, start=1):
        state = load_state(path)
        X = feature_matrix(state[origins], log_mean, log_sd)
        A += X.T@X
        for ell, step in enumerate(LAG_STEPS):
            Y = feature_matrix(state[origins+step], log_mean, log_sd)
            B[ell] += X.T@Y
            C[ell] += Y.T@Y
        print(f"{label}: stream {number}/{len(paths)}", flush=True)
    pairs = len(paths)*N_ORIGINS
    return Aggregate(A/pairs, B/pairs, C/pairs, pairs, len(paths))


def subset(matrix: np.ndarray, indices: np.ndarray) -> np.ndarray:
    return matrix[np.ix_(indices, indices)]


def fit(agg: Aggregate, dictionary: str, tau: float, cutoff: float,
        vectors: bool = True):
    indices = extended.DICT_INDICES[dictionary]
    ell = LAGS.index(tau)
    return base.solve_edmd(
        subset(agg.A, indices), subset(agg.B[ell], indices),
        tau, cutoff, vectors=vectors,
    )


def select_complex(result):
    candidates = [item for item in base.canonical_modes(result)
                  if 4.5 <= item[2].imag <= 6.0 and -20 <= item[2].real <= -0.05]
    return max(candidates, key=lambda item: item[2].real) if candidates else None


def select_real(result):
    candidates = [item for item in base.canonical_modes(result)
                  if -3 <= item[2].real <= -0.05 and abs(item[2].imag) <= 0.10
                  and 0 < abs(item[1]) <= 1.05]
    return max(candidates, key=lambda item: item[2].real) if candidates else None


def nearest(reference: complex, result):
    candidates = base.canonical_modes(result)
    if not candidates:
        return None
    chosen = min(candidates, key=lambda item: abs(item[2]-reference))
    tolerance = max(0.50, 0.35*abs(reference))
    return chosen if abs(chosen[2]-reference) <= tolerance else None


def stability(values: list[complex]) -> tuple[bool, dict]:
    if not values:
        return False, {"count": 0}
    reals = np.asarray([v.real for v in values])
    imags = np.asarray([abs(v.imag) for v in values])
    real_range = float(reals.max()-reals.min())
    imag_range = float(imags.max()-imags.min())
    real_tolerance = max(0.15, 0.20*abs(float(np.median(reals))))
    imag_tolerance = max(0.50, 0.20*float(np.median(imags)))
    passed = real_range <= real_tolerance and imag_range <= imag_tolerance
    return passed, {
        "count": len(values), "values": [[v.real, v.imag] for v in values],
        "real_range": real_range, "real_tolerance": real_tolerance,
        "imag_range": imag_range, "imag_tolerance": imag_tolerance,
    }


def validation_residual(agg: Aggregate, dictionary: str, tau: float,
                        mu: complex, coefficient: np.ndarray) -> float:
    indices = extended.DICT_INDICES[dictionary]
    ell = LAGS.index(tau)
    A = subset(agg.A, indices)
    B = subset(agg.B[ell], indices)
    C = subset(agg.C[ell], indices)
    denominator = float(np.real(np.vdot(coefficient, A@coefficient)))
    cross = np.vdot(coefficient, B@coefficient)
    numerator = float(np.real(
        np.vdot(coefficient, C@coefficient)
        - 2*np.conj(mu)*cross
        + abs(mu)**2*np.vdot(coefficient, A@coefficient)
    ))
    return math.sqrt(max(numerator, 0.0)/denominator)

