#!/usr/bin/env python3
"""Calibrate weak stationary-identity resolution on equilibrium streams."""

from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent.parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE/"stationary_density"))

from reduced_operator import WEAK_FUNCTION_NAMES, weak_generator_values  # noqa: E402

SPLITS = HERE/"DIRECT_STATE_SPLITS.csv"
OUT_CAL = HERE/"analysis"/"weak_identity_calibration.json"
OUT_RAW = HERE/"analysis"/"weak_identity_stream_means.csv"
OUT_TEST = HERE/"analysis"/"weak_identity_blind_tests.csv"
N_SNAP = 100_001
N_COL = 5
BOOTSTRAPS = 2000
BOOTSTRAP_SEED = 2026091402
DIAGNOSTIC_ONLY = {
    "sin_theta1", "cos_theta1", "sin_theta3", "cos_theta3",
    "I2_sin_theta1", "I2_sin_theta3",
}
GATE_INDICES = np.array([
    index for index, name in enumerate(WEAK_FUNCTION_NAMES)
    if name not in DIAGNOSTIC_ONLY
], dtype=int)


def load_state(path: Path) -> np.ndarray:
    return np.fromfile(path, dtype="<f8").reshape(N_SNAP, N_COL)


def stream_means(row: dict[str, str], chunk: int = 20_000) -> np.ndarray:
    state = load_state(REPO/row["relative_path"])
    total = np.zeros(len(WEAK_FUNCTION_NAMES), dtype=np.longdouble)
    count = 0
    t1 = 5.0 if row["case"].startswith("equilibrium") else 2.0
    t3 = 5.0 if row["case"].startswith("equilibrium") else 8.0
    for start in range(0, N_SNAP, chunk):
        values = weak_generator_values(state[start:start+chunk], t1, t3, 0.1)
        total += values.sum(axis=0, dtype=np.longdouble)
        count += len(values)
    return np.asarray(total/count, dtype=np.float64)


def mean_se(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return values.mean(axis=0), values.std(axis=0, ddof=1)/np.sqrt(values.shape[0])


def centered_max_z_quantile(values: np.ndarray, rng: np.random.Generator) -> tuple[float, np.ndarray]:
    centered = values-values.mean(axis=0, keepdims=True)
    maxima = np.empty(BOOTSTRAPS, dtype=np.float64)
    for index in range(BOOTSTRAPS):
        sample = centered[rng.integers(0, len(centered), size=len(centered))]
        mean, se = mean_se(sample)
        z = np.divide(mean, se, out=np.zeros_like(mean), where=se > 0)
        maxima[index] = np.max(np.abs(z))
    return float(np.quantile(maxima, 0.95)), maxima


def main() -> None:
    with SPLITS.open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle) if row["branch"] == "density"]
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["case"], row["split"])].append(row)

    matrices = {}
    raw_rows = []
    for key, members in sorted(grouped.items()):
        matrix = []
        for row in sorted(members, key=lambda item: int(item["stream_id"])):
            values = stream_means(row)
            matrix.append(values)
            for name, value in zip(WEAK_FUNCTION_NAMES, values):
                raw_rows.append({
                    "case": key[0], "split": key[1],
                    "stream_id": row["stream_id"], "function": name,
                    "mean_Lf": value,
                })
        matrices[key] = np.asarray(matrix)

    with OUT_RAW.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(raw_rows[0]))
        writer.writeheader(); writer.writerows(raw_rows)

    rng = np.random.default_rng(BOOTSTRAP_SEED)
    eq_coarse = matrices[("equilibrium_dt1e-3", "validation")]
    eq_fine = matrices[("equilibrium_dt2p5e-4", "validation")]
    coarse_mean, coarse_se = mean_se(eq_coarse)
    fine_mean, fine_se = mean_se(eq_fine)
    q_coarse, _ = centered_max_z_quantile(eq_coarse[:, GATE_INDICES], rng)
    q_fine, _ = centered_max_z_quantile(eq_fine[:, GATE_INDICES], rng)
    threshold = max(3.0, q_coarse, q_fine)
    floors = np.abs(coarse_mean-fine_mean)

    calibration = {
        "protocol_version": "section4-v1",
        "bootstrap_seed": BOOTSTRAP_SEED,
        "bootstrap_replicates": BOOTSTRAPS,
        "functions": list(WEAK_FUNCTION_NAMES),
        "gate_functions": [
            name for name in WEAK_FUNCTION_NAMES if name not in DIAGNOSTIC_ONLY
        ],
        "diagnostic_only_infinite_variance": sorted(DIAGNOSTIC_ONLY),
        "equilibrium_validation_max_z_q95": {
            "dt1e-3": q_coarse, "dt2p5e-4": q_fine,
        },
        "familywise_z_threshold": threshold,
        "per_function": [
            {
                "function": name,
                "coarse_mean": float(coarse_mean[i]), "coarse_se": float(coarse_se[i]),
                "fine_mean": float(fine_mean[i]), "fine_se": float(fine_se[i]),
                "coarse_fine_floor": float(floors[i]),
            }
            for i, name in enumerate(WEAK_FUNCTION_NAMES)
        ],
    }
    OUT_CAL.write_text(json.dumps(calibration, indent=2, sort_keys=True)+"\n")

    test_rows = []
    for case in ("equilibrium_dt1e-3", "equilibrium_dt2p5e-4",
                 "driven_dt1e-3", "driven_dt2p5e-4"):
        mean, se = mean_se(matrices[(case, "test")])
        for i, name in enumerate(WEAK_FUNCTION_NAMES):
            tolerance = max(threshold*se[i], floors[i])
            eligible = name not in DIAGNOSTIC_ONLY
            passed = bool(abs(mean[i]) <= tolerance) if eligible else None
            test_rows.append({
                "case": case, "function": name, "mean_Lf": mean[i], "se": se[i],
                "z": mean[i]/se[i] if se[i] > 0 else np.nan,
                "calibrated_floor": floors[i], "tolerance": tolerance,
                "gate_eligible": int(eligible),
                "status": "PASS" if passed else "FAIL" if passed is False else "DIAGNOSTIC_ONLY",
            })
    with OUT_TEST.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(test_rows[0]))
        writer.writeheader(); writer.writerows(test_rows)

    failed = [row for row in test_rows if row["status"] == "FAIL"]
    print(json.dumps({
        "familywise_z_threshold": threshold,
        "failures": failed,
        "all_blind_direct_trajectory_weak_identities_pass": not failed,
    }, indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
