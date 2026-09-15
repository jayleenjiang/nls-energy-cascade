#!/usr/bin/env python3
"""Compare R11 fine and R10 coarse three-moment ratios on both test supports."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import tensorflow as tf

from triple_calibrated_transport_ratio_model import load_bundle
from v4_common import evaluation_roles as v4_roles, load_rows as v4_load
from v6_common import evaluation_roles as v6_roles, load_rows as v6_load


def indices(length: int, count: int, salt: str) -> np.ndarray:
    seed = int(hashlib.sha256(salt.encode()).hexdigest()[:16], 16)
    return np.random.default_rng(seed).choice(length, min(length, count), replace=False)


def values(model, state: np.ndarray, batch: int = 8192) -> np.ndarray:
    return np.concatenate([
        model.log_ratio(tf.constant(state[start:start + batch], tf.float64)).numpy()
        for start in range(0, len(state), batch)
    ])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fine-model", type=Path, action="append", required=True)
    parser.add_argument("--coarse-model", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if len(args.fine_model) != 3 or len(args.coarse_model) != 3:
        raise ValueError("exactly three models per timestep are required")
    if args.output.exists():
        raise SystemExit(f"refusing to overwrite {args.output}")
    tf.keras.backend.set_floatx("float64")
    fine = [load_bundle(path) for path in args.fine_model]
    coarse = [load_bundle(path) for path in args.coarse_model]
    rows = []
    ensemble = {}
    for support in ("fine", "coarse"):
        if support == "fine":
            target_rows, _ = v6_roles("driven", "fine", "test")
            state, _ = v6_load(target_rows)
        else:
            target_rows, _ = v4_roles("driven", "coarse", "test")
            state, _ = v4_load(target_rows)
        state = state[indices(len(state), 50_000, f"r11-three-moment|{support}")]
        fv = np.asarray([values(model, state) for model in fine])
        cv = np.asarray([values(model, state) for model in coarse])
        fv -= fv.mean(axis=1, keepdims=True)
        cv -= cv.mean(axis=1, keepdims=True)
        for index in range(3):
            rms = float(np.sqrt(np.mean((fv[index] - cv[index]) ** 2)))
            rows.append({
                "support": support,
                "comparison": f"paired_seed_{index + 1}",
                "centered_log_ratio_rms": rms,
                "gate": 0.10,
                "pass": int(rms <= 0.10),
            })
        rms = float(np.sqrt(np.mean((fv.mean(axis=0) - cv.mean(axis=0)) ** 2)))
        ensemble[support] = rms
        rows.append({
            "support": support,
            "comparison": "ensemble_mean",
            "centered_log_ratio_rms": rms,
            "gate": 0.10,
            "pass": int(rms <= 0.10),
        })
    args.output.mkdir(parents=True)
    with (args.output / "timestep_comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    result = {
        "sample_count_per_support": 50_000,
        "ensemble_rms_by_support": ensemble,
        "maximum_ensemble_rms": max(ensemble.values()),
        "gate": 0.10,
        "complete_pass": bool(max(ensemble.values()) <= 0.10),
    }
    (args.output / "verdict.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

