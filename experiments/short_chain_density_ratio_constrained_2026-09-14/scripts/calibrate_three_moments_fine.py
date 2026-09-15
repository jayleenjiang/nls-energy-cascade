#!/usr/bin/env python3
"""Fit the frozen R11 three-moment correction on the original fine calibration rows."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import tensorflow as tf
from scipy.optimize import root
from scipy.special import logsumexp

from common_v2 import load_row as load_old_row, rows as old_rows
from transport_ratio_model import load_bundle as load_base
from triple_calibrated_transport_ratio_model import (
    TripleCalibratedTransportGibbsRatio,
)
from v4_common import load_rows as load_v4_rows, rows as v4_rows


def observables(state: np.ndarray) -> np.ndarray:
    i1, i2, i3, th1, th3 = state.T
    return np.column_stack((
        i1,
        i1 * i2 * np.sin(th1),
        i2 * i3 * np.sin(th3),
    ))


def log_ratio(model, state: np.ndarray, batch: int = 8192) -> np.ndarray:
    return np.concatenate([
        model.log_ratio(tf.constant(state[start:start + batch], tf.float64)).numpy()
        for start in range(0, len(state), batch)
    ])


def iter_opened(case: str):
    old_case = f"{case}_dt2p5e-4"
    for split in ("train", "validation", "test"):
        for row in old_rows(old_case, split):
            yield f"old:{split}:{row['stream_id']}", load_old_row(row)
    for role in ("validation", "test"):
        selected = v4_rows(f"{case}_fine", role)
        _, streams = load_v4_rows(selected)
        for row, state in zip(selected, streams):
            yield f"v4:{role}:{row['stream_id']}", state


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-model", type=Path, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise SystemExit(f"refusing to overwrite {args.output}")

    started = time.time()
    tf.keras.backend.set_floatx("float64")
    target_sum = np.zeros(3)
    target_count = 0
    target_streams = []
    for label, state in iter_opened("driven"):
        values = observables(state)
        target_sum += values.sum(axis=0)
        target_count += len(values)
        target_streams.append(label)
    target_mean = target_sum / target_count

    base = load_base(args.base_model)
    reference_observables = []
    reference_log_ratio = []
    reference_streams = []
    for label, state in iter_opened("equilibrium"):
        reference_observables.append(observables(state))
        reference_log_ratio.append(log_ratio(base, state))
        reference_streams.append(label)
    reference_observables = np.concatenate(reference_observables)
    reference_log_ratio = np.concatenate(reference_log_ratio)
    center = reference_observables.mean(axis=0)
    scale = reference_observables.std(axis=0, ddof=1)
    standardized = (reference_observables - center) / scale

    def weighted_mean(alpha: np.ndarray) -> np.ndarray:
        logw = reference_log_ratio + standardized @ alpha
        weights = np.exp(logw - logsumexp(logw))
        return weights @ reference_observables

    solution = root(lambda alpha: weighted_mean(alpha) - target_mean, np.zeros(3))
    if not solution.success:
        raise RuntimeError(f"three-moment calibration failed: {solution.message}")
    alpha = np.asarray(solution.x, np.float64)
    fitted_mean = weighted_mean(alpha)
    residual = fitted_mean - target_mean
    if np.max(np.abs(residual)) > 1e-10:
        raise RuntimeError(f"calibration residual too large: {residual}")

    base_config = json.loads((args.base_model / "model_config.json").read_text())
    model = TripleCalibratedTransportGibbsRatio(
        base_config["architecture"],
        np.asarray(base_config["action_scales"], np.float64),
        args.seed, alpha, center, scale,
    )
    _ = model(tf.constant([[1.0, 1.0, 1.0, 0.0, 0.0]], tf.float64))
    model.set_weights(base.get_weights())
    base_hash = hashlib.sha256((args.base_model / "model.weights.h5").read_bytes()).hexdigest()
    model.save_bundle(args.output, {
        "training_base_seed": base_config.get("training_base_seed"),
        "base_model": str(args.base_model),
        "base_model_sha256": base_hash,
        "calibration_observables": [
            "I1", "I1*I2*sin(theta1)", "I2*I3*sin(theta3)",
        ],
        "calibration_target_mean": target_mean.tolist(),
        "calibration_fitted_mean": fitted_mean.tolist(),
        "calibration_residual": residual.tolist(),
        "calibration_target_streams": target_streams,
        "calibration_reference_streams": reference_streams,
        "method": "R11 fine opened-development-only three-moment calibration",
        "timestep": "fine",
    })
    result = {
        "seed": args.seed,
        "base_model": str(args.base_model),
        "base_model_sha256": base_hash,
        "alpha": alpha.tolist(),
        "center": center.tolist(),
        "scale": scale.tolist(),
        "target_mean": target_mean.tolist(),
        "fitted_mean": fitted_mean.tolist(),
        "residual": residual.tolist(),
        "target_stream_count": len(target_streams),
        "reference_stream_count": len(reference_streams),
        "runtime_seconds": time.time() - started,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "calibration.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

