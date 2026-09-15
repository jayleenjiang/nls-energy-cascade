#!/usr/bin/env python3
"""Fit the frozen R10 three-moment correction on coarse fitting rows."""

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

import train_formal_v3
from common_v2 import load_row
from transport_ratio_model import load_bundle as load_base
from triple_calibrated_transport_ratio_model import (
    TripleCalibratedTransportGibbsRatio,
)


def observables(state):
    i1, i2, i3, th1, th3 = state.T
    return np.column_stack((
        i1,
        i1 * i2 * np.sin(th1),
        i2 * i3 * np.sin(th3),
    ))


def log_ratio(model, state, batch=8192):
    return np.concatenate([
        model.log_ratio(tf.constant(state[start:start + batch], tf.float64)).numpy()
        for start in range(0, len(state), batch)
    ])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-model", type=Path, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise SystemExit(f"refusing to overwrite {args.output}")
    started = time.time()
    tf.keras.backend.set_floatx("float64")
    target_rows, reference_rows = train_formal_v3.combined_roles("driven", "coarse")
    target_sum = np.zeros(3)
    target_count = 0
    for row in target_rows:
        values = observables(load_row(row))
        target_sum += values.sum(axis=0)
        target_count += len(values)
    target_mean = target_sum / target_count
    base = load_base(args.base_model)
    gs = []
    fs = []
    for row in reference_rows:
        state = load_row(row)
        gs.append(observables(state))
        fs.append(log_ratio(base, state))
    gs = np.concatenate(gs)
    fs = np.concatenate(fs)
    center = gs.mean(axis=0)
    scale = gs.std(axis=0, ddof=1)
    standardized = (gs - center) / scale

    def weighted_mean(alpha):
        logw = fs + standardized @ alpha
        return np.exp(logw - logsumexp(logw)) @ gs

    solution = root(lambda alpha: weighted_mean(alpha) - target_mean, np.zeros(3))
    if not solution.success:
        raise RuntimeError(solution.message)
    alpha = np.asarray(solution.x, np.float64)
    fitted = weighted_mean(alpha)
    residual = fitted - target_mean
    if np.max(np.abs(residual)) > 1e-10:
        raise RuntimeError(f"calibration residual too large: {residual}")
    config = json.loads((args.base_model / "model_config.json").read_text())
    model = TripleCalibratedTransportGibbsRatio(
        config["architecture"], np.asarray(config["action_scales"], np.float64),
        args.seed, alpha, center, scale,
    )
    _ = model(tf.constant([[1.0, 1.0, 1.0, 0.0, 0.0]], tf.float64))
    model.set_weights(base.get_weights())
    base_hash = hashlib.sha256((args.base_model / "model.weights.h5").read_bytes()).hexdigest()
    model.save_bundle(args.output, {
        "training_base_seed": config.get("training_base_seed"),
        "base_model": str(args.base_model),
        "base_model_sha256": base_hash,
        "calibration_observables": [
            "I1", "I1*I2*sin(theta1)", "I2*I3*sin(theta3)",
        ],
        "calibration_target_mean": target_mean.tolist(),
        "calibration_fitted_mean": fitted.tolist(),
        "calibration_residual": residual.tolist(),
        "calibration_target_streams": [int(row["stream_id"]) for row in target_rows],
        "calibration_reference_streams": [int(row["stream_id"]) for row in reference_rows],
        "method": "R10 coarse train+validation three-moment exponential calibration",
        "timestep": "coarse",
    })
    result = {
        "seed": args.seed,
        "alpha": alpha.tolist(),
        "center": center.tolist(),
        "scale": scale.tolist(),
        "target_mean": target_mean.tolist(),
        "fitted_mean": fitted.tolist(),
        "residual": residual.tolist(),
        "target_stream_count": len(target_rows),
        "reference_stream_count": len(reference_rows),
        "runtime_seconds": time.time() - started,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "calibration.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
