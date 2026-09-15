#!/usr/bin/env python3
"""Wrap frozen zero equilibrium models in the R11 three-moment model class."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import tensorflow as tf

from transport_ratio_model import load_bundle as load_base
from triple_calibrated_transport_ratio_model import TripleCalibratedTransportGibbsRatio


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-model", type=Path, action="append", required=True)
    parser.add_argument("--seed", type=int, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if len(args.base_model) != 3 or len(args.seed) != 3:
        raise ValueError("exactly three base models and three seeds are required")
    if args.output.exists() and any(args.output.iterdir()):
        raise SystemExit(f"refusing to overwrite {args.output}")
    tf.keras.backend.set_floatx("float64")
    manifest = []
    for base_path, seed in zip(args.base_model, args.seed):
        base = load_base(base_path)
        config = json.loads((base_path / "model_config.json").read_text())
        model = TripleCalibratedTransportGibbsRatio(
            config["architecture"], np.asarray(config["action_scales"], np.float64),
            seed, np.zeros(3), np.zeros(3), np.ones(3),
        )
        _ = model(tf.constant([[1.0, 1.0, 1.0, 0.0, 0.0]], tf.float64))
        model.set_weights(base.get_weights())
        directory = args.output / f"seed_{seed}"
        base_hash = hashlib.sha256((base_path / "model.weights.h5").read_bytes()).hexdigest()
        model.save_bundle(directory, {
            "training_base_seed": config.get("training_base_seed"),
            "base_zero_model": str(base_path),
            "base_zero_model_sha256": base_hash,
            "calibration_observables": [
                "I1", "I1*I2*sin(theta1)", "I2*I3*sin(theta3)",
            ],
            "calibration_shrink_factor": 0.625,
            "method": "R11 exact-equilibrium zero ratio in three-moment model class",
            "phase": "equilibrium",
            "timestep": "fine",
        })
        manifest.append({"seed": seed, "path": str(directory), "base_sha256": base_hash})
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "manifest.json").write_text(
        json.dumps({"models": manifest}, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()

