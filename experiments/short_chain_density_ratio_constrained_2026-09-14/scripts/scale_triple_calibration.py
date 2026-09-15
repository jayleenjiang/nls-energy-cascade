#!/usr/bin/env python3
"""Apply the fixed R10 common shrink factor to a three-moment bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import tensorflow as tf

from triple_calibrated_transport_ratio_model import (
    TripleCalibratedTransportGibbsRatio,
    load_bundle,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--factor", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.factor != 0.625:
        raise ValueError("R10 freezes the inherited factor at 0.625")
    if args.output.exists() and any(args.output.iterdir()):
        raise SystemExit(f"refusing to overwrite {args.output}")
    config = json.loads((args.source / "model_config.json").read_text())
    source = load_bundle(args.source)
    alpha = args.factor * np.asarray(config["calibration_alpha"], np.float64)
    model = TripleCalibratedTransportGibbsRatio(
        config["architecture"], np.asarray(config["action_scales"], np.float64),
        int(config["seed"]), alpha,
        np.asarray(config["calibration_center"], np.float64),
        np.asarray(config["calibration_scale"], np.float64),
    )
    _ = model(tf.constant([[1.0, 1.0, 1.0, 0.0, 0.0]], tf.float64))
    model.set_weights(source.get_weights())
    source_hash = hashlib.sha256((args.source / "model_config.json").read_bytes()).hexdigest()
    model.save_bundle(args.output, {
        "training_base_seed": config.get("training_base_seed"),
        "source_r10_bundle": str(args.source),
        "source_r10_config_sha256": source_hash,
        "calibration_shrink_factor": args.factor,
        "unshrunk_calibration_alpha": config["calibration_alpha"],
        "calibration_observables": config["calibration_observables"],
        "method": "R10 fixed common-factor three-moment calibration shrinkage",
    })


if __name__ == "__main__":
    main()
