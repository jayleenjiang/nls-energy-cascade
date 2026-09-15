#!/usr/bin/env python3
"""Boundary-regular ratio with non-saturating action-only corrections."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import tensorflow as tf

from regular_ratio_model import ARCHITECTURES, RegularGibbsRatio, energy_tf


class HybridGibbsRatio(RegularGibbsRatio):
    """Append log1p action terms while retaining regular angular gates."""

    def features(self, state: tf.Tensor) -> tf.Tensor:
        regular = super().features(state)
        action = tf.cast(state[:, :3], tf.float64)
        u = tf.math.log1p(action/self.action_scales[None, :])
        u1, u2, u3 = tf.unstack(u, axis=-1)
        action_correction = tf.stack((
            u1, u2, u3, u1*u1, u2*u2, u3*u3,
            u1*u2, u2*u3, u1*u3,
        ), axis=-1)
        return tf.concat((regular, action_correction), axis=-1)

    def log_prob(self, state: tf.Tensor) -> tf.Tensor:
        return -energy_tf(state)/tf.constant(5.0, tf.float64)+self.log_ratio(state)

    def save_bundle(self, directory: Path, extra: dict) -> None:
        directory.mkdir(parents=True, exist_ok=True)
        _ = self(tf.constant([[1.0, 1.0, 1.0, 0.0, 0.0]], tf.float64))
        config = {
            "model_type": "hybrid_regular_log1p_action",
            "architecture": self.architecture,
            "action_scales": self.action_scales.numpy().tolist(),
            "seed": self.seed,
            **extra,
        }
        (directory/"model_config.json").write_text(
            json.dumps(config, indent=2, sort_keys=True)+"\n"
        )
        self.save_weights(directory/"model.weights.h5")


def load_bundle(directory: Path) -> HybridGibbsRatio:
    directory = Path(directory)
    config = json.loads((directory/"model_config.json").read_text())
    model = HybridGibbsRatio(
        config["architecture"], np.asarray(config["action_scales"], np.float64),
        int(config["seed"]),
    )
    _ = model(tf.constant([[1.0, 1.0, 1.0, 0.0, 0.0]], tf.float64))
    model.load_weights(directory/"model.weights.h5")
    return model


def transplant_regular_bundle(source: Path, destination: Path, seed: int) -> None:
    source = Path(source); destination = Path(destination)
    config = json.loads((source/"model_config.json").read_text())
    regular = RegularGibbsRatio(
        config["architecture"], np.asarray(config["action_scales"], np.float64),
        int(config["seed"]),
    )
    dummy = tf.constant([[1.0, 1.0, 1.0, 0.0, 0.0]], tf.float64)
    _ = regular(dummy)
    regular.load_weights(source/"model.weights.h5")
    hybrid = HybridGibbsRatio(
        config["architecture"], np.asarray(config["action_scales"], np.float64), seed,
    )
    _ = hybrid(dummy)
    old_weights = regular.get_weights()
    new_weights = hybrid.get_weights()
    if old_weights[0].shape[0]+9 != new_weights[0].shape[0]:
        raise RuntimeError("unexpected first-layer dimensions")
    new_weights[0][:old_weights[0].shape[0], :] = old_weights[0]
    new_weights[0][old_weights[0].shape[0]:, :] = 0.0
    for index in range(1, len(old_weights)):
        if old_weights[index].shape != new_weights[index].shape:
            raise RuntimeError(f"weight shape mismatch at {index}")
        new_weights[index] = old_weights[index]
    hybrid.set_weights(new_weights)
    probe = tf.constant([
        [0.2, 0.7, 1.1, -1.0, 2.0],
        [3.0, 2.0, 4.0, 0.2, -0.4],
    ], tf.float64)
    maximum_error = float(tf.reduce_max(tf.abs(
        hybrid.log_ratio(probe)-regular.log_ratio(probe)
    )).numpy())
    if maximum_error > 1e-14:
        raise RuntimeError(f"transplant changed represented function: {maximum_error}")
    hybrid.save_bundle(destination, {
        "transplanted_from": str(source),
        "transplant_max_abs_log_ratio_error": maximum_error,
        "schedule": "initial_transplant",
        "best_epoch": 0,
        "development_logz": config["development_logz"],
    })
