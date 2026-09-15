#!/usr/bin/env python3
"""Hybrid density ratio with explicit smooth bond-transport features."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import tensorflow as tf

from hybrid_ratio_model import HybridGibbsRatio
from regular_ratio_model import ARCHITECTURES


class TransportGibbsRatio(HybridGibbsRatio):
    def features(self, state: tf.Tensor) -> tf.Tensor:
        hybrid = super().features(state)
        state = tf.cast(state, tf.float64)
        i1, i2, i3, th1, th3 = tf.unstack(state, axis=-1)
        s1, s2, s3 = tf.unstack(self.action_scales)
        transport = tf.stack((
            i1*i2*tf.sin(th1)/(s1*s2),
            i1*i2*tf.cos(th1)/(s1*s2),
            i2*i3*tf.sin(th3)/(s2*s3),
            i2*i3*tf.cos(th3)/(s2*s3),
        ), axis=-1)
        return tf.concat((hybrid, transport), axis=-1)

    def save_bundle(self, directory: Path, extra: dict) -> None:
        directory.mkdir(parents=True, exist_ok=True)
        _ = self(tf.constant([[1.0, 1.0, 1.0, 0.0, 0.0]], tf.float64))
        config = {
            "model_type": "transport_hybrid_four_cross_terms",
            "architecture": self.architecture,
            "action_scales": self.action_scales.numpy().tolist(),
            "seed": self.seed,
            **extra,
        }
        (directory/"model_config.json").write_text(
            json.dumps(config, indent=2, sort_keys=True)+"\n"
        )
        self.save_weights(directory/"model.weights.h5")


def load_bundle(directory: Path) -> TransportGibbsRatio:
    directory = Path(directory)
    config = json.loads((directory/"model_config.json").read_text())
    model = TransportGibbsRatio(
        config["architecture"], np.asarray(config["action_scales"], np.float64),
        int(config["seed"]),
    )
    _ = model(tf.constant([[1.0, 1.0, 1.0, 0.0, 0.0]], tf.float64))
    model.load_weights(directory/"model.weights.h5")
    return model


def transplant_hybrid_bundle(source: Path, destination: Path, seed: int) -> None:
    from hybrid_ratio_model import load_bundle as load_hybrid

    source = Path(source)
    old = load_hybrid(source)
    config = json.loads((source/"model_config.json").read_text())
    new = TransportGibbsRatio(
        config["architecture"], np.asarray(config["action_scales"], np.float64), seed
    )
    dummy = tf.constant([[1.0, 1.0, 1.0, 0.0, 0.0]], tf.float64)
    _ = new(dummy)
    old_weights = old.get_weights()
    new_weights = new.get_weights()
    if old_weights[0].shape[0]+4 != new_weights[0].shape[0]:
        raise RuntimeError("unexpected transport feature dimension")
    new_weights[0][:old_weights[0].shape[0], :] = old_weights[0]
    new_weights[0][old_weights[0].shape[0]:, :] = 0.0
    for index in range(1, len(old_weights)):
        if old_weights[index].shape != new_weights[index].shape:
            raise RuntimeError(f"weight shape mismatch at {index}")
        new_weights[index] = old_weights[index]
    new.set_weights(new_weights)
    probe = tf.constant([
        [0.2, 0.7, 1.1, -1.0, 2.0],
        [3.0, 2.0, 4.0, 0.2, -0.4],
    ], tf.float64)
    error = float(tf.reduce_max(tf.abs(
        new.log_ratio(probe)-old.log_ratio(probe)
    )).numpy())
    if error > 1e-14:
        raise RuntimeError(f"transport transplant changed function: {error}")
    new.save_bundle(destination, {
        "transplanted_from": str(source),
        "transplant_max_abs_log_ratio_error": error,
        "schedule": "initial_transport_transplant", "best_epoch": 0,
        "development_logz": config["development_logz"],
    })
