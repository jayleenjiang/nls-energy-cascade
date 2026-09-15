#!/usr/bin/env python3
"""Transport density ratio with a frozen three-moment correction."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import tensorflow as tf

from transport_ratio_model import TransportGibbsRatio


class TripleCalibratedTransportGibbsRatio(TransportGibbsRatio):
    def __init__(
        self, architecture: str, action_scales: np.ndarray, seed: int,
        calibration_alpha: np.ndarray, calibration_center: np.ndarray,
        calibration_scale: np.ndarray,
    ):
        super().__init__(architecture, action_scales, seed)
        self.calibration_alpha = tf.constant(calibration_alpha, tf.float64)
        self.calibration_center = tf.constant(calibration_center, tf.float64)
        self.calibration_scale = tf.constant(calibration_scale, tf.float64)

    def calibration_observables(self, state: tf.Tensor) -> tf.Tensor:
        state = tf.cast(state, tf.float64)
        i1, i2, i3, th1, th3 = tf.unstack(state, axis=-1)
        return tf.stack((
            i1,
            i1 * i2 * tf.sin(th1),
            i2 * i3 * tf.sin(th3),
        ), axis=-1)

    def log_ratio(self, state: tf.Tensor) -> tf.Tensor:
        base = super().log_ratio(state)
        g = self.calibration_observables(state)
        correction = tf.reduce_sum(
            self.calibration_alpha[None, :]
            * (g - self.calibration_center[None, :])
            / self.calibration_scale[None, :],
            axis=-1,
        )
        return base + correction

    def save_bundle(self, directory: Path, extra: dict) -> None:
        directory.mkdir(parents=True, exist_ok=True)
        _ = self(tf.constant([[1.0, 1.0, 1.0, 0.0, 0.0]], tf.float64))
        config = {
            "model_type": "transport_hybrid_three_moment_calibration",
            "architecture": self.architecture,
            "action_scales": self.action_scales.numpy().tolist(),
            "seed": self.seed,
            "calibration_alpha": self.calibration_alpha.numpy().tolist(),
            "calibration_center": self.calibration_center.numpy().tolist(),
            "calibration_scale": self.calibration_scale.numpy().tolist(),
            **extra,
        }
        (directory / "model_config.json").write_text(
            json.dumps(config, indent=2, sort_keys=True) + "\n"
        )
        self.save_weights(directory / "model.weights.h5")


def load_bundle(directory: Path) -> TripleCalibratedTransportGibbsRatio:
    directory = Path(directory)
    config = json.loads((directory / "model_config.json").read_text())
    model = TripleCalibratedTransportGibbsRatio(
        config["architecture"], np.asarray(config["action_scales"], np.float64),
        int(config["seed"]), np.asarray(config["calibration_alpha"], np.float64),
        np.asarray(config["calibration_center"], np.float64),
        np.asarray(config["calibration_scale"], np.float64),
    )
    _ = model(tf.constant([[1.0, 1.0, 1.0, 0.0, 0.0]], tf.float64))
    model.load_weights(directory / "model.weights.h5")
    return model
