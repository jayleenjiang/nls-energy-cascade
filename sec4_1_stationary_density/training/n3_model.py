#!/usr/bin/env python3
"""Boundary-regular periodic Gibbs density ratio for the three-mode chain (n = 3 port of n4_model.py,
feature version phase_regular_v4).  log rho(x) = -E(x)/5 + f(x) - log Z,  E = H/2."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import tensorflow as tf

from n3_common import energy_tf


class N3GibbsRatio(tf.keras.Model):
    def __init__(self, action_scales: np.ndarray, seed: int, width: int = 128, depth: int = 3,
                 anchor_temperature: float = 5.0):
        super().__init__(name="n3_gibbs_ratio")
        self.seed, self.width, self.depth = int(seed), int(width), int(depth)
        self.anchor_temperature = float(anchor_temperature)
        self.action_scales = tf.constant(np.asarray(action_scales, np.float64), tf.float64)
        self.hidden = [tf.keras.layers.Dense(
            width, activation="tanh", dtype=tf.float64,
            kernel_initializer=tf.keras.initializers.GlorotUniform(seed + 101 + j),
            bias_initializer="zeros", name=f"hidden_{j}") for j in range(self.depth)]
        self.output_layer = tf.keras.layers.Dense(1, dtype=tf.float64, kernel_initializer="zeros",
                                                  bias_initializer="zeros", name="log_ratio")

    def features(self, state: tf.Tensor) -> tf.Tensor:
        state = tf.cast(state, tf.float64)
        action = state[:, :3]
        h = action / (action + self.action_scales[None, :])
        h1, h2, h3 = tf.unstack(h, axis=-1)
        t2, t3 = tf.unstack(state[:, 3:], axis=-1)
        # theta_j = 2(phi_j - phi_{j-1}): a first harmonic is smooth in the Cartesian amplitudes only
        # when it vanishes linearly at both participating actions (as phase_regular_v4 for n = 4).
        g2, g3 = h1 * h2, h2 * h3
        values = [h1, h2, h3, h1 * h1, h2 * h2, h3 * h3, h1 * h2, h2 * h3, h1 * h3]
        for gate, angle in ((g2, t2), (g3, t3)):
            values += [gate * tf.sin(angle), gate * tf.cos(angle),
                       gate * gate * tf.sin(2.0 * angle), gate * gate * tf.cos(2.0 * angle)]
        # theta2 - theta3 = 2(2 phi2 - phi1 - phi3): order two in I2; theta2 + theta3 = 2(phi3 - phi1)
        minus_gate, plus_gate = h1 * h2 * h2 * h3, h1 * h3
        values += [minus_gate * tf.sin(t2 - t3), minus_gate * tf.cos(t2 - t3),
                   plus_gate * tf.sin(t2 + t3), plus_gate * tf.cos(t2 + t3)]
        regular = tf.stack(values, axis=-1)
        u = tf.math.log1p(action / self.action_scales[None, :])
        u1, u2, u3 = tf.unstack(u, axis=-1)
        action_features = tf.stack((u1, u2, u3, u1 * u1, u2 * u2, u3 * u3, u1 * u2, u2 * u3, u1 * u3), axis=-1)
        s1, s2, s3 = tf.unstack(self.action_scales)
        transport = tf.stack((
            action[:, 0] * action[:, 1] * tf.sin(t2) / (s1 * s2),
            action[:, 0] * action[:, 1] * tf.cos(t2) / (s1 * s2),
            action[:, 1] * action[:, 2] * tf.sin(t3) / (s2 * s3),
            action[:, 1] * action[:, 2] * tf.cos(t3) / (s2 * s3)), axis=-1)
        return tf.concat((regular, action_features, transport), axis=-1)

    def log_ratio(self, state: tf.Tensor) -> tf.Tensor:
        value = self.features(state)
        for layer in self.hidden:
            value = layer(value)
        return self.output_layer(value)[:, 0]

    def log_prob(self, state: tf.Tensor) -> tf.Tensor:
        return -energy_tf(state) / tf.constant(self.anchor_temperature, tf.float64) + self.log_ratio(state)

    def call(self, state: tf.Tensor) -> tf.Tensor:
        return self.log_ratio(state)

    def save_bundle(self, directory: Path, extra: dict) -> None:
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        _ = self(tf.constant([[1., 1., 1., 0., 0.]], tf.float64))
        config = {"action_scales": self.action_scales.numpy().tolist(), "seed": self.seed, "width": self.width,
                  "depth": self.depth, "anchor_temperature": self.anchor_temperature,
                  "feature_version": "phase_regular_v4_n3", **extra}
        (directory / "model_config.json").write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")
        self.save_weights(directory / "model.weights.h5")


def load_bundle(directory: Path) -> N3GibbsRatio:
    directory = Path(directory)
    config = json.loads((directory / "model_config.json").read_text())
    model = N3GibbsRatio(np.asarray(config["action_scales"]), int(config["seed"]), int(config["width"]),
                         int(config["depth"]), float(config.get("anchor_temperature", 5.0)))
    _ = model(tf.constant([[1., 1., 1., 0., 0.]], tf.float64))
    model.load_weights(directory / "model.weights.h5")
    return model
