#!/usr/bin/env python3
"""Gibbs-anchored ratio with boundary-regular periodic coordinates."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import tensorflow as tf


ARCHITECTURES = {"regular64": (64, 64), "regular128": (128, 128)}


def energy_tf(state: tf.Tensor) -> tf.Tensor:
    i1, i2, i3, th1, th3 = tf.unstack(tf.cast(state, tf.float64), axis=-1)
    mass = i1 + i2 + i3
    return (
        0.5 * mass * mass - 0.25 * (i1 * i1 + i2 * i2 + i3 * i3)
        + i1 * i2 * tf.cos(th1) + i2 * i3 * tf.cos(th3)
    )


class RegularGibbsRatio(tf.keras.Model):
    """A smooth log ratio with angular dependence gated at singular actions.

    theta1 is undefined at I1=0 and theta3 is undefined at I3=0.  The feature
    map therefore multiplies every theta1 harmonic by h1 and every theta3
    harmonic by h3; joint harmonics contain h1*h3.  No bare angular coordinate
    is supplied to the network.
    """

    def __init__(self, architecture: str, action_scales: np.ndarray, seed: int):
        super().__init__(name="regular_gibbs_ratio")
        if architecture not in ARCHITECTURES:
            raise ValueError(architecture)
        self.architecture = architecture
        self.widths = ARCHITECTURES[architecture]
        self.seed = int(seed)
        self.action_scales = tf.constant(np.asarray(action_scales, np.float64), tf.float64)
        self.hidden = [
            tf.keras.layers.Dense(
                width,
                activation="tanh",
                dtype=tf.float64,
                kernel_initializer=tf.keras.initializers.GlorotUniform(seed + 101 + j),
                bias_initializer="zeros",
                name=f"hidden_{j}",
            )
            for j, width in enumerate(self.widths)
        ]
        self.output_layer = tf.keras.layers.Dense(
            1,
            dtype=tf.float64,
            kernel_initializer="zeros",
            bias_initializer="zeros",
            name="log_ratio",
        )

    def features(self, state: tf.Tensor) -> tf.Tensor:
        state = tf.cast(state, tf.float64)
        action = state[:, :3]
        h = action / (action + self.action_scales[None, :])
        h1, h2, h3 = tf.unstack(h, axis=-1)
        th1, th3 = state[:, 3], state[:, 4]
        return tf.stack(
            (
                h1, h2, h3,
                h1*h1, h2*h2, h3*h3,
                h1*h2, h2*h3, h1*h3,
                h1*tf.sin(th1), h1*tf.cos(th1),
                h1*h1*tf.sin(2.0*th1), h1*h1*tf.cos(2.0*th1),
                h3*tf.sin(th3), h3*tf.cos(th3),
                h3*h3*tf.sin(2.0*th3), h3*h3*tf.cos(2.0*th3),
                h1*h3*tf.sin(th1-th3), h1*h3*tf.cos(th1-th3),
                h1*h3*tf.sin(th1+th3), h1*h3*tf.cos(th1+th3),
            ),
            axis=-1,
        )

    def log_ratio(self, state: tf.Tensor) -> tf.Tensor:
        value = self.features(state)
        for layer in self.hidden:
            value = layer(value)
        return self.output_layer(value)[:, 0]

    def log_prob(self, state: tf.Tensor) -> tf.Tensor:
        return -energy_tf(state) / tf.constant(5.0, tf.float64) + self.log_ratio(state)

    def call(self, state: tf.Tensor) -> tf.Tensor:
        return self.log_ratio(state)

    def save_bundle(self, directory: Path, extra: dict) -> None:
        directory.mkdir(parents=True, exist_ok=True)
        _ = self(tf.constant([[1.0, 1.0, 1.0, 0.0, 0.0]], tf.float64))
        config = {
            "architecture": self.architecture,
            "action_scales": self.action_scales.numpy().tolist(),
            "seed": self.seed,
            **extra,
        }
        (directory / "model_config.json").write_text(
            json.dumps(config, indent=2, sort_keys=True) + "\n"
        )
        self.save_weights(directory / "model.weights.h5")


def load_bundle(directory: Path) -> RegularGibbsRatio:
    directory = Path(directory)
    config = json.loads((directory / "model_config.json").read_text())
    model = RegularGibbsRatio(
        config["architecture"], np.asarray(config["action_scales"], np.float64),
        int(config["seed"]),
    )
    _ = model(tf.constant([[1.0, 1.0, 1.0, 0.0, 0.0]], tf.float64))
    model.load_weights(directory / "model.weights.h5")
    return model
