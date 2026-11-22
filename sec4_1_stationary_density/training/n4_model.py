#!/usr/bin/env python3
"""Boundary-regular periodic Gibbs density ratio for the four-mode chain."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import tensorflow as tf

from n4_common import energy_tf


class N4GibbsRatio(tf.keras.Model):
    def __init__(self, action_scales: np.ndarray, seed: int, width: int = 128,
                 feature_version: str = "transport_v2", depth: int = 2,
                 linear_skip: bool = False):
        super().__init__(name="n4_gibbs_ratio")
        self.seed = int(seed)
        self.width = int(width)
        self.depth = int(depth)
        if feature_version not in ("regular_v1", "transport_v2", "transport_v3",
                                   "phase_regular_v4"):
            raise ValueError(feature_version)
        self.feature_version = feature_version
        self.linear_skip = bool(linear_skip)
        self.action_scales = tf.constant(np.asarray(action_scales, np.float64), tf.float64)
        self.hidden = [
            tf.keras.layers.Dense(
                width, activation="tanh", dtype=tf.float64,
                kernel_initializer=tf.keras.initializers.GlorotUniform(seed + 101 + j),
                bias_initializer="zeros", name=f"hidden_{j}",
            ) for j in range(self.depth)
        ]
        self.output_layer = tf.keras.layers.Dense(
            1, dtype=tf.float64, kernel_initializer="zeros",
            bias_initializer="zeros", name="log_ratio",
        )
        if self.linear_skip:
            # A linear action tilt remains integrable against the quadratic
            # Gibbs tail. Bounded phase terms retain the required endpoint
            # vanishing orders. This bypasses tanh saturation at large action.
            self.linear_skip_weights = self.add_weight(
                name="linear_skip_weights", shape=(10,), initializer="zeros",
                dtype=tf.float64, trainable=True)

    def features(self, state: tf.Tensor) -> tf.Tensor:
        state = tf.cast(state, tf.float64)
        action = state[:, :4]
        h = action / (action + self.action_scales[None, :])
        h1, h2, h3, h4 = tf.unstack(h, axis=-1)
        t2, t3, t4 = tf.unstack(state[:, 4:], axis=-1)
        # Endpoint harmonics vanish where their stochastic polar phase is
        # singular.  The middle phase is multiplied by the two internal gates.
        if self.feature_version == "phase_regular_v4":
            # theta_j = 2(phi_j-phi_{j-1}).  A first harmonic is smooth in
            # the Cartesian amplitudes only when it vanishes linearly at
            # both participating actions; a second harmonic needs order two.
            g2, g3, g4 = h1*h2, h2*h3, h3*h4
        else:
            g2 = h1
            g3 = tf.ones_like(h1) if self.feature_version == "transport_v3" else h2*h3
            g4 = h4
        values = [h1, h2, h3, h4]
        values += [h1*h1, h2*h2, h3*h3, h4*h4]
        values += [h1*h2, h2*h3, h3*h4, h1*h3, h2*h4, h1*h4]
        for gate, angle in ((g2, t2), (g3, t3), (g4, t4)):
            values += [gate*tf.sin(angle), gate*tf.cos(angle),
                       gate*gate*tf.sin(2.0*angle), gate*gate*tf.cos(2.0*angle)]
        if self.feature_version == "phase_regular_v4":
            # For combinations, use the charge of every original phase.
            # theta2+theta3, for example, is 2(phi3-phi1), so I2 must not
            # enter its vanishing order; theta2-theta3 needs I2 squared.
            pair_terms = (
                (h1*h2*h2*h3, h1*h3, t2, t3),
                (h2*h3*h3*h4, h2*h4, t3, t4),
                (h1*h2*h3*h4, h1*h2*h3*h4, t2, t4),
            )
            for minus_gate, plus_gate, a, b in pair_terms:
                values += [minus_gate*tf.sin(a-b), minus_gate*tf.cos(a-b),
                           plus_gate*tf.sin(a+b), plus_gate*tf.cos(a+b)]
        else:
            for ga, gb, a, b in ((g2, g3, t2, t3), (g3, g4, t3, t4),
                                  (g2, g4, t2, t4)):
                values += [ga*gb*tf.sin(a-b), ga*gb*tf.cos(a-b),
                           ga*gb*tf.sin(a+b), ga*gb*tf.cos(a+b)]
        regular = tf.stack(values, axis=-1)
        if self.feature_version == "regular_v1":
            return regular
        u = tf.math.log1p(action / self.action_scales[None, :])
        u1, u2, u3, u4 = tf.unstack(u, axis=-1)
        action_features = tf.stack((
            u1, u2, u3, u4, u1*u1, u2*u2, u3*u3, u4*u4,
            u1*u2, u2*u3, u3*u4, u1*u3, u2*u4, u1*u4,
        ), axis=-1)
        s1, s2, s3, s4 = tf.unstack(self.action_scales)
        transport = tf.stack((
            action[:, 0]*action[:, 1]*tf.sin(t2)/(s1*s2),
            action[:, 0]*action[:, 1]*tf.cos(t2)/(s1*s2),
            action[:, 1]*action[:, 2]*tf.sin(t3)/(s2*s3),
            action[:, 1]*action[:, 2]*tf.cos(t3)/(s2*s3),
            action[:, 2]*action[:, 3]*tf.sin(t4)/(s3*s4),
            action[:, 2]*action[:, 3]*tf.cos(t4)/(s3*s4),
        ), axis=-1)
        return tf.concat((regular, action_features, transport), axis=-1)

    def log_ratio(self, state: tf.Tensor) -> tf.Tensor:
        value = self.features(state)
        for layer in self.hidden:
            value = layer(value)
        output = self.output_layer(value)[:, 0]
        if self.linear_skip:
            state = tf.cast(state, tf.float64)
            action = state[:, :4]
            h = action / (action + self.action_scales[None, :])
            t2, t3, t4 = tf.unstack(state[:, 4:], axis=-1)
            g2, g3, g4 = h[:, 0]*h[:, 1], h[:, 1]*h[:, 2], h[:, 2]*h[:, 3]
            skip = tf.stack((action[:, 0]/self.action_scales[0],
                             action[:, 1]/self.action_scales[1],
                             action[:, 2]/self.action_scales[2],
                             action[:, 3]/self.action_scales[3],
                             g2*tf.sin(t2), g2*tf.cos(t2),
                             g3*tf.sin(t3), g3*tf.cos(t3),
                             g4*tf.sin(t4), g4*tf.cos(t4)), axis=-1)
            output += tf.reduce_sum(skip*self.linear_skip_weights[None, :], axis=-1)
        return output

    def log_prob(self, state: tf.Tensor) -> tf.Tensor:
        return -energy_tf(state) / tf.constant(5.0, tf.float64) + self.log_ratio(state)

    def call(self, state: tf.Tensor) -> tf.Tensor:
        return self.log_ratio(state)

    def save_bundle(self, directory: Path, extra: dict) -> None:
        directory.mkdir(parents=True, exist_ok=True)
        _ = self(tf.constant([[1., 1., 1., 1., 0., 0., 0.]], tf.float64))
        config = {
            "action_scales": self.action_scales.numpy().tolist(),
            "seed": self.seed, "width": self.width,
            "depth": self.depth, "feature_version": self.feature_version,
            "linear_skip": self.linear_skip, **extra,
        }
        (directory / "model_config.json").write_text(
            json.dumps(config, indent=2, sort_keys=True) + "\n")
        self.save_weights(directory / "model.weights.h5")


def load_bundle(directory: Path) -> N4GibbsRatio:
    directory = Path(directory)
    config = json.loads((directory / "model_config.json").read_text())
    model = N4GibbsRatio(np.asarray(config["action_scales"]),
                         int(config["seed"]), int(config["width"]),
                         config.get("feature_version", "regular_v1"),
                         int(config.get("depth", 2)),
                         bool(config.get("linear_skip", False)))
    _ = model(tf.constant([[1., 1., 1., 1., 0., 0., 0.]], tf.float64))
    model.load_weights(directory / "model.weights.h5")
    return model
