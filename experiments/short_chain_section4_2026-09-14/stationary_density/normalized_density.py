#!/usr/bin/env python3
"""Normalized neural density for positive actions and two periodic angles."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import tensorflow as tf


LOG_TWO_PI = float(np.log(2.0*np.pi))


def inverse_softplus(value: np.ndarray) -> np.ndarray:
    value = np.asarray(value, dtype=np.float64)
    return np.log(np.expm1(value))


class NormalizedNLS3Density(tf.keras.Model):
    """GMM log-action density times a neural conditional von-Mises mixture."""

    def __init__(
        self,
        components: int,
        widths: tuple[int, ...],
        u_mean: np.ndarray,
        u_std: np.ndarray,
        gmm_weights: np.ndarray,
        gmm_means: np.ndarray,
        gmm_covariances: np.ndarray,
        seed: int,
    ) -> None:
        super().__init__(name="normalized_nls3_density")
        self.components = int(components)
        self.widths = tuple(int(item) for item in widths)
        self.seed = int(seed)
        self.u_mean = tf.constant(np.asarray(u_mean, dtype=np.float64), tf.float64)
        self.u_std = tf.constant(np.asarray(u_std, dtype=np.float64), tf.float64)

        chol = np.linalg.cholesky(np.asarray(gmm_covariances, dtype=np.float64))
        raw_chol = np.tril(chol)
        diag = np.diagonal(raw_chol, axis1=1, axis2=2)
        diag_raw = inverse_softplus(np.maximum(diag-1e-4, 1e-6))
        for k in range(self.components):
            np.fill_diagonal(raw_chol[k], diag_raw[k])

        self.gmm_logits = self.add_weight(
            name="gmm_logits", shape=(self.components,), dtype=tf.float64,
            initializer=tf.constant_initializer(np.log(np.maximum(gmm_weights, 1e-12))),
            trainable=True,
        )
        self.gmm_means = self.add_weight(
            name="gmm_means", shape=(self.components, 3), dtype=tf.float64,
            initializer=tf.constant_initializer(gmm_means), trainable=True,
        )
        self.gmm_raw_chol = self.add_weight(
            name="gmm_raw_chol", shape=(self.components, 3, 3), dtype=tf.float64,
            initializer=tf.constant_initializer(raw_chol), trainable=True,
        )

        self.hidden = []
        for index, width in enumerate(self.widths):
            self.hidden.append(tf.keras.layers.Dense(
                width, activation="tanh", dtype=tf.float64,
                kernel_initializer=tf.keras.initializers.GlorotUniform(seed+self.seed+index),
                name=f"conditional_hidden_{index}",
            ))
        kappa_bias = float(np.log((1.0/50.0)/(1.0-1.0/50.0)))
        output_bias = np.tile([0.0, 0.0, kappa_bias, 0.0, kappa_bias], self.components)
        self.conditional_output = tf.keras.layers.Dense(
            5*self.components, dtype=tf.float64,
            kernel_initializer=tf.keras.initializers.RandomNormal(
                mean=0.0, stddev=0.01, seed=seed+7919
            ),
            bias_initializer=tf.constant_initializer(output_bias),
            name="conditional_parameters",
        )

    def get_cholesky(self) -> tf.Tensor:
        lower = tf.linalg.band_part(self.gmm_raw_chol, -1, 0)
        diagonal = tf.nn.softplus(tf.linalg.diag_part(lower))+tf.constant(1e-4, tf.float64)
        return tf.linalg.set_diag(lower, diagonal)

    def conditional_parameters(self, z: tf.Tensor) -> tuple[tf.Tensor, ...]:
        h = tf.cast(z, tf.float64)
        for layer in self.hidden:
            h = layer(h)
        raw = tf.reshape(self.conditional_output(h), (-1, self.components, 5))
        logits = raw[..., 0]
        mu1 = np.pi*tf.tanh(raw[..., 1])
        kappa1 = 50.0*tf.sigmoid(raw[..., 2])
        mu3 = np.pi*tf.tanh(raw[..., 3])
        kappa3 = 50.0*tf.sigmoid(raw[..., 4])
        return logits, mu1, kappa1, mu3, kappa3

    @staticmethod
    def log_i0(kappa: tf.Tensor) -> tf.Tensor:
        return tf.math.log(tf.math.bessel_i0e(kappa))+tf.math.abs(kappa)

    def log_prob_z(self, z: tf.Tensor) -> tf.Tensor:
        z = tf.cast(z, tf.float64)
        difference = z[:, None, :]-self.gmm_means[None, :, :]
        chol = self.get_cholesky()
        solved = tf.linalg.triangular_solve(
            chol[None, :, :, :], difference[:, :, :, None], lower=True
        )[..., 0]
        mahalanobis = tf.reduce_sum(solved*solved, axis=-1)
        logdet = tf.reduce_sum(tf.math.log(tf.linalg.diag_part(chol)), axis=-1)
        component = -0.5*(3.0*LOG_TWO_PI+mahalanobis)-logdet[None, :]
        return tf.reduce_logsumexp(
            tf.nn.log_softmax(self.gmm_logits)[None, :]+component, axis=-1
        )

    def log_prob_angles(self, z: tf.Tensor, theta: tf.Tensor) -> tf.Tensor:
        logits, mu1, kappa1, mu3, kappa3 = self.conditional_parameters(z)
        th1 = tf.cast(theta[:, 0:1], tf.float64)
        th3 = tf.cast(theta[:, 1:2], tf.float64)
        log_vm1 = kappa1*tf.cos(th1-mu1)-LOG_TWO_PI-self.log_i0(kappa1)
        log_vm3 = kappa3*tf.cos(th3-mu3)-LOG_TWO_PI-self.log_i0(kappa3)
        return tf.reduce_logsumexp(
            tf.nn.log_softmax(logits, axis=-1)+log_vm1+log_vm3, axis=-1
        )

    def log_prob(self, state: tf.Tensor) -> tf.Tensor:
        state = tf.cast(state, tf.float64)
        actions = state[:, :3]
        # Positive support is enforced by the model.  Invalid evaluation points
        # receive -infinity rather than being projected into the support.
        safe_actions = tf.where(actions > 0, actions, tf.ones_like(actions))
        u = tf.math.log(safe_actions)
        z = (u-self.u_mean)/self.u_std
        log_prob = self.log_prob_z(z)-tf.reduce_sum(tf.math.log(self.u_std), axis=-1)
        log_prob += self.log_prob_angles(z, state[:, 3:])
        log_prob -= tf.reduce_sum(tf.math.log(safe_actions), axis=-1)
        valid = tf.reduce_all(actions > 0, axis=-1)
        return tf.where(valid, log_prob,
                        tf.fill(tf.shape(log_prob), tf.constant(-np.inf, tf.float64)))

    def call(self, state: tf.Tensor) -> tf.Tensor:
        return self.log_prob(state)

    def sample_numpy(self, count: int, seed: int, batch: int = 65_536) -> np.ndarray:
        """Draw independent normalized model samples in reduced coordinates."""
        rng = np.random.default_rng(seed)
        result = np.empty((int(count), 5), dtype=np.float64)
        weights = tf.nn.softmax(self.gmm_logits).numpy()
        means = self.gmm_means.numpy()
        chol = self.get_cholesky().numpy()
        u_mean = self.u_mean.numpy()
        u_std = self.u_std.numpy()
        for start in range(0, int(count), int(batch)):
            stop = min(int(count), start+int(batch))
            n = stop-start
            gmm_component = rng.choice(self.components, size=n, p=weights)
            normal = rng.normal(size=(n, 3))
            z = means[gmm_component] + np.einsum(
                "nij,nj->ni", chol[gmm_component], normal
            )
            u = u_mean+u_std*z
            result[start:stop, :3] = np.exp(u)
            logits, mu1, kappa1, mu3, kappa3 = self.conditional_parameters(
                tf.constant(z, tf.float64)
            )
            conditional_weights = tf.nn.softmax(logits, axis=-1).numpy()
            uniforms = rng.random(n)
            conditional_component = np.sum(
                uniforms[:, None] > np.cumsum(conditional_weights, axis=1), axis=1
            )
            rows = np.arange(n)
            result[start:stop, 3] = rng.vonmises(
                mu1.numpy()[rows, conditional_component],
                kappa1.numpy()[rows, conditional_component],
            )
            result[start:stop, 4] = rng.vonmises(
                mu3.numpy()[rows, conditional_component],
                kappa3.numpy()[rows, conditional_component],
            )
        return result

    def config_dict(self) -> dict:
        return {
            "components": self.components,
            "widths": list(self.widths),
            "seed": self.seed,
            "u_mean": self.u_mean.numpy().tolist(),
            "u_std": self.u_std.numpy().tolist(),
        }

    def save_bundle(self, directory: Path, extra: dict | None = None) -> None:
        directory.mkdir(parents=True, exist_ok=True)
        config = self.config_dict()
        if extra:
            config.update(extra)
        (directory/"model_config.json").write_text(json.dumps(config, indent=2, sort_keys=True)+"\n")
        # Keras subclass bookkeeping is not set when training invokes only the
        # component log-probability methods.  Build the unchanged call graph
        # once before serialization; this does not modify trainable values.
        dummy = tf.constant([[1.0, 1.0, 1.0, 0.0, 0.0]], tf.float64)
        _ = self(dummy)
        self.save_weights(directory/"model.weights.h5")


def load_bundle(directory: Path) -> NormalizedNLS3Density:
    """Reconstruct a saved model without refitting any preprocessing."""
    directory = Path(directory)
    config = json.loads((directory/"model_config.json").read_text())
    initial = np.load(directory/"initial_gmm_and_scaler.npz")
    model = NormalizedNLS3Density(
        components=int(config["components"]),
        widths=tuple(int(x) for x in config["widths"]),
        u_mean=np.asarray(initial["u_mean"], dtype=np.float64),
        u_std=np.asarray(initial["u_std"], dtype=np.float64),
        gmm_weights=np.asarray(initial["weights"], dtype=np.float64),
        gmm_means=np.asarray(initial["means"], dtype=np.float64),
        gmm_covariances=np.asarray(initial["covariances"], dtype=np.float64),
        seed=int(config["seed"]),
    )
    dummy = tf.constant([[1.0, 1.0, 1.0, 0.0, 0.0]], tf.float64)
    _ = model.log_prob(dummy)
    model.load_weights(directory/"model.weights.h5")
    return model
