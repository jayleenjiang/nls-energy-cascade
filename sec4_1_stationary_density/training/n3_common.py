#!/usr/bin/env python3
"""Three-mode reduced data, Hamiltonian, drift and FP residual (n = 3 port of n4_common.py).

State (long-chain convention, as in the n = 4 code):
    x = (I1, I2, I3, theta2, theta3),  theta_j = 2(phi_j - phi_{j-1}).
The paper's three-mode angle is theta1 = 2(phi1 - phi2) = -theta2.
"""
from __future__ import annotations

import struct
from pathlib import Path

import numpy as np
import tensorflow as tf

GAMMA = 0.1
D = 5


def read_nls5d(path) -> np.ndarray:
    with open(path, "rb") as f:
        if f.read(8) != b"NLS5D01\0":
            raise ValueError(f"bad header: {path}")
        n, d, b = struct.unpack("<QII", f.read(16))
        x = np.fromfile(f, "<f4")
    if d != D or b != 4 or x.size != n * D:
        raise ValueError(f"bad layout: {path}")
    x = x.reshape(n, D).astype(np.float64)
    if not np.isfinite(x).all() or np.any(x[:, :3] <= 0.0):
        raise ValueError(f"invalid state: {path}")
    return x


def write_nls5d(path, x) -> None:
    with open(path, "wb") as f:
        f.write(b"NLS5D01\0")
        f.write(struct.pack("<QII", len(x), D, 4))
        f.write(np.asarray(x, "<f4").tobytes())


def read_kde(path):
    with open(path, "rb") as f:
        if f.read(8) != b"NLS5KDE1":
            raise ValueError(f"bad header: {path}")
        q, n = struct.unpack("<QQ", f.read(16))
        d = np.fromfile(f, "<f8")
    if d.size != q:
        raise ValueError(f"truncated: {path}")
    return d, int(n)


def read_summary(path):
    vals = {}
    for line in Path(path).read_text().splitlines():
        toks = line.split()
        if toks and toks[0] in ("mean_I", "mean_current"):
            vals[toks[0]] = [float(v) for v in toks[1:]]
            continue
        for k, v in zip(toks[::2], toks[1::2]):
            vals[k] = v
    return vals


def energy_tf(state: tf.Tensor) -> tf.Tensor:
    """H/2, whose equilibrium density is exp[-energy/T]."""
    state = tf.cast(state, tf.float64)
    actions = state[:, :3]
    th2, th3 = tf.unstack(state[:, 3:], axis=-1)
    mass = tf.reduce_sum(actions, axis=-1)
    onsite = 0.5 * mass * mass - 0.25 * tf.reduce_sum(actions * actions, axis=-1)
    bonds = actions[:, 0] * actions[:, 1] * tf.cos(th2) + actions[:, 1] * actions[:, 2] * tf.cos(th3)
    return onsite + bonds


def energy_np(state: np.ndarray) -> np.ndarray:
    state = np.asarray(state, np.float64)
    a = state[:, :3]
    m = a.sum(axis=1)
    return (0.5 * m * m - 0.25 * np.sum(a * a, axis=1)
            + a[:, 0] * a[:, 1] * np.cos(state[:, 3]) + a[:, 1] * a[:, 2] * np.cos(state[:, 4]))


def drift_tf(state: tf.Tensor, t_left: float, t_right: float, gamma: float = GAMMA) -> tf.Tensor:
    state = tf.cast(state, tf.float64)
    i1, i2, i3, th2, th3 = tf.unstack(state, axis=-1)
    s2, s3, c2, c3 = tf.sin(th2), tf.sin(th3), tf.cos(th2), tf.cos(th3)
    mass = i1 + i2 + i3
    g = tf.cast(gamma, tf.float64)
    tl = tf.cast(t_left, tf.float64)
    tr = tf.cast(t_right, tf.float64)
    p1 = 2.0 * mass - i1 + 2.0 * i2 * c2
    p2 = 2.0 * mass - i2 + 2.0 * i1 * c2 + 2.0 * i3 * c3
    p3 = 2.0 * mass - i3 + 2.0 * i2 * c3
    b1 = -4.0 * i1 * i2 * s2 + 2.0 * g * (2.0 * tl - i1 * p1)
    b2 = 4.0 * i2 * (i1 * s2 - i3 * s3)
    b3 = 4.0 * i3 * i2 * s3 + 2.0 * g * (2.0 * tr - i3 * p3)
    p1 = p1 - 2.0 * g * i2 * s2
    p3 = p3 + 2.0 * g * i2 * s3
    return tf.stack((b1, b2, b3, 2.0 * (p2 - p1), 2.0 * (p3 - p2)), axis=-1)


def half_covariance(state: tf.Tensor, t_left: float, t_right: float, gamma: float = GAMMA) -> tf.Tensor:
    i1, _, i3, _, _ = tf.unstack(tf.cast(state, tf.float64), axis=-1)
    g = tf.cast(gamma, tf.float64)
    tl = tf.cast(t_left, tf.float64)
    tr = tf.cast(t_right, tf.float64)
    zero = tf.zeros_like(i1)
    return tf.stack((4.0 * g * tl * i1, zero, 4.0 * g * tr * i3, 4.0 * g * tl / i1, 4.0 * g * tr / i3), axis=-1)


def relative_adjoint_residual(model, state: tf.Tensor, t_left: float, t_right: float,
                              gamma: float = GAMMA) -> tf.Tensor:
    """Return (L* rho)/rho using automatic derivatives in all 5 variables (as n4_common)."""
    state = tf.cast(state, tf.float64)
    with tf.GradientTape(persistent=True) as second:
        second.watch(state)
        with tf.GradientTape(persistent=True) as first:
            first.watch(state)
            logp = model.log_prob(state)
            drift = drift_tf(state, t_left, t_right, gamma)
            logp_sum = tf.reduce_sum(logp)
            drift_columns = [tf.reduce_sum(drift[:, j]) for j in range(D)]
        grad = first.gradient(logp_sum, state)
        grad_columns = [tf.reduce_sum(grad[:, j]) for j in range(D)]
    divergence_terms = [first.gradient(drift_columns[j], state)[:, j] for j in range(D)]
    hdiag = [second.gradient(grad_columns[j], state)[:, j] for j in range(D)]
    del first, second
    divb = tf.add_n(divergence_terms)
    hdiag = tf.stack(hdiag, axis=-1)
    g = tf.cast(gamma, tf.float64)
    a = half_covariance(state, t_left, t_right, gamma)
    residual = -divb - tf.reduce_sum(drift * grad, axis=-1)
    residual += tf.reduce_sum(a * (hdiag + grad * grad), axis=-1)
    # d/dI of the action diffusion 4 gamma T I gives 2 * 4 gamma T * d log rho / dI
    residual += 8.0 * g * tf.cast(t_left, tf.float64) * grad[:, 0] + 8.0 * g * tf.cast(t_right, tf.float64) * grad[:, 2]
    return residual


class ExactGibbs(tf.Module):
    def __init__(self, temperature: float):
        super().__init__()
        self.temperature = tf.constant(temperature, tf.float64)

    def log_prob(self, state: tf.Tensor) -> tf.Tensor:
        return -energy_tf(state) / self.temperature
