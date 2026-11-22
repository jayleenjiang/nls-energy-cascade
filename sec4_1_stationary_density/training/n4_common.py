#!/usr/bin/env python3
"""Four-mode reduced data, Hamiltonian, drift and FP residual."""

from __future__ import annotations

import os
import struct
from pathlib import Path

import numpy as np
import tensorflow as tf


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REPO = ROOT.parents[1]
SOURCE = Path(os.environ.get(
    "NLS_N4_RAW_SOURCE",
    REPO / "experiments/n4_7d_section4_1_2026-09-25/raw",
))
GAMMA = 0.1


def read_binary(path: Path) -> np.ndarray:
    with path.open("rb") as handle:
        if handle.read(8) != b"NLS7D01\0":
            raise ValueError(f"bad NLS7D header: {path}")
        count, dimension, float_bytes = struct.unpack("<QII", handle.read(16))
        if dimension != 7 or float_bytes != 4:
            raise ValueError(f"unsupported NLS7D layout: {path}")
        values = np.fromfile(handle, dtype="<f4")
    if values.size != count * dimension:
        raise ValueError(f"truncated NLS7D file: {path}")
    state = values.reshape(count, dimension).astype(np.float64)
    if not np.isfinite(state).all() or np.any(state[:, :4] <= 0.0):
        raise ValueError(f"invalid state: {path}")
    return state


def data_path(condition: str, timestep: str, seed: int) -> Path:
    if condition not in ("eq", "neq"):
        raise ValueError(condition)
    if timestep == "fine":
        tag = "fine_burn20"
    elif timestep == "coarse":
        tag = "coarse"
    else:
        raise ValueError(timestep)
    return SOURCE / f"{condition}_{tag}_seed{seed}.bin"


def load(condition: str, timestep: str, seed: int) -> np.ndarray:
    return read_binary(data_path(condition, timestep, seed))


def trajectories(state: np.ndarray) -> list[np.ndarray]:
    """Recover the 2 OpenMP-thread x 16 SIMD-lane trajectories."""
    if state.shape != (1_000_000, 7):
        raise ValueError(f"unexpected state shape {state.shape}")
    streams = []
    for thread_block in np.split(state, 2):
        shaped = thread_block.reshape(-1, 16, 7)
        streams.extend([shaped[:, lane, :] for lane in range(16)])
    return streams


def deterministic_indices(length: int, count: int, seed: int) -> np.ndarray:
    return np.random.default_rng(seed).choice(length, min(length, count), replace=False)


def energy_tf(state: tf.Tensor) -> tf.Tensor:
    """H/2, whose equilibrium density is exp[-energy/T]."""
    state = tf.cast(state, tf.float64)
    actions = state[:, :4]
    th2, th3, th4 = tf.unstack(state[:, 4:], axis=-1)
    mass = tf.reduce_sum(actions, axis=-1)
    onsite = 0.5 * mass * mass - 0.25 * tf.reduce_sum(actions * actions, axis=-1)
    bonds = (actions[:, 0] * actions[:, 1] * tf.cos(th2)
             + actions[:, 1] * actions[:, 2] * tf.cos(th3)
             + actions[:, 2] * actions[:, 3] * tf.cos(th4))
    return onsite + bonds


def energy_np(state: np.ndarray) -> np.ndarray:
    state = np.asarray(state, np.float64)
    actions = state[:, :4]
    mass = actions.sum(axis=1)
    return (0.5 * mass * mass - 0.25 * np.sum(actions * actions, axis=1)
            + actions[:, 0] * actions[:, 1] * np.cos(state[:, 4])
            + actions[:, 1] * actions[:, 2] * np.cos(state[:, 5])
            + actions[:, 2] * actions[:, 3] * np.cos(state[:, 6]))


def drift_tf(state: tf.Tensor, t_left: float, t_right: float,
             gamma: float = GAMMA) -> tf.Tensor:
    state = tf.cast(state, tf.float64)
    i1, i2, i3, i4, th2, th3, th4 = tf.unstack(state, axis=-1)
    s2, s3, s4 = tf.sin(th2), tf.sin(th3), tf.sin(th4)
    c2, c3, c4 = tf.cos(th2), tf.cos(th3), tf.cos(th4)
    mass = i1 + i2 + i3 + i4
    g = tf.cast(gamma, tf.float64)
    tl = tf.cast(t_left, tf.float64)
    tr = tf.cast(t_right, tf.float64)

    p1 = 2.0 * mass - i1 + 2.0 * i2 * c2
    p2 = 2.0 * mass - i2 + 2.0 * i1 * c2 + 2.0 * i3 * c3
    p3 = 2.0 * mass - i3 + 2.0 * i2 * c3 + 2.0 * i4 * c4
    p4 = 2.0 * mass - i4 + 2.0 * i3 * c4

    b1 = -4.0 * i1 * i2 * s2 + 2.0 * g * (2.0 * tl - i1 * p1)
    b2 = 4.0 * i2 * (i1 * s2 - i3 * s3)
    b3 = 4.0 * i3 * (i2 * s3 - i4 * s4)
    b4 = 4.0 * i4 * i3 * s4 + 2.0 * g * (2.0 * tr - i4 * p4)
    p1 = p1 - 2.0 * g * i2 * s2
    p4 = p4 + 2.0 * g * i3 * s4
    return tf.stack((b1, b2, b3, b4,
                     2.0 * (p2 - p1), 2.0 * (p3 - p2), 2.0 * (p4 - p3)),
                    axis=-1)


def relative_adjoint_residual(model, state: tf.Tensor, t_left: float,
                              t_right: float, gamma: float = GAMMA) -> tf.Tensor:
    """Return (L* rho)/rho using automatic derivatives in all 7 variables."""
    state = tf.cast(state, tf.float64)
    with tf.GradientTape(persistent=True) as second:
        second.watch(state)
        with tf.GradientTape(persistent=True) as first:
            first.watch(state)
            logp = model.log_prob(state)
            drift = drift_tf(state, t_left, t_right, gamma)
            logp_sum = tf.reduce_sum(logp)
            drift_columns = [tf.reduce_sum(drift[:, j]) for j in range(7)]
        grad = first.gradient(logp_sum, state)
        grad_columns = [tf.reduce_sum(grad[:, j]) for j in range(7)]
    divergence_terms = [first.gradient(drift_columns[j], state)[:, j]
                        for j in range(7)]
    hdiag = [second.gradient(grad_columns[j], state)[:, j] for j in range(7)]
    del first, second
    divb = tf.add_n(divergence_terms)
    hdiag = tf.stack(hdiag, axis=-1)

    i1, _, _, i4, _, _, _ = tf.unstack(state, axis=-1)
    g = tf.cast(gamma, tf.float64)
    tl = tf.cast(t_left, tf.float64)
    tr = tf.cast(t_right, tf.float64)
    zero = tf.zeros_like(i1)
    half_covariance = tf.stack((
        4.0 * g * tl * i1, zero, zero, 4.0 * g * tr * i4,
        4.0 * g * tl / i1, zero, 4.0 * g * tr / i4,
    ), axis=-1)
    residual = -divb - tf.reduce_sum(drift * grad, axis=-1)
    residual += tf.reduce_sum(half_covariance * (hdiag + grad * grad), axis=-1)
    residual += 8.0 * g * tl * grad[:, 0] + 8.0 * g * tr * grad[:, 3]
    return residual


class ExactGibbs(tf.Module):
    def __init__(self, temperature: float):
        super().__init__()
        self.temperature = tf.constant(temperature, tf.float64)

    def log_prob(self, state: tf.Tensor) -> tf.Tensor:
        return -energy_tf(state) / self.temperature
