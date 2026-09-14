#!/usr/bin/env python3
"""Small differentiability smoke test for the normalized density model."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import tensorflow as tf


HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE/"stationary_density"))
from fp_tensorflow import relative_adjoint_residual  # noqa: E402
from normalized_density import NormalizedNLS3Density  # noqa: E402


def main() -> None:
    tf.keras.backend.set_floatx("float64")
    rng = np.random.default_rng(1234)
    components = 2
    model = NormalizedNLS3Density(
        components=components,
        widths=(8,),
        u_mean=np.zeros(3),
        u_std=np.ones(3),
        gmm_weights=np.full(components, 1/components),
        gmm_means=np.array([[-0.5, 0.0, 0.5], [0.5, 0.0, -0.5]]),
        gmm_covariances=np.repeat(np.eye(3)[None, :, :], components, axis=0),
        seed=7,
    )
    state = np.column_stack((
        np.exp(rng.normal(size=(32, 3))),
        rng.uniform(-np.pi, np.pi, size=(32, 2)),
    ))
    tensor = tf.constant(state, tf.float64)
    # Build all Keras variables before differentiating with respect to state.
    _ = model.log_prob(tensor[:1])
    optimizer = tf.keras.optimizers.Adam(1e-4)
    with tf.GradientTape() as tape:
        logp = model.log_prob(tensor)
        residual = relative_adjoint_residual(model, tensor, 5.0, 5.0, 0.1)
        loss = -tf.reduce_mean(logp)+1e-6*tf.reduce_mean(residual*residual)
    gradients = tape.gradient(loss, model.trainable_variables)
    finite_gradients = all(g is not None and bool(tf.reduce_all(tf.math.is_finite(g))) for g in gradients)
    optimizer.apply_gradients(zip(gradients, model.trainable_variables))
    sample = model.sample_numpy(10_000, seed=991)
    result = {
        "finite_log_prob": bool(tf.reduce_all(tf.math.is_finite(logp))),
        "finite_residual": bool(tf.reduce_all(tf.math.is_finite(residual))),
        "finite_complete_gradients": finite_gradients,
        "trainable_variables": len(model.trainable_variables),
        "loss": float(loss),
        "residual_rms": float(tf.sqrt(tf.reduce_mean(residual*residual))),
        "finite_samples": bool(np.isfinite(sample).all()),
        "positive_sampled_actions": bool(np.all(sample[:, :3] > 0)),
        "sampled_angles_in_principal_interval": bool(
            np.all(sample[:, 3:] >= -np.pi) and np.all(sample[:, 3:] <= np.pi)
        ),
    }
    result["pass"] = all(result[key] for key in
                         ("finite_log_prob", "finite_residual", "finite_complete_gradients",
                          "finite_samples", "positive_sampled_actions",
                          "sampled_angles_in_principal_interval"))
    out = HERE/"stationary_density"/"MODEL_SMOKE.json"
    out.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    if not result["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
