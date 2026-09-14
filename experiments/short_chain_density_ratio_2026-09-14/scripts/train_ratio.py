#!/usr/bin/env python3
"""Train one preregistered Gibbs-anchored density-ratio fit."""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
import tensorflow as tf
from scipy.special import logsumexp
from sklearn.metrics import roc_auc_score

from density_ratio_common import fixed_indices, load_rows, target_reference_rows
from gibbs_ratio_model import ARCHITECTURES, GibbsRatioDensity

BASE = Path(__file__).resolve().parent.parent
OLD = BASE.parent/"short_chain_section4_2026-09-14"
import sys
sys.path.insert(0, str(OLD/"stationary_density"))
from fp_tensorflow import relative_adjoint_residual  # noqa: E402


SEEDS = (5201, 5202, 5203)


def args_parser():
    p = argparse.ArgumentParser()
    p.add_argument("--phase", choices=("equilibrium", "driven"), required=True)
    p.add_argument("--timestep", choices=("fine", "coarse"), default="fine")
    p.add_argument("--architecture", choices=tuple(ARCHITECTURES), required=True)
    p.add_argument("--seed", type=int, choices=SEEDS, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--epochs", type=int, default=100)
    p.add_argument("--steps-per-epoch", type=int, default=32)
    p.add_argument("--batch", type=int, default=1024)
    p.add_argument("--fp-batch", type=int, default=256)
    p.add_argument("--warmup-epochs", type=int, default=20)
    p.add_argument("--patience", type=int, default=15)
    p.add_argument("--validation-points", type=int, default=100_000)
    p.add_argument("--fp-points", type=int, default=10_000)
    p.add_argument("--smoke", action="store_true")
    return p.parse_args()


def bce(model, target, reference, batch=8192):
    total = 0.0
    count = 0
    for start in range(0, len(target), batch):
        value = model.log_ratio(tf.constant(target[start:start+batch], tf.float64))
        total += float(tf.reduce_sum(tf.nn.softplus(-value)).numpy())
        count += len(value)
    for start in range(0, len(reference), batch):
        value = model.log_ratio(tf.constant(reference[start:start+batch], tf.float64))
        total += float(tf.reduce_sum(tf.nn.softplus(value)).numpy())
        count += len(value)
    return total/count*2.0


def log_ratios(model, state, batch=8192):
    out = []
    for start in range(0, len(state), batch):
        out.append(model.log_ratio(tf.constant(state[start:start+batch], tf.float64)).numpy())
    return np.concatenate(out)


def residuals(model, state, t1, t3, batch=256):
    out = []
    for start in range(0, len(state), batch):
        out.append(relative_adjoint_residual(
            model, tf.constant(state[start:start+batch], tf.float64), t1, t3, 0.1
        ).numpy())
    return np.concatenate(out)


def main():
    a = args_parser()
    started = time.time()
    tf.keras.backend.set_floatx("float64")
    tf.random.set_seed(a.seed)
    rng = np.random.default_rng(a.seed)
    target_rows, reference_rows = target_reference_rows(a.phase, a.timestep, "train")
    val_target_rows, val_reference_rows = target_reference_rows(a.phase, a.timestep, "validation")
    # A smoke test checks one file per role and is explicitly excluded from the
    # scientific pipeline.  Truncate the file lists before I/O so the test does
    # not hydrate or scan the full frozen dataset.
    if a.smoke:
        target_rows = target_rows[:1]; reference_rows = reference_rows[:1]
        val_target_rows = val_target_rows[:1]; val_reference_rows = val_reference_rows[:1]
    target, _ = load_rows(target_rows)
    reference, _ = load_rows(reference_rows)
    val_target, _ = load_rows(val_target_rows)
    val_reference, _ = load_rows(val_reference_rows)
    if a.smoke:
        target = target[:4096]; reference = reference[:4096]
        val_target = val_target[:2048]; val_reference = val_reference[:2048]
        a.epochs = 1; a.steps_per_epoch = 1; a.warmup_epochs = 0
        a.validation_points = 1024; a.fp_points = 256

    scale_index = fixed_indices(len(reference), min(200_000, len(reference)),
                                f"scale|{a.phase}|{a.timestep}")
    action_scales = np.median(reference[scale_index, :3], axis=0)
    model = GibbsRatioDensity(a.architecture, action_scales, a.seed)
    _ = model(tf.constant(target[:2], tf.float64))
    optimizer = tf.keras.optimizers.Adam(1e-3, global_clipnorm=10.0)
    t1, t3 = ((5.0, 5.0) if a.phase == "equilibrium" else (2.0, 8.0))

    val_target = val_target[fixed_indices(
        len(val_target), a.validation_points,
        f"validation-target|{a.phase}|{a.timestep}")]
    val_reference = val_reference[fixed_indices(
        len(val_reference), a.validation_points,
        f"validation-reference|{a.phase}|{a.timestep}")]

    @tf.function(reduce_retracing=True)
    def step(p, q, xfp, physics_weight):
        with tf.GradientTape() as tape:
            fp = model.log_ratio(p); fq = model.log_ratio(q)
            classification = tf.reduce_mean(tf.nn.softplus(-fp))+tf.reduce_mean(tf.nn.softplus(fq))
            rr = relative_adjoint_residual(model, xfp, t1, t3, 0.1)
            physics = tf.reduce_mean(rr*rr)
            regularization = tf.add_n([tf.reduce_sum(v*v) for v in model.trainable_variables])
            loss = classification+physics_weight*physics+tf.constant(1e-6, tf.float64)*regularization
        gradients = tape.gradient(loss, model.trainable_variables)
        optimizer.apply_gradients(zip(gradients, model.trainable_variables))
        return classification, physics, loss

    best_validation = bce(model, val_target, val_reference)
    best_weights = model.get_weights()
    best_epoch = 0
    stale = 0
    history = [{"epoch": 0, "validation_bce": best_validation}]
    for epoch in range(1, a.epochs+1):
        sums = np.zeros(3)
        weight = 0.0 if epoch <= a.warmup_epochs else 0.01
        for _ in range(a.steps_per_epoch):
            p = target[rng.integers(0, len(target), a.batch)]
            q = reference[rng.integers(0, len(reference), a.batch)]
            xfp = target[rng.integers(0, len(target), a.fp_batch)]
            values = step(tf.constant(p, tf.float64), tf.constant(q, tf.float64),
                          tf.constant(xfp, tf.float64), tf.constant(weight, tf.float64))
            sums += [float(v.numpy()) for v in values]
        validation = bce(model, val_target, val_reference)
        history.append({
            "epoch": epoch, "classification": sums[0]/a.steps_per_epoch,
            "physics": sums[1]/a.steps_per_epoch, "loss": sums[2]/a.steps_per_epoch,
            "validation_bce": validation,
        })
        if validation < best_validation-1e-7:
            best_validation = validation; best_epoch = epoch
            best_weights = model.get_weights(); stale = 0
        else:
            stale += 1
        if epoch >= a.warmup_epochs and stale >= a.patience:
            break
    model.set_weights(best_weights)

    ft = log_ratios(model, val_target)
    fq = log_ratios(model, val_reference)
    logz = float(logsumexp(fq)-math.log(len(fq)))
    normalized_reference = fq-logz
    weights = np.exp(normalized_reference)
    ess = float(weights.sum()**2/np.sum(weights*weights))
    auc = float(roc_auc_score(
        np.r_[np.ones(len(ft)), np.zeros(len(fq))], np.r_[ft, fq]
    ))
    fp_state = val_target[fixed_indices(
        len(val_target), a.fp_points, f"fp|{a.phase}|{a.timestep}")]
    rr = residuals(model, fp_state, t1, t3)
    normalized_target = ft-logz
    metrics = {
        "protocol_version": "section4-density-ratio-v1",
        "phase": a.phase, "timestep": a.timestep,
        "architecture": a.architecture, "seed": a.seed,
        "action_scales": action_scales.tolist(),
        "best_epoch": best_epoch, "epochs_completed": len(history)-1,
        "validation_bce": best_validation,
        "validation_auc": auc,
        "validation_logz": logz,
        "centered_log_ratio_rms": float(np.sqrt(np.mean(normalized_target**2))),
        "reference_weight_ess": ess,
        "reference_weight_ess_fraction": ess/len(weights),
        "max_weight_share": float(weights.max()/weights.sum()),
        "fp_abs_median": float(np.median(np.abs(rr))),
        "fp_abs_p90": float(np.quantile(np.abs(rr), 0.90)),
        "finite": bool(np.isfinite(np.r_[ft, fq, rr, weights]).all()),
        "target_train_streams": [int(x["stream_id"]) for x in target_rows],
        "reference_train_streams": [int(x["stream_id"]) for x in reference_rows],
        "target_validation_streams": [int(x["stream_id"]) for x in val_target_rows],
        "reference_validation_streams": [int(x["stream_id"]) for x in val_reference_rows],
        "runtime_seconds": time.time()-started,
        "smoke": bool(a.smoke),
    }
    if a.phase == "equilibrium":
        metrics["validation_admissible"] = bool(
            metrics["centered_log_ratio_rms"] <= 0.05
            and 0.48 <= auc <= 0.52
            and metrics["reference_weight_ess_fraction"] >= 0.50
            and metrics["fp_abs_median"] <= 0.10
            and metrics["fp_abs_p90"] <= 0.50
            and metrics["finite"]
        )
    a.output.mkdir(parents=True, exist_ok=True)
    model.save_bundle(a.output, {"phase": a.phase, "timestep": a.timestep,
                                 "validation_logz": logz})
    (a.output/"history.json").write_text(json.dumps(history, indent=2)+"\n")
    (a.output/"validation_metrics.json").write_text(json.dumps(metrics, indent=2, sort_keys=True)+"\n")
    print(json.dumps(metrics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
