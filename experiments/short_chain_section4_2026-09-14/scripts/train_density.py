#!/usr/bin/env python3
"""Train one frozen normalized-density candidate on trajectory-level splits."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import tensorflow as tf
from sklearn.mixture import GaussianMixture


HERE = Path(__file__).resolve().parent.parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE/"stationary_density"))

from fp_tensorflow import relative_adjoint_residual  # noqa: E402
from normalized_density import NormalizedNLS3Density  # noqa: E402
from reduced_operator import energy  # noqa: E402

N_SNAP = 100_001
N_COL = 5
SPLITS = HERE/"DIRECT_STATE_SPLITS.csv"
GLOBAL_SELECTION_SEED = 2026091403


def case_parameters(case: str) -> tuple[float, float]:
    if case.startswith("equilibrium"):
        return 5.0, 5.0
    if case.startswith("driven"):
        return 2.0, 8.0
    raise ValueError(case)


def stable_case_seed(case: str) -> int:
    return int(hashlib.sha256(case.encode()).hexdigest()[:8], 16)^GLOBAL_SELECTION_SEED


def load_stream(path: Path) -> np.ndarray:
    return np.fromfile(path, dtype="<f8").reshape(N_SNAP, N_COL)


def load_split(case: str, split: str) -> np.ndarray:
    with SPLITS.open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle)
                if row["branch"] == "density" and row["case"] == case and row["split"] == split]
    if not rows:
        raise RuntimeError(f"empty split {case}/{split}")
    arrays = [load_stream(REPO/row["relative_path"])
              for row in sorted(rows, key=lambda item: int(item["stream_id"]))]
    return np.concatenate(arrays, axis=0)


def iter_log_prob(model, state: np.ndarray, batch: int = 8192) -> np.ndarray:
    values = []
    for start in range(0, len(state), batch):
        values.append(model.log_prob(tf.constant(state[start:start+batch], tf.float64)).numpy())
    return np.concatenate(values)


def iter_residual(model, state: np.ndarray, t1: float, t3: float,
                  batch: int = 512) -> np.ndarray:
    values = []
    for start in range(0, len(state), batch):
        x = tf.constant(state[start:start+batch], tf.float64)
        values.append(relative_adjoint_residual(model, x, t1, t3, 0.1).numpy())
    return np.concatenate(values)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("--case", required=True)
    result.add_argument("--components", type=int, choices=(8, 16), required=True)
    result.add_argument("--widths", choices=("64,64", "128,128"), required=True)
    result.add_argument("--lambda-fp", type=float, choices=(0.01, 0.1), required=True)
    result.add_argument("--seed", type=int, choices=(4101, 4102, 4103), required=True)
    result.add_argument("--output", type=Path, required=True)
    result.add_argument("--pretrain-epochs", type=int, default=50)
    result.add_argument("--joint-epochs", type=int, default=450)
    result.add_argument("--steps-per-epoch", type=int, default=128)
    result.add_argument("--nll-batch", type=int, default=4096)
    result.add_argument("--fp-batch", type=int, default=512)
    result.add_argument("--gmm-sample", type=int, default=200_000)
    result.add_argument("--validation-sample", type=int, default=200_000)
    result.add_argument("--validation-fp-sample", type=int, default=20_000)
    return result


def main() -> None:
    args = parser().parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    started = time.time()
    tf.keras.backend.set_floatx("float64")
    tf.random.set_seed(args.seed)
    rng = np.random.default_rng(args.seed)
    selection_rng = np.random.default_rng(stable_case_seed(args.case))
    t1, t3 = case_parameters(args.case)

    train = load_split(args.case, "train")
    validation = load_split(args.case, "validation")
    if not np.isfinite(train).all() or np.any(train[:, :3] <= 0):
        raise RuntimeError("invalid training state")
    u = np.log(train[:, :3])
    u_mean = u.mean(axis=0)
    u_std = u.std(axis=0, ddof=1)
    z = (u-u_mean)/u_std
    gmm_n = min(args.gmm_sample, len(z))
    gmm_index = selection_rng.choice(len(z), size=gmm_n, replace=False)
    gmm = GaussianMixture(
        n_components=args.components, covariance_type="full", reg_covar=1e-5,
        max_iter=300, n_init=1, random_state=args.seed,
    ).fit(z[gmm_index])

    widths = tuple(int(item) for item in args.widths.split(","))
    model = NormalizedNLS3Density(
        args.components, widths, u_mean, u_std,
        gmm.weights_, gmm.means_, gmm.covariances_, args.seed,
    )
    _ = model.log_prob(tf.constant(train[:2], tf.float64))
    optimizer = tf.keras.optimizers.Adam(learning_rate=1e-3, global_clipnorm=10.0)

    lower = np.quantile(u, 0.001, axis=0)
    upper = np.quantile(u, 0.999, axis=0)
    val_n = min(args.validation_sample, len(validation))
    val_index = selection_rng.choice(len(validation), size=val_n, replace=False)
    val_state = validation[val_index]
    val_fp_n = min(args.validation_fp_sample, val_n)
    val_fp_index = selection_rng.choice(val_n, size=val_fp_n, replace=False)
    val_fp_state = val_state[val_fp_index]

    @tf.function(reduce_retracing=True)
    def train_step(nll_state: tf.Tensor, fp_state: tf.Tensor, weight: tf.Tensor):
        with tf.GradientTape() as tape:
            nll = -tf.reduce_mean(model.log_prob(nll_state))
            residual = relative_adjoint_residual(model, fp_state, t1, t3, 0.1)
            fp_loss = tf.reduce_mean(residual*residual)
            loss = nll+weight*fp_loss
        gradients = tape.gradient(loss, model.trainable_variables)
        optimizer.apply_gradients(zip(gradients, model.trainable_variables))
        return loss, nll, fp_loss

    def draw_fp_batch() -> np.ndarray:
        empirical_n = int(round(0.8*args.fp_batch))
        proposal_n = args.fp_batch-empirical_n
        empirical = train[rng.integers(0, len(train), size=empirical_n)]
        proposal = np.empty((proposal_n, 5), dtype=np.float64)
        proposal[:, :3] = np.exp(rng.uniform(lower, upper, size=(proposal_n, 3)))
        proposal[:, 3:] = rng.uniform(-np.pi, np.pi, size=(proposal_n, 2))
        return np.concatenate((empirical, proposal), axis=0)

    history = []
    best_nll = np.inf
    best_weights = None
    stale = 0
    current_lr = 1e-3
    total_epochs = args.pretrain_epochs+args.joint_epochs
    for epoch in range(total_epochs):
        physics_weight = 0.0 if epoch < args.pretrain_epochs else args.lambda_fp
        sums = np.zeros(3, dtype=np.float64)
        for _ in range(args.steps_per_epoch):
            nll_batch = train[rng.integers(0, len(train), size=args.nll_batch)]
            fp_batch = draw_fp_batch()
            losses = train_step(
                tf.constant(nll_batch, tf.float64), tf.constant(fp_batch, tf.float64),
                tf.constant(physics_weight, tf.float64),
            )
            sums += np.array([float(item) for item in losses])
        validation_nll = float(-np.mean(iter_log_prob(model, val_state)))
        history.append({
            "epoch": epoch+1,
            "phase": "pretrain" if physics_weight == 0 else "joint",
            "learning_rate": current_lr,
            "training_loss": sums[0]/args.steps_per_epoch,
            "training_nll": sums[1]/args.steps_per_epoch,
            "training_fp_loss": sums[2]/args.steps_per_epoch,
            "validation_nll": validation_nll,
        })
        if epoch+1 == args.pretrain_epochs:
            best_nll = np.inf; best_weights = None; stale = 0
        if epoch >= args.pretrain_epochs:
            if validation_nll < best_nll-1e-6:
                best_nll = validation_nll
                best_weights = [variable.numpy().copy() for variable in model.weights]
                stale = 0
            else:
                stale += 1
                if stale > 0 and stale % 15 == 0 and current_lr > 1e-5:
                    current_lr = max(1e-5, current_lr*0.5)
                    optimizer.learning_rate.assign(current_lr)
                if stale >= 30:
                    break

    if best_weights is None:
        best_weights = [variable.numpy().copy() for variable in model.weights]
        best_nll = float(-np.mean(iter_log_prob(model, val_state)))
    for variable, value in zip(model.weights, best_weights):
        variable.assign(value)

    val_logp = iter_log_prob(model, val_state)
    residual = iter_residual(model, val_fp_state, t1, t3)
    metrics = {
        "case": args.case,
        "components": args.components,
        "widths": list(widths),
        "lambda_fp": args.lambda_fp,
        "seed": args.seed,
        "training_rows": len(train),
        "validation_rows": len(validation),
        "validation_sample": val_n,
        "validation_fp_sample": val_fp_n,
        "best_validation_nll": best_nll,
        "fp_abs_median": float(np.median(np.abs(residual))),
        "fp_abs_p90": float(np.quantile(np.abs(residual), 0.90)),
        "fp_abs_p99": float(np.quantile(np.abs(residual), 0.99)),
        "epochs_completed": len(history),
        "runtime_seconds": time.time()-started,
    }
    if args.case.startswith("equilibrium"):
        exact_shape = -energy(val_state)/t1
        bulk_lo, bulk_hi = np.quantile(exact_shape, (0.005, 0.995))
        bulk = (exact_shape >= bulk_lo)&(exact_shape <= bulk_hi)
        slope, intercept = np.polyfit(exact_shape[bulk], val_logp[bulk], 1)
        additive_constant = float(np.mean(val_logp[bulk]-exact_shape[bulk]))
        centered = (val_logp[bulk]-exact_shape[bulk])-additive_constant
        metrics.update({
            "gibbs_bulk_log_slope": float(slope),
            "gibbs_bulk_log_intercept": float(intercept),
            "gibbs_validation_additive_constant": additive_constant,
            "gibbs_bulk_centered_log_rmse": float(np.sqrt(np.mean(centered*centered))),
            "gibbs_bulk_bounds": [float(bulk_lo), float(bulk_hi)],
            "gibbs_validation_percentiles": {
                str(q): float(np.quantile(exact_shape, q/100.0))
                for q in (1, 5, 25, 75, 95, 99)
            },
        })
    metrics["selection_admissible"] = bool(
        args.case == "equilibrium_dt2p5e-4"
        and 0.95 <= metrics.get("gibbs_bulk_log_slope", np.nan) <= 1.05
        and metrics["fp_abs_median"] <= 0.20
        and metrics["fp_abs_p90"] <= 1.00
    )

    with (args.output/"training_history.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(history[0]))
        writer.writeheader(); writer.writerows(history)
    (args.output/"validation_metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True)+"\n"
    )
    model.save_bundle(args.output, extra={
        "case": args.case, "T1": t1, "T3": t3,
        "lambda_fp": args.lambda_fp,
    })
    np.savez(args.output/"initial_gmm_and_scaler.npz",
             u_mean=u_mean, u_std=u_std, weights=gmm.weights_, means=gmm.means_,
             covariances=gmm.covariances_)
    print(json.dumps(metrics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
