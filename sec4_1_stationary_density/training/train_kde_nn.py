#!/usr/bin/env python3
"""Fit the n=4 Gibbs-anchored network to streamed KDE labels plus the FP equation.

Model (unchanged from n4_7d_neural_fp_2026-09-25/scripts/n4_model.py):
    log rho(x) = -E(x)/5 + f(x) - log Z,   E = H/2,
f = N4GibbsRatio network.  The exact T=5 Gibbs density is f = 0.

Data term ("smoothing-consistent"): the streamed KDE at query q estimates
    E[KDE(q)] = int rho(|y_I|, y_theta) phi_h(q - y) 1{|q-y|_h <= r} dy,
i.e. the reflected, truncated Gaussian smoothing of rho.  We apply the same
smoothing to the MODEL (Monte Carlo over M kernel offsets) and fit
    log S_h[rho_f](q) + c  ~  log KDE(q),
with c the (weighted) mean offset.  The bandwidth bias then appears on both
sides and cancels instead of being learned.

Physics term: Huber loss on the relative stationary FP residual (L* rho)/rho,
the audited 7D operator of n4_common.relative_adjoint_residual.

Usage:
  python3 train_kde_nn.py --case eq  --labels eq_h0.2_..._kde.bin  --queries queries_eq.bin  ...
  python3 train_kde_nn.py --case neq --labels neq_h0.2_..._kde.bin --queries queries_neq.bin ...
See --help.  Writes OUT/{model_config.json, model.weights.h5, history.json, result.json}.
"""
from __future__ import annotations

import argparse
import json
import math
import struct
import time
from pathlib import Path

import numpy as np
import tensorflow as tf
from scipy.special import logsumexp

from n4_common import energy_np, relative_adjoint_residual
from n4_model import N4GibbsRatio


# ------------------------------------------------------------------ I/O
def read_nls7d(path) -> np.ndarray:
    with open(path, "rb") as f:
        if f.read(8) != b"NLS7D01\0":
            raise ValueError(f"bad header: {path}")
        n, d, b = struct.unpack("<QII", f.read(16))
        x = np.fromfile(f, "<f4")
    if d != 7 or b != 4 or x.size != n * 7:
        raise ValueError(f"bad layout: {path}")
    return x.reshape(n, 7).astype(np.float64)


def read_kde(path):
    with open(path, "rb") as f:
        if f.read(8) != b"NLS7KDE1":
            raise ValueError(f"bad header: {path}")
        q, n = struct.unpack("<QQ", f.read(16))
        d = np.fromfile(f, "<f8")
    if d.size != q:
        raise ValueError(f"truncated: {path}")
    return d, int(n)


def read_summary(path):
    """hI, hT, rcut, reflect from the streamed run's summary file."""
    vals = {}
    for line in Path(path).read_text().splitlines():
        toks = line.split()
        for k, v in zip(toks[::2], toks[1::2]):
            vals[k] = v
    return float(vals["hI"]), float(vals["hT"]), float(vals["rcut"]), int(vals["reflect"])


def observables_np(x):
    i1, i2, i3, i4, t2, t3, t4 = x.T
    return np.column_stack((
        i1, i2, i3, i4, i1 + i2 + i3 + i4, energy_np(x),
        np.sin(t2), np.cos(t2), np.sin(t3), np.cos(t3), np.sin(t4), np.cos(t4),
        i1 - i4, i1 * i2 * np.sin(t2), i2 * i3 * np.sin(t3), i3 * i4 * np.sin(t4)))


OBS_NAMES = ["I1", "I2", "I3", "I4", "M", "E", "sin2", "cos2", "sin3", "cos3",
             "sin4", "cos4", "I1-I4", "I1I2sin2", "I2I3sin3", "I3I4sin4"]


def args_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--case", choices=("eq", "neq"), required=True)
    p.add_argument("--labels", type=Path, required=True, help="streamed *_kde.bin")
    p.add_argument("--summary", type=Path, help="streamed *_summary.txt (default: next to labels)")
    p.add_argument("--queries", type=Path, required=True, help="queries_*.bin used by that run")
    p.add_argument("--reference", type=Path, required=True,
                   help="equilibrium (5,5) trajectory or Gibbs .bin: proposal for ESS/moments")
    p.add_argument("--validation", type=Path, required=True,
                   help="independent trajectory .bin of the SAME case (moments, FP check); "
                        "for --case eq a direct-Gibbs .bin is also used for the known answer")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--feature-version", default="phase_regular_v4",
                   choices=("regular_v1", "transport_v2", "transport_v3", "phase_regular_v4"))
    p.add_argument("--width", type=int, default=128)
    p.add_argument("--depth", type=int, default=2)
    p.add_argument("--epochs", type=int, default=120)
    p.add_argument("--steps", type=int, default=50, help="optimizer steps per epoch")
    p.add_argument("--batch", type=int, default=256, help="labelled queries per step")
    p.add_argument("--kernel-mc", type=int, default=128, help="kernel offsets M per query")
    p.add_argument("--physics-batch", type=int, default=128)
    p.add_argument("--pde-weight", type=float, default=1.0)
    p.add_argument("--warmup", type=int, default=20, help="data-only epochs before FP ramps in")
    p.add_argument("--ramp", type=int, default=30, help="epochs to ramp FP weight 0 -> pde-weight")
    p.add_argument("--learning-rate", type=float, default=1e-3)
    p.add_argument("--holdout", type=float, default=0.1, help="fraction of labels held out")
    p.add_argument("--min-weight", type=float, default=0.05)
    p.add_argument("--threads", type=int, default=0, help="TensorFlow CPU threads (0 = all)")
    return p


def main():
    a = args_parser().parse_args()
    if a.output.exists() and any(a.output.iterdir()):
        raise SystemExit(f"refusing to overwrite {a.output}")
    a.output.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    if a.threads > 0:
        tf.config.threading.set_intra_op_parallelism_threads(a.threads)
        tf.config.threading.set_inter_op_parallelism_threads(2)
    tf.keras.backend.set_floatx("float64")
    tf.random.set_seed(a.seed)
    rng = np.random.default_rng(a.seed)
    T1, T4 = (5.0, 5.0) if a.case == "eq" else (2.0, 8.0)

    # ------------------------------------------------------------ data
    q = read_nls7d(a.queries)
    kde, n_samples = read_kde(a.labels)
    summary = a.summary or Path(str(a.labels).replace("_kde.bin", "_summary.txt"))
    hI, hT, rcut, reflect = read_summary(summary)
    if not reflect:
        raise SystemExit("labels were computed without reflection; this script assumes --reflect 1")
    if len(q) != len(kde):
        raise SystemExit("queries and labels have different lengths")
    good = kde > 0
    q, y = q[good], np.log(kde[good])
    w = np.clip(np.exp(y - np.median(y)), a.min_weight, 1.0)  # ~ inverse label variance
    perm = rng.permutation(len(q))
    n_hold = int(round(a.holdout * len(q)))
    hold, train = perm[:n_hold], perm[n_hold:]

    reference = read_nls7d(a.reference)
    validation = read_nls7d(a.validation)
    ref_eval = reference[rng.choice(len(reference), min(100_000, len(reference)), replace=False)]
    val_obs = observables_np(validation)
    val_mean, val_sd = val_obs.mean(0), val_obs.std(0, ddof=1)
    fp_eval = validation[rng.choice(len(validation), 1024, replace=False)]
    action_scales = np.median(reference[:, :4], axis=0)

    # Fixed kernel offsets for evaluation; fresh ones each training step.
    def kernel_offsets(m, gen):
        out = np.empty((0, 7))
        while len(out) < m:
            z = gen.standard_normal((2 * m, 7))
            out = np.vstack([out, z[np.sum(z * z, 1) <= rcut * rcut]])
        return out[:m] * np.array([hI] * 4 + [hT] * 3)

    z_eval = tf.constant(kernel_offsets(4 * a.kernel_mc, np.random.default_rng(a.seed + 99)))

    # ------------------------------------------------------------ model
    model = N4GibbsRatio(action_scales, a.seed, width=a.width,
                         feature_version=a.feature_version, depth=a.depth)
    _ = model(tf.constant(q[:2]))
    opt = tf.keras.optimizers.Adam(a.learning_rate, global_clipnorm=10.0)
    pi = tf.constant(math.pi, tf.float64)

    def log_smoothed(xq, z):
        """log of (reflected, truncated) kernel smoothing of exp(log_prob), up to a constant."""
        b, m = tf.shape(xq)[0], tf.shape(z)[0]
        pts = xq[:, None, :] + z[None, :, :]
        acts = tf.abs(pts[..., :4])
        angs = tf.math.floormod(pts[..., 4:] + pi, 2.0 * pi) - pi
        flat = tf.reshape(tf.concat([acts, angs], -1), (-1, 7))
        lp = tf.reshape(model.log_prob(flat), (b, m))
        return tf.reduce_logsumexp(lp, axis=1)

    def data_loss(xq, yq, wq, z):
        r = log_smoothed(xq, z) - yq
        c = tf.reduce_sum(wq * r) / tf.reduce_sum(wq)
        return tf.reduce_sum(wq * tf.square(r - c)) / tf.reduce_sum(wq)

    @tf.function(reduce_retracing=True)
    def step(xq, yq, wq, z, xfp, pde_w):
        with tf.GradientTape() as tape:
            ld = data_loss(xq, yq, wq, z)
            if pde_w > 0:
                rr = relative_adjoint_residual(model, xfp, T1, T4)
                ar = tf.abs(rr)
                lp = tf.reduce_mean(tf.where(ar <= 0.5, 0.5 * rr * rr, 0.5 * (ar - 0.25)))
            else:
                lp = tf.constant(0.0, tf.float64)
            reg = tf.add_n([tf.reduce_sum(v * v) for v in model.trainable_variables])
            loss = ld + pde_w * lp + 1e-7 * reg
        g = tape.gradient(loss, model.trainable_variables)
        opt.apply_gradients(zip(g, model.trainable_variables))
        return ld, lp

    def batched(fn, x, bs=4096):
        return np.concatenate([fn(tf.constant(x[j:j + bs])).numpy() for j in range(0, len(x), bs)])

    def validate(epoch):
        out = {"epoch": epoch}
        # held-out label fit (same smoothing, fixed offsets)
        ls = batched(lambda t: log_smoothed(t, z_eval), q[hold], 256)
        r = ls - y[hold]; wh = w[hold]
        c = np.sum(wh * r) / np.sum(wh)
        out["heldout_label_rmse_w"] = float(np.sqrt(np.sum(wh * (r - c) ** 2) / np.sum(wh)))
        out["heldout_label_rmse"] = float(np.sqrt(np.mean((r - c) ** 2)))
        # FP residual on independent trajectory states
        rr = np.concatenate([relative_adjoint_residual(model, tf.constant(fp_eval[j:j + 256]), T1, T4).numpy()
                             for j in range(0, len(fp_eval), 256)])
        out["fp_abs_median"] = float(np.median(np.abs(rr)))
        out["fp_abs_p90"] = float(np.quantile(np.abs(rr), 0.9))
        # importance weights from the equilibrium proposal -> ESS and moments
        f = batched(model.log_ratio, ref_eval)
        lw = f - logsumexp(f)
        wt = np.exp(lw)
        out["ess_fraction"] = float(1.0 / np.sum(wt * wt) / len(wt))
        mom = np.sum(wt[:, None] * observables_np(ref_eval), 0)
        z = (mom - val_mean) / val_sd
        out["normalized_moment_rmse"] = float(np.sqrt(np.mean(z * z)))
        out["worst_moment"] = OBS_NAMES[int(np.argmax(np.abs(z)))]
        out["model_mean_I"] = mom[:4].tolist()
        out["model_currents"] = (4 * mom[13:16]).tolist()
        if a.case == "eq":
            fv = batched(model.log_ratio, validation[:20000])
            fv = fv - fv.mean()
            out["eq_known_answer_f_rms"] = float(np.sqrt(np.mean(fv * fv)))  # exact answer: f = 0
        return out

    history = [validate(0)]
    print(json.dumps(history[-1]), flush=True)
    best, best_w = None, None
    for epoch in range(1, a.epochs + 1):
        frac = 0.0 if epoch <= a.warmup else min(1.0, (epoch - a.warmup) / max(1, a.ramp))
        pde_w = tf.constant(a.pde_weight * frac, tf.float64)
        sums = np.zeros(2)
        for _ in range(a.steps):
            bi = train[rng.integers(0, len(train), a.batch)]
            fi = rng.integers(0, len(q), a.physics_batch)
            vals = step(tf.constant(q[bi]), tf.constant(y[bi]), tf.constant(w[bi]),
                        tf.constant(kernel_offsets(a.kernel_mc, rng)),
                        tf.constant(q[fi]), pde_w)
            sums += [float(v.numpy()) for v in vals]
        row = {"epoch": epoch, "pde_weight": float(pde_w.numpy()),
               "train_data": sums[0] / a.steps, "train_physics": sums[1] / a.steps}
        if epoch % 5 == 0 or epoch == a.epochs:
            row.update(validate(epoch))
            # selection uses held-out labels and the FP residual only (no trajectory moments)
            score = row["heldout_label_rmse_w"] + 0.1 * row["fp_abs_median"]
            if epoch > a.warmup and (best is None or score < best):
                best, best_w = score, [v.numpy() for v in model.trainable_variables]
                row["selected"] = True
            print(json.dumps(row), flush=True)
        history.append(row)

    if best_w is not None:
        for v, val in zip(model.trainable_variables, best_w):
            v.assign(val)
    final = validate("selected")
    model.save_bundle(a.output, {"case": a.case, "labels": str(a.labels), "queries": str(a.queries),
                                 "label_samples": n_samples, "hI": hI, "hT": hT, "rcut": rcut})
    (a.output / "history.json").write_text(json.dumps(history, indent=1) + "\n")
    result = {"args": {k: str(v) for k, v in vars(a).items()}, "final": final,
              "validation_trajectory_mean_I": val_mean[:4].tolist(),
              "validation_trajectory_currents": (4 * val_mean[13:16]).tolist(),
              "n_labels_train": int(len(train)), "n_labels_holdout": int(n_hold),
              "runtime_seconds": time.time() - t_start}
    (a.output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
