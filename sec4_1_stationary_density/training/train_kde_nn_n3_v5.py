#!/usr/bin/env python3
"""v5 (2026-10-11) = previous trainer + two options, defaults reproduce the previous trainer exactly:
  --weak-set bounded : adds bounded angular stationarity test functions to the weak loss,
        psi = P_k sin(theta_k)/(c + P_k),  P_k cos(theta_k)/(c + P_k),  P_k = I_{k-1} I_k,  c in {0.25, 1}.
        They are smooth in Cartesian coordinates and weight low-action states, where the previous models
        underestimate the bond-angle asymmetry <sin theta_k> (diagnosis 2026-10-11: E[L psi] off by 3-5 SE).
        Like all weak terms they impose stationarity E_rho[L psi] = 0, not agreement with any data moment,
        so the untrained check moments stay independent.
  --val-layout streams : validation files are nls_cartesian dumps with 256 independent streams written in
        (checkpoint=4, batch=16, snapshot, lane=16) order (threads = 1); the check-moment standard error is
        then computed across stream means instead of 20 contiguous time blocks.
------------------------------------------------------------------------------
n = 3 port of train_kde_nn_v4.py (Section 4.2 pipeline), same losses, estimators and defaults.

    log rho(x) = -E(x)/5 + f(x) - log Z,   E = H/2,   f = N3GibbsRatio
Loss = data + pde_w * FP + weak_w * weak + mom_w * moments  (warm-up, ramp, cosine LR, as v4):
  data     smoothing-consistent fit of log S_h[rho_f](q) + c to the streamed KDE labels log v(q)
  FP       Huber loss (threshold 1/2) on the relative residual (L* rho)/rho at the query points
  weak     E_rho[L phi_k] = 0 for phi_k in {I_j, I_j^2, I1 I2 cos th2, I2 I3 cos th3, I1 I2 sin th2,
           I2 I3 sin th3}, self-normalised importance weights on the Gibbs (T = 5) Metropolis pool,
           unbiased two-half-batch product estimator
  moments  <I_j> and the bond currents J2, J3 from the labelling run's summary, same estimator
Validation (not trained on): held-out labels, FP residual on independent trajectory states, and 11
UNCONSTRAINED moments vs independent trajectories with block standard errors (check_z).
State: x = (I1, I2, I3, theta2, theta3), theta_j = 2(phi_j - phi_{j-1}) (paper: theta1 = -theta2).
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
import tensorflow as tf
from scipy.special import logsumexp

from n3_common import GAMMA, D, drift_tf, energy_np, half_covariance, read_kde, read_nls5d, read_summary, \
    relative_adjoint_residual
from n3_model import N3GibbsRatio

WEAK_NAMES = ["I1", "I2", "I3", "I1^2", "I2^2", "I3^2", "I1I2cos2", "I2I3cos3", "I1I2sin2", "I2I3sin3"]
OBS_NAMES = ["I1", "I2", "I3", "M", "E", "sin2", "cos2", "sin3", "cos3", "I1-I3", "I1I2sin2", "I2I3sin3"]
CHECK_NAMES = ["I1^2", "I2^2", "I3^2", "M^2", "I1I3", "cos2", "cos3", "sin2", "sin3", "I1I2cos2", "I2I3cos3"]


def _base_test_functions(x):
    i1, i2, i3, t2, t3 = tf.unstack(x, axis=-1)
    return tf.stack([i1, i2, i3, i1 * i1, i2 * i2, i3 * i3, i1 * i2 * tf.cos(t2), i2 * i3 * tf.cos(t3),
                     i1 * i2 * tf.sin(t2), i2 * i3 * tf.sin(t3)], axis=-1)


def observables_np(x):
    i1, i2, i3, t2, t3 = x.T
    return np.column_stack((i1, i2, i3, i1 + i2 + i3, energy_np(x), np.sin(t2), np.cos(t2), np.sin(t3),
                            np.cos(t3), i1 - i3, i1 * i2 * np.sin(t2), i2 * i3 * np.sin(t3)))


def moment_obs(x):
    return np.column_stack([x[:, :3], 4 * x[:, 0] * x[:, 1] * np.sin(x[:, 3]), 4 * x[:, 1] * x[:, 2] * np.sin(x[:, 4])])


def check_obs(x):  # NOT constrained: independent validation
    i1, i2, i3, t2, t3 = x.T
    m = i1 + i2 + i3
    return np.column_stack([i1 * i1, i2 * i2, i3 * i3, m * m, i1 * i3, np.cos(t2), np.cos(t3), np.sin(t2),
                            np.sin(t3), i1 * i2 * np.cos(t2), i2 * i3 * np.cos(t3)])


BOUNDED_C = (0.25, 1.0)
USE_BOUNDED = False


def enable_bounded():
    global USE_BOUNDED
    USE_BOUNDED = True
    for c in BOUNDED_C:
        for kind in ("sin", "cos"):
            for k in range(2, 3 + 1):
                WEAK_NAMES.append(f"P{k}{kind}{k}/(c{c}+P)")


def test_functions(x):
    base = _base_test_functions(x)
    if not USE_BOUNDED:
        return base
    cols = tf.unstack(x, axis=-1)
    acts, angs = cols[:3], cols[3:]
    extra = []
    for c in BOUNDED_C:
        for fn in (tf.sin, tf.cos):
            for k in range(3 - 1):
                p = acts[k] * acts[k + 1]
                extra.append(p * fn(angs[k]) / (c + p))
    return tf.concat([base, tf.stack(extra, axis=-1)], axis=-1)


def stream_means(x, obs):
    v = x.reshape(4, 16, -1, 16, 5)
    v = np.transpose(v, (1, 3, 0, 2, 4)).reshape(256, -1, 5)
    return np.array([obs(s).mean(0) for s in v])


def generator_values(x_np, t_left, t_right, bs=50_000):
    """(L phi_k)(x) = b . grad phi_k + sum_i a_ii d_ii phi_k, the same operator as the FP residual."""
    out, K = [], len(WEAK_NAMES)
    for j in range(0, len(x_np), bs):
        x = tf.constant(x_np[j:j + bs], tf.float64)
        with tf.GradientTape(persistent=True) as t2:
            t2.watch(x)
            with tf.GradientTape(persistent=True) as t1:
                t1.watch(x)
                phi = test_functions(x)
                cols = [tf.reduce_sum(phi[:, k]) for k in range(K)]
            grads = [t1.gradient(cols[k], x) for k in range(K)]
            gsum = [[tf.reduce_sum(grads[k][:, i]) for i in range(D)] for k in range(K)]
        hdiag = [tf.stack([t2.gradient(gsum[k][i], x)[:, i] for i in range(D)], -1) for k in range(K)]
        del t1, t2
        a = half_covariance(x, t_left, t_right)
        b = drift_tf(x, t_left, t_right)
        L = tf.stack([tf.reduce_sum(b * grads[k], -1) + tf.reduce_sum(a * hdiag[k], -1) for k in range(K)], -1)
        out.append(L.numpy())
    return np.concatenate(out)


def args_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--case", choices=("eq", "neq"), required=True)
    p.add_argument("--labels", type=Path, required=True)
    p.add_argument("--summary", type=Path)
    p.add_argument("--queries", type=Path, required=True)
    p.add_argument("--reference", type=Path, required=True, help="Gibbs T=5 Metropolis .bin (proposal pool)")
    p.add_argument("--validation", type=Path, nargs="+", required=True, help="independent trajectory .bin, same case")
    p.add_argument("--moment-weight", type=float, default=30.0)
    p.add_argument("--moment-score", type=float, default=3.0)
    p.add_argument("--jit", type=int, default=0)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--width", type=int, default=128)
    p.add_argument("--depth", type=int, default=3)
    p.add_argument("--epochs", type=int, default=600)
    p.add_argument("--steps", type=int, default=50)
    p.add_argument("--batch", type=int, default=256)
    p.add_argument("--kernel-mc", type=int, default=128)
    p.add_argument("--physics-batch", type=int, default=256)
    p.add_argument("--pde-weight", type=float, default=0.3)
    p.add_argument("--weak-weight", type=float, default=30.0)
    p.add_argument("--weak-set", choices=("base", "bounded"), default="base")
    p.add_argument("--val-layout", choices=("blocks", "streams"), default="blocks")
    p.add_argument("--weak-batch", type=int, default=32768)
    p.add_argument("--weak-pool", type=int, default=400_000)
    p.add_argument("--warmup", type=int, default=20)
    p.add_argument("--ramp", type=int, default=40)
    p.add_argument("--learning-rate", type=float, default=1e-3)
    p.add_argument("--lr-final-frac", type=float, default=0.02)
    p.add_argument("--holdout", type=float, default=0.1)
    p.add_argument("--min-weight", type=float, default=0.05)
    p.add_argument("--ref-eval", type=int, default=300_000)
    p.add_argument("--validate-every", type=int, default=10)
    p.add_argument("--weak-score", type=float, default=3.0)
    p.add_argument("--threads", type=int, default=0)
    return p


def main():
    a = args_parser().parse_args()
    if a.weak_set == "bounded":
        enable_bounded()
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
    T1, T3 = (5.0, 5.0) if a.case == "eq" else (2.0, 8.0)

    # ------------------------------------------------------------ labels
    q = read_nls5d(a.queries)
    kde, n_samples = read_kde(a.labels)
    summary = a.summary or Path(str(a.labels).replace("_kde.bin", "_summary.txt"))
    sv = read_summary(summary)
    hI, hT, rcut, reflect = float(sv["hI"]), float(sv["hT"]), float(sv["rcut"]), int(sv["reflect"])
    if float(sv["T1"]) != T1 or float(sv["T3"]) != T3:
        raise SystemExit("label temperatures do not match --case")
    target = np.array(sv["mean_I"] + sv["mean_current"])
    if not reflect:
        raise SystemExit("labels were computed without reflection")
    if len(q) != len(kde):
        raise SystemExit("queries and labels have different lengths")
    good = kde > 0
    q, y = q[good], np.log(kde[good])
    w = np.clip(np.exp(y - np.median(y)), a.min_weight, 1.0)
    perm = rng.permutation(len(q))
    n_hold = int(round(a.holdout * len(q)))
    hold, train = perm[:n_hold], perm[n_hold:]

    # ------------------------------------------------------------ trajectories
    reference = read_nls5d(a.reference)
    val_files = [read_nls5d(v) for v in a.validation]
    validation = np.concatenate(val_files)
    idx = rng.permutation(len(reference))
    n_eval = min(a.ref_eval, len(reference) // 2)
    ref_eval = reference[idx[:n_eval]]
    pool = reference[idx[n_eval:n_eval + a.weak_pool]]
    val_obs = observables_np(validation)
    val_mean, val_sd = val_obs.mean(0), val_obs.std(0, ddof=1)
    fp_eval = validation[rng.choice(len(validation), 2048, replace=False)]
    fp_scale = (1.0 + fp_eval[:, :3].sum(1)) ** 2
    action_scales = np.median(reference[:, :3], axis=0)

    blocks = []
    for vf in val_files:
        if a.val_layout == "streams":
            blocks += list(stream_means(vf, check_obs))
        else:
            blocks += [check_obs(c).mean(0) for c in np.array_split(vf, 20)]
    blocks = np.array(blocks)
    check_truth = check_obs(validation).mean(0)
    check_se = blocks.std(0, ddof=1) / np.sqrt(len(blocks))
    M_pool = moment_obs(pool) - target
    M_scale = M_pool.std(0)
    M_pool_tf = tf.constant(M_pool / M_scale)
    O_eval, C_eval = moment_obs(ref_eval), check_obs(ref_eval)

    t0 = time.time()
    G_pool = generator_values(pool, T1, T3)
    G_scale = G_pool.std(0)
    G_pool /= G_scale
    G_eval = generator_values(ref_eval, T1, T3) / G_scale
    G_val = generator_values(validation[rng.choice(len(validation), min(300_000, len(validation)), replace=False)],
                             T1, T3) / G_scale
    weak_truth = G_val.mean(0)
    print(json.dumps({"weak_setup_seconds": time.time() - t0, "weak_truth_rms": float(np.sqrt(np.mean(weak_truth ** 2))),
                      "weak_truth": dict(zip(WEAK_NAMES, np.round(weak_truth, 4).tolist()))}), flush=True)
    G_pool_tf = tf.constant(G_pool)

    def kernel_offsets(m, gen):
        out = np.empty((0, D))
        while len(out) < m:
            z = gen.standard_normal((2 * m, D))
            out = np.vstack([out, z[np.sum(z * z, 1) <= rcut * rcut]])
        return out[:m] * np.array([hI] * 3 + [hT] * 2)

    z_eval = tf.constant(kernel_offsets(4 * a.kernel_mc, np.random.default_rng(a.seed + 99)))

    # ------------------------------------------------------------ model
    model = N3GibbsRatio(action_scales, a.seed, width=a.width, depth=a.depth)
    _ = model(tf.constant(q[:2]))
    total_steps = a.epochs * a.steps
    lr = tf.keras.optimizers.schedules.CosineDecay(a.learning_rate, total_steps, alpha=a.lr_final_frac)
    opt = tf.keras.optimizers.Adam(lr, global_clipnorm=10.0)
    pi = tf.constant(math.pi, tf.float64)

    def log_smoothed(xq, z):
        b, m = tf.shape(xq)[0], tf.shape(z)[0]
        pts = xq[:, None, :] + z[None, :, :]
        acts = tf.abs(pts[..., :3])
        angs = tf.math.floormod(pts[..., 3:] + pi, 2.0 * pi) - pi
        flat = tf.reshape(tf.concat([acts, angs], -1), (-1, D))
        lp = tf.reshape(model.log_prob(flat), (b, m))
        return tf.reduce_logsumexp(lp, axis=1)

    def data_loss(xq, yq, wq, z):
        r = log_smoothed(xq, z) - yq
        c = tf.reduce_sum(wq * r) / tf.reduce_sum(wq)
        return tf.reduce_sum(wq * tf.square(r - c)) / tf.reduce_sum(wq)

    def weak_mean(xw, gw):
        lw = tf.nn.log_softmax(model.log_ratio(xw))
        return tf.reduce_sum(tf.exp(lw)[:, None] * gw, 0)

    def make_step(use_pde, use_wm):
        @tf.function(reduce_retracing=True, jit_compile=bool(a.jit))
        def step(xq, yq, wq, z, xfp, pde_w, xa, ga, xb, gb, weak_w, ma, mb, mom_w):
            zero = tf.constant(0.0, tf.float64)
            with tf.GradientTape() as tape:
                ld = data_loss(xq, yq, wq, z)
                lp, lwk, lmo = zero, zero, zero
                if use_pde:
                    rr = relative_adjoint_residual(model, xfp, T1, T3)
                    ar = tf.abs(rr)
                    lp = tf.reduce_mean(tf.where(ar <= 0.5, 0.5 * rr * rr, 0.5 * (ar - 0.25)))
                if use_wm:
                    ra = weak_mean(xa, tf.concat([ga, ma], 1))
                    rb = weak_mean(xb, tf.concat([gb, mb], 1))
                    prod = ra * rb
                    nk = ga.shape[1]
                    lwk = tf.reduce_sum(prod[:nk])
                    lmo = tf.reduce_sum(prod[nk:])
                reg = tf.add_n([tf.reduce_sum(v * v) for v in model.trainable_variables])
                loss = ld + pde_w * lp + weak_w * lwk + mom_w * lmo + 1e-7 * reg
            g = tape.gradient(loss, model.trainable_variables)
            opt.apply_gradients(zip(g, model.trainable_variables))
            return ld, lp, lwk, lmo
        return step

    steps_cache = {}
    logratio_c = tf.function(model.log_ratio, reduce_retracing=True)
    smooth_c = tf.function(lambda t: log_smoothed(t, z_eval), reduce_retracing=True)
    resid_c = tf.function(lambda t: relative_adjoint_residual(model, t, T1, T3), reduce_retracing=True)

    def batched(fn, x, bs=8192):
        return np.concatenate([fn(tf.constant(x[j:j + bs])).numpy() for j in range(0, len(x), bs)])

    def validate(epoch):
        out = {"epoch": epoch}
        ls = batched(smooth_c, q[hold], 256)
        r = ls - y[hold]; wh = w[hold]
        c = np.sum(wh * r) / np.sum(wh)
        out["heldout_label_rmse_w"] = float(np.sqrt(np.sum(wh * (r - c) ** 2) / np.sum(wh)))
        out["heldout_label_rmse"] = float(np.sqrt(np.mean((r - c) ** 2)))
        rr = np.concatenate([resid_c(tf.constant(fp_eval[j:j + 256])).numpy() for j in range(0, len(fp_eval), 256)])
        out["fp_abs_median"] = float(np.median(np.abs(rr)))
        out["fp_abs_p90"] = float(np.quantile(np.abs(rr), 0.9))
        out["fp_scaled_median"] = float(np.median(np.abs(rr) / fp_scale))
        out["fp_scaled_p90"] = float(np.quantile(np.abs(rr) / fp_scale, 0.9))
        f = batched(logratio_c, ref_eval)
        wt = np.exp(f - logsumexp(f))
        out["ess_fraction"] = float(1.0 / np.sum(wt * wt) / len(wt))
        mom = np.sum(wt[:, None] * observables_np(ref_eval), 0)
        zz = (mom - val_mean) / val_sd
        out["normalized_moment_rmse"] = float(np.sqrt(np.mean(zz * zz)))
        out["worst_moment"] = OBS_NAMES[int(np.argmax(np.abs(zz)))]
        out["model_mean_I"] = mom[:3].tolist()
        out["model_currents"] = (4 * mom[10:12]).tolist()
        wk = wt @ G_eval
        out["weak_rms"] = float(np.sqrt(np.mean(wk * wk)))
        out["weak_worst"] = WEAK_NAMES[int(np.argmax(np.abs(wk)))]
        out["weak_dI_dt"] = (wk[:3] * G_scale[:3]).tolist()
        md = wt @ O_eval - target
        out["moment_dev"] = md.tolist()
        out["moment_rms"] = float(np.sqrt(np.mean((md / M_scale) ** 2)))
        cz = (wt @ C_eval - check_truth) / check_se
        out["check_z"] = dict(zip(CHECK_NAMES, np.round(cz, 2).tolist()))
        out["check_z_rms"] = float(np.sqrt(np.mean(cz * cz)))
        if a.case == "eq":
            fv = batched(logratio_c, validation[:20000])
            fv = fv - fv.mean()
            out["eq_known_answer_f_rms"] = float(np.sqrt(np.mean(fv * fv)))
        return out

    history = [validate(0)]
    print(json.dumps(history[-1]), flush=True)
    best, best_w = None, None
    nw = a.weak_batch // 2
    for epoch in range(1, a.epochs + 1):
        frac = 0.0 if epoch <= a.warmup else min(1.0, (epoch - a.warmup) / max(1, a.ramp))
        pde_w = tf.constant(a.pde_weight * frac, tf.float64)
        weak_w = tf.constant(a.weak_weight * frac, tf.float64)
        mom_w = tf.constant(a.moment_weight * frac, tf.float64)
        sums = np.zeros(4)
        for _ in range(a.steps):
            bi = train[rng.integers(0, len(train), a.batch)]
            fi = rng.integers(0, len(q), a.physics_batch)
            wi = rng.integers(0, len(pool), 2 * nw)
            key = (float(pde_w.numpy()) > 0, float(weak_w.numpy()) > 0 or float(mom_w.numpy()) > 0)
            if key not in steps_cache:
                steps_cache[key] = make_step(*key)
            vals = steps_cache[key](tf.constant(q[bi]), tf.constant(y[bi]), tf.constant(w[bi]),
                                    tf.constant(kernel_offsets(a.kernel_mc, rng)), tf.constant(q[fi]), pde_w,
                                    tf.constant(pool[wi[:nw]]), tf.gather(G_pool_tf, wi[:nw]),
                                    tf.constant(pool[wi[nw:]]), tf.gather(G_pool_tf, wi[nw:]), weak_w,
                                    tf.gather(M_pool_tf, wi[:nw]), tf.gather(M_pool_tf, wi[nw:]), mom_w)
            sums += [float(v.numpy()) for v in vals]
        row = {"epoch": epoch, "pde_weight": float(pde_w.numpy()), "weak_weight": float(weak_w.numpy()),
               "train_data": sums[0] / a.steps, "train_physics": sums[1] / a.steps,
               "train_weak": sums[2] / a.steps, "train_moment": sums[3] / a.steps,
               "moment_weight": float(mom_w.numpy()), "elapsed": time.time() - t_start}
        if epoch % a.validate_every == 0 or epoch == a.epochs:
            tv = time.time(); row.update(validate(epoch)); row["validate_seconds"] = time.time() - tv
            score = (row["heldout_label_rmse_w"] + 0.1 * row["fp_abs_median"]
                     + a.weak_score * row["weak_rms"] + a.moment_score * row["moment_rms"])
            row["score"] = score
            if epoch > a.warmup + a.ramp and (best is None or score < best):
                best, best_w = score, [v.numpy() for v in model.trainable_variables]
                row["selected"] = True
            print(json.dumps(row), flush=True)
        history.append(row)

    if best_w is not None:
        for v, val in zip(model.trainable_variables, best_w):
            v.assign(val)
    final = validate("selected")
    last = next(r for r in reversed(history) if "weak_rms" in r)
    model.save_bundle(a.output, {"case": a.case, "labels": str(a.labels), "queries": str(a.queries),
                                 "label_samples": n_samples, "hI": hI, "hT": hT, "rcut": rcut,
                                 "script": "train_kde_nn_n3.py"})
    (a.output / "history.json").write_text(json.dumps(history, indent=1) + "\n")
    result = {"args": {k: str(v) for k, v in vars(a).items()}, "final": final, "last_epoch": last,
              "validation_trajectory_mean_I": val_mean[:3].tolist(),
              "validation_trajectory_currents": (4 * val_mean[10:12]).tolist(),
              "weak_truth_rms": float(np.sqrt(np.mean(weak_truth ** 2))),
              "moment_target": target.tolist(), "moment_scale": M_scale.tolist(),
              "check_truth": dict(zip(CHECK_NAMES, check_truth.tolist())),
              "check_se": dict(zip(CHECK_NAMES, check_se.tolist())),
              "weak_names": WEAK_NAMES, "weak_scale": G_scale.tolist(),
              "n_labels_train": int(len(train)), "n_labels_holdout": int(n_hold),
              "runtime_seconds": time.time() - t_start}
    (a.output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
