#!/usr/bin/env python3
"""Where does the n=4 NESS network get <sin theta_j> wrong, and do candidate weak test functions see it?

1. <sin theta_j> conditioned on quantile bins of P_j = I_{j-1} I_j: network (IS on Gibbs) vs long trajectories.
2. Weak residual E[L psi] for bounded test functions psi = P sin(theta)/(c+P), P cos(theta)/(c+P), c = 0.25, 1:
   network (IS) vs long trajectories (truth 0, SE across streams).
usage: diag_sin_n4.py RUN_DIR [RUN_DIR ...]   (RES dir with n4_model.py, gibbs_n4_T5.bin given by env RES)
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import tensorflow as tf

tf.keras.backend.set_floatx("float64")
RES = Path(os.environ["RES"]); VAL = Path(os.environ["VAL"])
sys.path.insert(0, str(RES)); sys.path.insert(0, str(Path(__file__).parent))
from n4_common import drift_tf, GAMMA  # noqa: E402
from n4_model import N4GibbsRatio  # noqa: E402
from eval_check_long import read_states  # noqa: E402

T1, T4 = 2.0, 8.0
CS = (0.25, 1.0)
NAMES = [f"{k}{j}_c{c}" for c in CS for k in ("Psin", "Pcos") for j in (2, 3, 4)]


def psi(x):
    i = [x[:, k] for k in range(4)]; t = [x[:, 4 + k] for k in range(3)]
    cols = []
    for c in CS:
        for f in (tf.sin, tf.cos):
            for j in range(3):
                P = i[j] * i[j + 1]
                cols.append(P * f(t[j]) / (c + P))
    return tf.stack(cols, -1)


def generator(x_np, bs=20000):
    out, g = [], GAMMA
    K = len(NAMES)
    for s in range(0, len(x_np), bs):
        x = tf.constant(x_np[s:s + bs], tf.float64)
        with tf.GradientTape(persistent=True) as t2:
            t2.watch(x)
            with tf.GradientTape(persistent=True) as t1:
                t1.watch(x)
                ph = psi(x)
                cols = [tf.reduce_sum(ph[:, k]) for k in range(K)]
            grads = [t1.gradient(cols[k], x) for k in range(K)]
            gsum = [[tf.reduce_sum(grads[k][:, i]) for i in range(7)] for k in range(K)]
        hd = [tf.stack([t2.gradient(gsum[k][i], x)[:, i] for i in range(7)], -1) for k in range(K)]
        del t1, t2
        i1, i4 = x[:, 0], x[:, 3]; z = tf.zeros_like(i1)
        a = tf.stack((4 * g * T1 * i1, z, z, 4 * g * T4 * i4, 4 * g * T1 / i1, z, 4 * g * T4 / i4), -1)
        b = drift_tf(x, T1, T4)
        out.append(tf.stack([tf.reduce_sum(b * grads[k], -1) + tf.reduce_sum(a * hd[k], -1) for k in range(K)], -1).numpy())
    return np.concatenate(out)


def streams(path, thin):
    x = read_states(path).reshape(4, 16, -1, 16, 7)
    x = np.transpose(x, (1, 3, 0, 2, 4)).reshape(256, -1, 7)
    return x[:, ::thin]


# ---- reference (long trajectories, both sets)
V = np.concatenate([streams(VAL / "n4_neq_valA.bin", 1), streams(VAL / "n4_neq_valB.bin", 1)])   # 512 streams
Vt = V[:, ::25]                                                                                    # for L psi
LV = generator(Vt.reshape(-1, 7)).reshape(Vt.shape[0], Vt.shape[1], -1).mean(1)
lv_mean, lv_se = LV.mean(0), LV.std(0, ddof=1) / np.sqrt(len(LV))
Vflat = V.reshape(-1, 7)
ref = read_states(RES / "gibbs_n4_T5.bin")
Lref = generator(ref)
out = {"names": NAMES, "traj_L_mean": lv_mean.tolist(), "traj_L_se": lv_se.tolist()}
print("trajectories: E[L psi] / SE:", np.round(lv_mean / lv_se, 2).tolist(), flush=True)

for rd in sys.argv[1:]:
    rd = Path(rd)
    cfg = json.loads((rd / "model_config.json").read_text())
    m = N4GibbsRatio(np.array(cfg["action_scales"]), cfg["seed"], width=cfg["width"],
                     feature_version=cfg["feature_version"], depth=cfg["depth"])
    _ = m(tf.constant(ref[:2])); m.load_weights(str(rd / "model.weights.h5"))
    f = np.concatenate([m.log_ratio(tf.constant(ref[s:s + 20000])).numpy().ravel() for s in range(0, len(ref), 20000)])
    w = np.exp(f - f.max()); w /= w.sum()
    Lm = w @ Lref
    blocks = [(wb @ lb) / wb.sum() for wb, lb in zip(np.array_split(w, 50), np.array_split(Lref, 50))]
    lm_se = np.std(blocks, 0, ddof=1) / np.sqrt(50)
    z = (Lm - lv_mean) / np.hypot(lm_se, lv_se)
    rec = {"weak_L": Lm.tolist(), "weak_L_se": lm_se.tolist(), "z": dict(zip(NAMES, np.round(z, 2).tolist()))}
    # conditional sin theta
    cond = {}
    for j in range(3):
        Pv = Vflat[:, j] * Vflat[:, j + 1]; Pr = ref[:, j] * ref[:, j + 1]
        qs = np.quantile(Pv, [0, .2, .4, .6, .8, 1.0]); qs[-1] = np.inf
        rows = []
        for lo, hi in zip(qs[:-1], qs[1:]):
            sv = (Pv >= lo) & (Pv < hi); sr = (Pr >= lo) & (Pr < hi)
            tv = np.sin(Vflat[sv, 4 + j]).mean()
            tm = (w[sr] @ np.sin(ref[sr, 4 + j])) / w[sr].sum()
            mass_v = sv.mean(); mass_m = w[sr].sum()
            rows.append([float(lo), float(tv), float(tm), float(mass_v), float(mass_m)])
        cond[f"theta{j+2}"] = rows
    rec["cond_sin_by_P_quintile"] = cond
    out[str(rd)] = rec
    print("==", rd.name, "weak z:", rec["z"], flush=True)
    for k, rows in cond.items():
        print("  ", k, " | ".join("P>=%.2f: traj %.4f net %.4f (mass %.3f/%.3f)" % tuple(r) for r in rows), flush=True)
Path("diag_sin_n4.json").write_text(json.dumps(out, indent=1))
