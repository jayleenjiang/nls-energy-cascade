#!/usr/bin/env python3
"""Re-evaluate the untrained-moment check (check_z) of trained density networks against long,
independent validation trajectories.  No retraining: model selection never used check_z.

Validation reference: nls_cartesian runs with 256 independent streams (16 batches x 16 lanes,
threads = 1 so the dump order is (checkpoint, batch, snapshot, lane)), burn-in 500, measurement 10000
sampled every 0.5.  The standard error is computed across the 256 per-stream means, so time
correlation inside a stream does not enter.  Two independent sets (A, B) give the noise reference.

Network moments: self-normalised importance sampling with weights exp(f) on the full Gibbs (T = 5)
reference; their standard error from 50 contiguous blocks of the reference (ratio estimator per block).

z = (network - validation) / sqrt(se_val^2 + se_IS^2).  Also reported: the old short validation
trajectories against the long ones, and (equilibrium) everything against the exact Gibbs moments.

usage: eval_check_long.py RESULTS_DIR VAL_DIR OUT.json
"""
import json
import struct
import sys
from pathlib import Path

import numpy as np
import tensorflow as tf

tf.keras.backend.set_floatx("float64")


def read_states(path):
    with open(path, "rb") as f:
        magic = f.read(8)
        n, d, b = struct.unpack("<QII", f.read(16))
        x = np.fromfile(f, "<f4")
    if not magic.startswith(b"NLS") or b != 4 or x.size != n * d:
        raise ValueError(path)
    return x.reshape(n, d).astype(np.float64)


def check_obs(x, n):
    I = x[:, :n]
    th = x[:, n:]
    m = I.sum(1)
    cols = [I[:, j] ** 2 for j in range(n)] + [m * m, I[:, 0] * I[:, n - 1]]
    cols += [np.cos(th[:, k]) for k in range(n - 1)] + [np.sin(th[:, k]) for k in range(n - 1)]
    cols += [I[:, k] * I[:, k + 1] * np.cos(th[:, k]) for k in range(n - 1)]
    return np.column_stack(cols)


def check_names(n):
    return ([f"I{j}^2" for j in range(1, n + 1)] + ["M^2", f"I1I{n}"] + [f"cos{k}" for k in range(2, n + 1)]
            + [f"sin{k}" for k in range(2, n + 1)] + [f"I{k-1}I{k}cos{k}" for k in range(2, n + 1)])


def extra_obs(x, n):  # profile and currents, for information
    I, th = x[:, :n], x[:, n:]
    return np.column_stack([I] + [4 * I[:, k] * I[:, k + 1] * np.sin(th[:, k]) for k in range(n - 1)])


def stream_means(x, n, obs):
    D = 2 * n - 1
    v = x.reshape(4, 16, -1, 16, D)               # checkpoint, batch, snapshot, lane, coordinate
    v = np.transpose(v, (1, 3, 0, 2, 4)).reshape(256, -1, D)
    return np.array([obs(s, n).mean(0) for s in v])


def val_stats(path, n, obs):
    sm = stream_means(read_states(path), n, obs)
    return sm.mean(0), sm.std(0, ddof=1) / np.sqrt(len(sm))


def is_moments(model, ref, n, obs, nblocks=50):
    f = np.concatenate([model.log_ratio(tf.constant(ref[j:j + 20000])).numpy().ravel() for j in range(0, len(ref), 20000)])
    w = np.exp(f - f.max())
    g = obs(ref, n)
    mean = (w[:, None] * g).sum(0) / w.sum()
    ess = w.sum() ** 2 / (w * w).sum() / len(w)
    blocks = [(wb[:, None] * gb).sum(0) / wb.sum() for wb, gb in zip(np.array_split(w, nblocks), np.array_split(g, nblocks))]
    se = np.std(blocks, axis=0, ddof=1) / np.sqrt(nblocks)
    return mean, se, ess


def zrms(z):
    return float(np.sqrt(np.mean(np.square(z))))


def main():
    RES, VAL, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    sys.path.insert(0, str(RES))
    global N3GibbsRatio, N4GibbsRatio
    from n3_model import N3GibbsRatio
    from n4_model import N4GibbsRatio
    out = {}
    for n in (3, 4):
        ref = read_states(RES / f"gibbs_n{n}_T5.bin")
        names = check_names(n)
        for case in ("eq", "neq"):
            mA, sA = val_stats(VAL / f"n{n}_{case}_valA.bin", n, check_obs)
            mB, sB = val_stats(VAL / f"n{n}_{case}_valB.bin", n, check_obs)
            mV, sV = (mA + mB) / 2, np.hypot(sA, sB) / 2
            eA, _ = val_stats(VAL / f"n{n}_{case}_valA.bin", n, extra_obs)
            eB, _ = val_stats(VAL / f"n{n}_{case}_valB.bin", n, extra_obs)
            rec = {"names": names, "val_long_mean": mV.tolist(), "val_long_se": sV.tolist(),
                   "val_long_profile_currents": ((eA + eB) / 2).tolist(),
                   "z_rms_long_A_vs_B": zrms((mA - mB) / np.hypot(sA, sB))}
            # the old short validation trajectories (as used in training), against the long reference
            old = np.concatenate([read_states(RES / f"n{n}_{case}_val{k}.bin") for k in (1, 2)])
            mo = check_obs(old, n).mean(0)
            rec["old_val_mean"] = mo.tolist()
            rec["old_val_minus_long_over_long_se"] = ((mo - mV) / sV).round(2).tolist()
            if case == "eq":  # exact answer: the Gibbs measure itself
                gm, gse, _ = is_moments(type("U", (), {"log_ratio": staticmethod(lambda s: tf.zeros(tf.shape(s)[0], tf.float64))}), ref, n, check_obs)
                rec["gibbs_exact_mean"] = gm.tolist()
                rec["z_rms_long_vs_gibbs"] = zrms((mV - gm) / np.hypot(sV, gse))
                rec["z_long_vs_gibbs"] = dict(zip(names, ((mV - gm) / np.hypot(sV, gse)).round(2).tolist()))
            runs = [f"n{n}_eq"] if case == "eq" else [f"n{n}_neq_s1", f"n{n}_neq_s2"]
            for run in runs:
                cfg = json.loads((RES / run / "model_config.json").read_text())
                if n == 3:
                    model = N3GibbsRatio(np.array(cfg["action_scales"]), cfg["seed"], width=cfg["width"], depth=cfg["depth"])
                else:
                    model = N4GibbsRatio(np.array(cfg["action_scales"]), cfg["seed"], width=cfg["width"],
                                         feature_version=cfg["feature_version"], depth=cfg["depth"])
                _ = model(tf.constant(ref[:2]))
                model.load_weights(str(RES / run / "model.weights.h5"))
                mm, mse, ess = is_moments(model, ref, n, check_obs)
                z = (mm - mV) / np.hypot(sV, mse)
                z_valonly = (mm - mV) / sV
                em, _, _ = is_moments(model, ref, n, extra_obs)
                old = json.loads((RES / run / "result.json").read_text())["final"]
                rec[run] = {"model_mean": mm.tolist(), "is_se": mse.tolist(), "ess_fraction_full_ref": ess,
                            "z": dict(zip(names, z.round(2).tolist())), "z_rms": zrms(z),
                            "z_rms_val_se_only": zrms(z_valonly),
                            "max_rel_dev": float(np.max(np.abs(mm - mV) / np.abs(mV).clip(1e-2))),
                            "model_profile_currents": em.tolist(),
                            "old_check_z_rms": old["check_z_rms"]}
                print(n, case, run, "z_rms(long) = %.2f  [val-se only %.2f]  old %.2f  ESS %.2f"
                      % (zrms(z), zrms(z_valonly), old["check_z_rms"], ess), flush=True)
            print(n, case, "long A vs B z_rms = %.2f" % rec["z_rms_long_A_vs_B"],
                  ("  long vs Gibbs exact z_rms = %.2f" % rec["z_rms_long_vs_gibbs"]) if case == "eq" else "", flush=True)
            out[f"n{n}_{case}"] = rec
    OUT.write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
