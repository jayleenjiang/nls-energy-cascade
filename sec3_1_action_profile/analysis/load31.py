"""Load Sec. 3.1 v2 summaries; merge replicates (stream-weighted means, SE from per-replicate SE)."""
import glob
import os
import re

import numpy as np

ROOT = os.environ.get("S31", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))   # folder with atlas/, closed/, dense/, dtcheck/


def read_summary(path):
    d = {}
    for line in open(path):
        k, *v = line.split()
        d[k] = v
    out = {k: d[k][0] for k in ("bc", "n", "T1", "Tn", "closed", "dt", "streams", "midpoint_failures", "nonfinite")}
    out["n"] = int(out["n"]); out["streams"] = int(out["streams"])
    for k in ("mean_I", "se_I", "mean_current"):
        out[k] = np.array([float(x) for x in d[k]])
    out["checkpoints"] = [np.array([float(x) for x in d[k]]) for k in sorted(d) if k.startswith("checkpoint_") and k.endswith("_mean_I")]
    out["ok"] = out["midpoint_failures"] == "0" and out["nonfinite"] == "0" and np.all(np.isfinite(out["mean_I"]))
    return out


def merge(reps):
    reps = [r for r in reps if r["ok"]]
    if not reps:
        return None
    w = np.array([r["streams"] for r in reps], float); w /= w.sum()
    m = sum(wi * r["mean_I"] for wi, r in zip(w, reps))
    se = np.sqrt(sum((wi * r["se_I"]) ** 2 for wi, r in zip(w, reps)))
    J = sum(wi * r["mean_current"] for wi, r in zip(w, reps))
    return {"mean_I": m, "se_I": se, "mean_current": J, "streams": int(sum(r["streams"] for r in reps)),
            "n_reps": len(reps), "reps": reps}


def load(group):
    files = sorted(glob.glob(f"{ROOT}/{group}/*_summary.txt"))
    out = {}
    for f in files:
        lab = os.path.basename(f)[:-len("_summary.txt")]
        key = re.sub(r"_rep\d+$", "", lab)
        out.setdefault(key, []).append(read_summary(f))
    return out
