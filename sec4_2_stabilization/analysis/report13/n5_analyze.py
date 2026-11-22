"""4.3 extension, A4: n = 5, (2,8), gamma = 0.1, start from I2 = E0 or I3 = E0 (others 0.01), 128 streams.
Per snapshot: the shortest contiguous block of modes that holds >= 90% of M ("cluster"), and whether all bonds inside
are anti-phase (cos th < -0.9).  tau_M: first t with M <= 2 <M>_NESS(n = 5, (2,8)) (from n5_ness).
Leak: finite-difference dM/dt over consecutive snapshots in which the cluster is the anti-phase pair 2-3 or 3-4 (t >= 5,
M > 15), compared with the rule "only bath modes adjacent to the pair drain it": pair 2-3 -> -2 gamma f T1,
pair 3-4 -> -2 gamma f T5, f = 0.391 (gamma = 0.1, equal anti-phase pair, theory_trap.py).
Usage: python3 n5_analyze.py EXTDIR OUTDIR"""
import json, sys, glob, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
X, O = sys.argv[1:3]
n, D = 5, 9


def load(p, iv):
    a = np.fromfile(p, dtype=np.float32).reshape(8, -1, 16, D).transpose(0, 2, 1, 3).reshape(128, -1, D).astype(float)
    return a, (np.arange(a.shape[1]) + 1) * iv


ness, _ = load(f"{X}/n5/n5_ness_dump.f32", 0.1)
Mn = float(ness[:, :, :n].sum(-1).mean())
res = {"n5_ness_mean_M": Mn, "n5_ness_mean_I": ness[:, :, :n].mean((0, 1)).tolist()}
blocks = [(i, j) for L in range(1, n + 1) for i in range(n - L + 1) for j in [i + L - 1]]   # by length


def cluster(I, th):
    """returns block index (into blocks) of the shortest block with >= 90% of M, and anti-phase flag"""
    M = I.sum(-1); out = np.full(M.shape, -1); anti = np.zeros(M.shape, bool)
    for b, (i, j) in enumerate(blocks):
        frac = I[..., i:j + 1].sum(-1) / M
        new = (out < 0) & (frac >= 0.9)
        if j > i:
            ap = np.all(np.cos(th[..., i:j]) < -0.9, axis=-1)       # bonds i+1..j  -> th index i..j-1 (th_{k+1} at k)
        else:
            ap = np.ones(M.shape, bool)
        out[new] = b; anti[new] = ap[new]
    return out, anti


fig, ax = plt.subplots(1, 3, figsize=(16, 4.3))
for p in sorted(glob.glob(f"{X}/n5/n5_I*_E*_dump.f32")):
    name = os.path.basename(p).replace("_dump.f32", "")
    a, t = load(p, 0.05)
    I, th = a[:, :, :n], a[:, :, n:]
    M = I.sum(-1)
    tau = np.array([t[np.argmax(r)] if r.any() else np.inf for r in (M <= 2 * Mn)])
    cl, anti = cluster(I, th)
    k5 = int(5 / 0.05) - 1
    hot = (M > 20) & (t >= 5)[None, :]
    lab = lambda b: "-".join(str(k + 1) for k in range(blocks[b][0], blocks[b][1] + 1))
    types5 = {}
    for b in np.unique(cl[:, k5]):
        sel = cl[:, k5] == b
        types5[lab(b)] = {
            "fraction": float(sel.mean()), "anti_phase_fraction": float(anti[sel, k5].mean()),
            "tau_M_median": float(np.median(tau[sel]))}
    hotd = {}
    if hot.any():
        for b in np.unique(cl[hot]):
            m = hot & (cl == b)
            hotd[lab(b)] = {"time_fraction": float(m.sum() / hot.sum()), "anti_phase": float(anti[m].mean())}
    fd = (M[:, 1:] - M[:, :-1]) / 0.05
    leak = {}
    for lbl, bidx, Tb in [("2-3", blocks.index((1, 2)), 2.0), ("3-4", blocks.index((2, 3)), 8.0)]:
        sel = (cl[:, :-1] == bidx) & (cl[:, 1:] == bidx) & anti[:, :-1] & anti[:, 1:] & (t[None, :-1] >= 5) & (M[:, :-1] > 15)
        if sel.sum() > 100:
            gm = [fd[g][sel[g]].mean() for g in np.array_split(np.arange(128), 16) if sel[g].sum() > 50]
            leak[lbl] = {"dMdt_fd": float(fd[sel].mean()), "se": float(np.std(gm, ddof=1) / np.sqrt(len(gm))), "pred": -2 * 0.1 * 0.390625 * Tb, "snapshots": int(sel.sum())}
    e = {"leak_by_pair": leak, "tau_M_median": float(np.median(tau)), "tau_M_q75": float(np.quantile(tau, .75)), "tau_M_q90": float(np.quantile(tau, .9)),
         "not_relaxed": int(np.isinf(tau).sum()), "T_run": float(t[-1]), "cluster_at_t5": types5, "cluster_while_M>20_after_t5": hotd}
    res[name] = e
    print(name, json.dumps(e), flush=True)
    k = 0 if "_I2_" in name else 1
    ax[k].plot(t, M.mean(0), label=name)
    s = np.sort(tau); fin = np.isfinite(s)
    ax[2].step(np.append(0, s[fin]), 1 - np.arange(fin.sum() + 1) / len(s), where="post", label=name)
for a_, ttl in [(ax[0], "n=5, start I2=E0: <M>(t)"), (ax[1], "n=5, start I3=E0: <M>(t)"), (ax[2], "survival P(tau_M > t)")]:
    a_.set_xscale("log"); a_.set_title(ttl, fontsize=9); a_.legend(fontsize=6); a_.set_xlabel("t")
ax[0].set_yscale("log"); ax[1].set_yscale("log")
for a_ in ax[:2]:
    a_.axhline(Mn, color="grey", ls=":")
fig.tight_layout(); fig.savefig(f"{O}/fig_n5.png", dpi=150)
json.dump(res, open(f"{O}/n5.json", "w"), indent=1)
