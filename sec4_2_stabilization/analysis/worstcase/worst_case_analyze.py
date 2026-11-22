"""batch_07 analysis: worst-case start (I_1 = I_n = 0.01, middle modes share M0 equally), n = 3..6.
Usage: python3 worst_case_analyze.py OUTDIR_OF_BATCH TASKS_TSV RESULTDIR
  OUTDIR_OF_BATCH: the batch's out/ folder with worst/, hires/, ness/ (rk4a_chain_v2 runs)
Per run (128 trajectories, sampled every 0.1):
  - M(t) = sum_j I_j(t): mean, median, 10% / 90% quantiles over trajectories;
  - last t with M(t) >= L for L = 30, 50, M0/2 (per trajectory; if M >= L at the last sample the trajectory is counted
    as not finished): median, 90% quantile, max, number not finished;
  - first t with M(t) <= 2 <M>_NESS (as in the report);
  - at t = 5 and t = 50: fraction of trajectories in which two neighbouring middle modes hold >= 90% of M with
    cos(angle between them) < -0.9.
hires runs (M0 = 100, first time unit, every 1e-3): spread of the endpoint angles around 4pi/3, and M(t) at early times.
<M>_NESS at (2,8), gamma = 0.1 (rk4a_chain runs): n=3 4.686 and n=4 5.237 (batch_06 ness/), n=5 5.727 (batch_06 n5/),
n=6 from this batch (ness/n6_T2_8)."""
import json, sys, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
B, TASKS, O = sys.argv[1:4]
os.makedirs(O, exist_ok=True)
args = {l.split("\t")[0] + "/" + l.split("\t")[1]: l.rstrip("\n").split("\t")[2].split() for l in open(TASKS) if l.strip()}
gv = lambda name, k: args[name][args[name].index(k) + 1]
MN = {3: 4.686, 4: 5.237, 5: 5.727}
p6 = f"{B}/ness/n6_T2_8_dump.f32"
if os.path.exists(p6):
    MN[6] = float(np.fromfile(p6, dtype=np.float32).reshape(-1, 11)[:, :6].astype(float).sum(-1).mean())


def load(folder, name):
    key = f"{folder}/{name}"
    n = int(gv(key, "--n")); D = 2 * n - 1; S = int(gv(key, "--streams")); iv = float(gv(key, "--interval"))
    a = np.fromfile(f"{B}/{folder}/{name}_dump.f32", dtype=np.float32).reshape(S // 16, -1, 16, D).transpose(0, 2, 1, 3).reshape(S, -1, D)
    return n, a.astype(float), (np.arange(a.shape[1]) + 1) * iv


def last_passage(M, t, L):
    above = M >= L
    out = np.full(M.shape[0], np.nan); notfin = above[:, -1]
    for s in range(M.shape[0]):
        if notfin[s]:
            out[s] = np.inf
        elif above[s].any():
            out[s] = t[np.where(above[s])[0][-1]]
        else:
            out[s] = 0.0                       # never >= L after the first sample
    return out


def stat(x):
    """median and 90% quantile (nearest sample; inf = not finished within the run), largest finite value, count not finished"""
    fin = np.isfinite(x); y = np.where(fin, x, 1e300)
    q = lambda p: (lambda v: float(v) if v < 1e299 else float("inf"))(np.quantile(y, p, method="nearest"))
    return {"median": q(0.5), "q90": q(0.9), "max_finite": float(x[fin].max()) if fin.any() else None, "not_finished": int((~fin).sum())}


def neighbour_pair_fraction(I, L_th, n):
    """fraction of trajectories with I_j + I_{j+1} >= 0.9 M and cos(theta_{j+1}) < -0.9 for some middle pair j, j+1"""
    M = I.sum(-1); hit = np.zeros(I.shape[0], bool)
    for j in range(1, n - 2):                  # 0-based middle modes j, j+1 in 1..n-2
        hit |= ((I[:, j] + I[:, j + 1]) / M >= 0.9) & (np.cos(L_th[:, j]) < -0.9)   # long-chain theta_{j+2} -> index j
    return float(hit.mean())


res = {"M_NESS": MN}; curves = {}
for key in args:
    if not key.startswith("worst/"):
        continue
    name = key.split("/", 1)[1]
    if not os.path.exists(f"{B}/worst/{name}_dump.f32"):
        continue
    n, a, t = load("worst", name)
    I = a[:, :, :n]; Lth = a[:, :, n:]; M = I.sum(-1); M0 = float(name.split("_")[1][1:])
    e = {"n": n, "M0": M0, "phase": name.split("_")[2], "T_run": float(t[-1])}
    for L in (30.0, 50.0, M0 / 2):
        e[f"last_t_M_ge_{L:g}"] = stat(last_passage(M, t, L))
    fp = np.array([t[np.argmax(m <= 2 * MN[n])] if (m <= 2 * MN[n]).any() else np.inf for m in M]) if n in MN else None
    if fp is not None:
        e["first_t_M_le_2MNESS"] = stat(fp)
    for tt in (5, 50):
        k = int(round(tt / 0.1)) - 1
        if k < M.shape[1]:
            e[f"frac_neighbour_middle_pair_t{tt}"] = neighbour_pair_fraction(I[:, k], Lth[:, k], n) if n >= 4 else 0.0
    for tt in (1, 10, 50, 100, 200, 500, 1000):
        k = int(round(tt / 0.1)) - 1
        if k < M.shape[1]:
            e[f"M_mean_t{tt}"] = float(M[:, k].mean())
    res[name] = e
    curves[name] = (t, M.mean(0), np.median(M, 0), np.quantile(M, 0.1, 0), np.quantile(M, 0.9, 0))
    print(name, json.dumps({k: v for k, v in e.items() if k.startswith(("last_t_M_ge_30", "first", "frac"))}), flush=True)
# high-resolution start
hi = {}
for key in args:
    if not key.startswith("hires/"):
        continue
    name = key.split("/", 1)[1]
    if not os.path.exists(f"{B}/hires/{name}_dump.f32"):
        continue
    n, a, t = load("hires", name)
    th1 = -a[:, :, n]; thn = a[:, :, 2 * n - 2]
    cs = lambda th: np.sqrt(-2 * np.log(np.abs(np.exp(1j * th).mean(0))))       # circular std over trajectories
    dev1 = np.abs(np.angle(np.exp(1j * (th1 - 4 * np.pi / 3))))
    hi[name] = {"t": t.tolist(), "circ_std_theta1": cs(th1).tolist(), "circ_std_thetan": cs(thn).tolist(),
                "median_abs_dev_theta1_from_4pi3": np.median(dev1, 0).tolist(), "M_mean": a[:, :, :n].sum(-1).mean(0).tolist(),
                "I1_mean": a[:, :, 0].mean(0).tolist()}
    k = np.argmax(np.median(dev1, 0) > 0.5) if (np.median(dev1, 0) > 0.5).any() else -1
    res[key] = {"t_median_dev_theta1_gt_0.5": float(t[k]) if k >= 0 else None,
                            "circ_std_theta1_at_t0.01": float(cs(th1)[9]), "M_mean_t1": float(a[:, -1, :n].sum(-1).mean())}
json.dump(res, open(f"{O}/worst_case.json", "w"), indent=1)
json.dump(hi, open(f"{O}/worst_case_hires.json", "w"))
np.savez_compressed(f"{O}/worst_case_curves.npz", **{f"{k}__{q}": v for k, c in curves.items() for q, v in zip(("t", "mean", "median", "q10", "q90"), c)})

# ---- figures
C = {"worst": "#eb6834", "anti": "#2a78d6", "random": "#1baf7a"}
ns = sorted({v["n"] for k, v in res.items() if isinstance(v, dict) and "n" in v})
if ns:
    fig, ax = plt.subplots(1, len(ns), figsize=(4.2 * len(ns), 3.6), squeeze=False)
    for i, n in enumerate(ns):
        for ph in ("worst", "anti", "random"):
            k = f"n{n}_M200_{ph}"
            if k in curves:
                t, m, med, q1, q9 = curves[k]
                ax[0, i].plot(t, m, color=C[ph], label=ph)
                ax[0, i].fill_between(t, q1, q9, color=C[ph], alpha=0.12, lw=0)
        ax[0, i].axhline(30, color="grey", ls=":", lw=0.8)
        ax[0, i].set_title(f"n = {n}, M0 = 200", fontsize=9); ax[0, i].set_xlabel("t")
        if ax[0, i].get_legend_handles_labels()[0]: ax[0, i].legend(fontsize=7)
    ax[0, 0].set_ylabel("M(t): mean, band 10-90%")
    fig.tight_layout(); fig.savefig(f"{O}/fig_worst_M_t.png", dpi=150)
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    for j, M0 in enumerate((50, 100, 200)):
        for ph in ("worst", "anti", "random"):
            x, y, q = [], [], []
            for n in ns:
                k = f"n{n}_M{M0}_{ph}"
                if k in res:
                    s = res[k]["last_t_M_ge_30"]; x.append(n); y.append(s["median"]); q.append(s["q90"])
            if x:
                ax[j].plot(x, y, "o-", color=C[ph], label=f"{ph} (median)"); ax[j].plot(x, q, "s--", color=C[ph], alpha=0.6, label=f"{ph} (90%)")
        ax[j].set_title(f"last t with M >= 30, M0 = {M0}", fontsize=9); ax[j].set_xlabel("n"); ax[j].set_xticks(ns)
        if ax[j].get_legend_handles_labels()[0]: ax[j].legend(fontsize=6)
    fig.tight_layout(); fig.savefig(f"{O}/fig_worst_last_passage.png", dpi=150)
if hi:
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    for name, h in hi.items():
        ls = "-" if "worst" in name else "--"
        ax[0].plot(h["t"], h["median_abs_dev_theta1_from_4pi3"], ls, label=name)
        ax[1].plot(h["t"], np.array(h["M_mean"]) / h["M_mean"][0], ls, label=name)
    ax[0].set_xscale("log"); ax[0].set_xlabel("t"); ax[0].set_ylabel(r"median $|\theta_1-4\pi/3|$"); ax[0].legend(fontsize=6)
    ax[1].set_xlabel("t"); ax[1].set_ylabel("mean M(t) / M(first sample)"); ax[1].legend(fontsize=6)
    fig.tight_layout(); fig.savefig(f"{O}/fig_worst_hires.png", dpi=150)
print("done")
