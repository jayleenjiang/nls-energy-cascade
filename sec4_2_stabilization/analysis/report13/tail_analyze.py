"""4.3.4b: tail of the total action M in the NESS (T1,Tn) = (2,8) vs Gibbs laws exp(-H/2T) at T = 2, 5, 8.
Gibbs action density (angles integrated exactly):  exp(-(M^2 - sum I^2 / 2) / (2T)) * prod_j I0(I_{j-1} I_j / T).
n = 3: 3D quadrature (cross-checked by Metropolis at T = 5); n = 4: Metropolis on log-actions (4096 chains).
Usage: python3 tail_analyze.py HIDIR OUTDIR"""
import json, sys
import numpy as np
from scipy.special import i0e
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
H, O = sys.argv[1:3]
mgrid = np.linspace(0, 25, 501)


def load(path, n, nb=16):
    D = 2 * n - 1
    return np.fromfile(path, dtype=np.float32).reshape(-1, D).astype(np.float64)


def logw(I, T):
    M = I.sum(-1)
    lw = -(M ** 2 - 0.5 * (I ** 2).sum(-1)) / (2 * T)
    for j in range(1, I.shape[-1]):
        x = I[..., j - 1] * I[..., j] / T
        lw += np.log(i0e(x)) + x
    return lw


def gibbs_ccdf_n3(T, h=0.025, L=22.0):
    J = np.arange(h / 2, L, h)
    A, C = np.meshgrid(J, J, indexing="ij")
    hist = np.zeros(len(mgrid)); tot = 0.0
    # log I0(x) <= x gives logw <= -(sum I^2/2 + 2 I1 I3)/(2T) <= 0, so exp(logw) needs no rescaling
    for i2 in J:
        I = np.stack([A, np.full_like(A, i2), C], -1)
        w = np.exp(logw(I, T))
        hist += np.histogram(A + i2 + C, bins=np.append(mgrid, 1e9), weights=w)[0]
        tot += w.sum()
    ccdf = 1 - np.cumsum(hist) / tot + hist / tot          # P(M >= mgrid[k]) approximately
    return np.clip(ccdf, 0, 1)


def gibbs_samples(T, n=4, chains=4096, burn=4000, steps=24000, thin=8, seed=3):
    rng = np.random.default_rng(seed)
    x = np.log(np.full((chains, n), np.sqrt(T) * 0.5))      # log actions
    step = 0.35
    lp = logw(np.exp(x), T) + x.sum(-1)                      # Jacobian prod I
    out = []
    for k in range(burn + steps):
        y = x + step * rng.standard_normal(x.shape)
        lq = logw(np.exp(y), T) + y.sum(-1)
        acc = np.log(rng.random(chains)) < lq - lp
        x[acc] = y[acc]; lp[acc] = lq[acc]
        if k >= burn and k % thin == 0:
            out.append(np.exp(x).sum(-1))
    return np.concatenate(out)


def emp_ccdf(M):
    M = np.sort(M)
    return np.array([1 - np.searchsorted(M, m, side="left") / len(M) for m in mgrid])


res = {}
fig, axes = plt.subplots(1, 2, figsize=(12, 4.3))
for ax, n in zip(axes, [3, 4]):
    a = load(f"{H}/ness/n{n}_T2_8_dump.f32", n)
    Mn = a[:, :n].sum(-1)
    c = emp_ccdf(Mn)
    ax.semilogy(mgrid, c, "k", lw=2, label="NESS (2,8)")
    res[f"n{n}_NESS"] = {"mean_M": float(Mn.mean()), "q99": float(np.quantile(Mn, .99)), "q999": float(np.quantile(Mn, .999)), "max": float(Mn.max()), "samples": int(len(Mn))}
    for T in [2.0, 5.0, 8.0]:
        if n == 3:
            g = gibbs_ccdf_n3(T)
            q = lambda p: float(np.interp(p, 1 - g, mgrid))   # quantile: CDF = 1 - ccdf = p
            res[f"n3_Gibbs_T{T:g}"] = {"q99": q(.99), "q999": q(.999)}
            if T == 5.0:   # cross-check quadrature vs Metropolis
                s3 = gibbs_samples(T, 3)
                res["n3_Gibbs_T5_metropolis"] = {"q99": float(np.quantile(s3, .99)), "q999": float(np.quantile(s3, .999)),
                                                 "max_abs_ccdf_diff": float(np.max(np.abs(emp_ccdf(s3) - g)))}
        else:
            s = gibbs_samples(T, 4)
            g = emp_ccdf(s)
            res[f"n4_Gibbs_T{T:g}"] = {"mean_M": float(s.mean()), "q99": float(np.quantile(s, .99)), "q999": float(np.quantile(s, .999)), "samples": int(len(s))}
        ax.semilogy(mgrid, g, ls="--", label=f"Gibbs T={T:g}")
    ax.set_ylim(1e-7, 1.2); ax.set_xlabel("M"); ax.set_ylabel(r"$P(M\geq m)$"); ax.set_title(f"n={n}: tail of the total action", fontsize=9); ax.legend(fontsize=7)
    print(n, res, flush=True)
fig.tight_layout(); fig.savefig(f"{O}/fig_M_tail.png", dpi=150)
json.dump(res, open(f"{O}/tail.json", "w"), indent=1)
