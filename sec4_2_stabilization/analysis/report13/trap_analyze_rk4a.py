"""[rk4a copy 2026-10-27: n = 4 only (n = 3 return times come from the drain runs); MN from rk4a NESS]
4.3.4c: long relaxation from a large excursion of I2 (others 0.01), interval 0.05, 64 streams.
tau_M = first time M <= 2 <M>_NESS.  For n = 4 also: state of the trapped energy (I2, I3, the angle between them)
and the endpoint growth rates 4 I_mid sin(angle) while M > 0.5 E0.
Usage: python3 trap_analyze.py HIDIR OUTDIR"""
import json, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
H, O = sys.argv[1:3]
MN = {n_: float(np.fromfile(f"{H}/ness/n{n_}_T2_8_dump.f32", dtype=np.float32).reshape(-1, 2 * n_ - 1)[:, :n_].astype(float).sum(-1).mean()) for n_ in (3, 4)}   # rk4a NESS means


def load(p, n, nb=4):
    D = 2 * n - 1
    a = np.fromfile(p, dtype=np.float32).reshape(nb, -1, 16, D)
    return a.transpose(0, 2, 1, 3).reshape(nb * 16, -1, D).astype(float)


res = {}
fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
for n, E0s, T in [(4, [20, 30, 50, 100, 200], 200)]:
    for E0 in E0s:
        a = load(f"{H}/trap/n{n}_E{E0}_dump.f32", n)
        t = (np.arange(a.shape[1]) + 1) * 0.05
        I = a[:, :, :n]; M = I.sum(-1)
        rel = M <= 2 * MN[n]
        tau = np.array([t[np.argmax(r)] if r.any() else np.inf for r in rel])
        e = {"n": n, "E0": E0, "streams": len(tau), "tau_M_median": float(np.median(tau)), "tau_M_q90": float(np.quantile(tau, .9)),
             "not_relaxed_by_T": int(np.isinf(tau).sum()), "T": T, "M_mean_at": {str(tt): float(M[:, int(tt / 0.05) - 1].mean()) for tt in [1, 10, 50, T]}}
        if n == 4:
            hot = M > 0.5 * E0                     # trapped stage
            if hot.any():
                th3 = a[:, :, 4 + 1]                # long-chain theta_3 (between I2 and I3)
                left2 = -a[:, :, 4]                 # angle 2(phi1 - phi2)
                right3 = a[:, :, 4 + 2]             # theta_4 = 2(phi4 - phi3)
                e["trapped_fraction_of_M_in_I2_I3"] = float(((I[:, :, 1] + I[:, :, 2]) / M)[hot].mean())
                e["trapped_mean_cos_theta3"] = float(np.cos(th3[hot]).mean())
                e["trapped_endpoint1_rate"] = float((4 * I[:, :, 1] * np.sin(left2))[hot].mean())      # d log I1/dt
                e["trapped_endpoint4_rate"] = float((-4 * I[:, :, 2] * np.sin(right3))[hot].mean())    # d log I4/dt
                e["trapped_I2_minus_I3_rms"] = float(np.sqrt(((I[:, :, 1] - I[:, :, 2]) ** 2)[hot].mean()))
                # slow leak: fit log M on [5, 50] for streams still hot over that window
                w = (t >= 5) & (t <= 50)
                keep = hot[:, w].all(1)
                if keep.sum() >= 4:
                    sl = np.array([np.polyfit(t[w], np.log(M[s, w]), 1)[0] for s in np.where(keep)[0]])
                    e["leak_rate_dlogM_dt_hot_streams"] = [float(sl.mean()), float(sl.std(ddof=1) / np.sqrt(len(sl))), int(keep.sum())]
        res[f"n{n}_E{E0}"] = e
        print(e, flush=True)
        ax = axes[0 if n == 4 else 2]
        ax.plot(t, M.mean(0), label=f"E0={E0}")
        if n == 4:
            s = np.sort(tau)
            axes[1].step(np.append(0, s[np.isfinite(s)]), 1 - np.arange(np.isfinite(s).sum() + 1) / len(s), where="post", label=f"E0={E0}")
for ax, title in [(axes[0], "n=4: mean total action M(t)"), (axes[2], "n=3: mean total action M(t)")]:
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("t"); ax.set_ylabel(r"$\langle M\rangle$"); ax.set_title(title, fontsize=9); ax.legend(fontsize=7)
axes[0].axhline(MN[4], color="grey", ls=":"); axes[2].axhline(MN[3], color="grey", ls=":")
axes[1].set_xscale("log"); axes[1].set_xlabel("t"); axes[1].set_ylabel(r"fraction with $M>2\langle M\rangle_{NESS}$"); axes[1].set_title("n=4: survival of the excursion", fontsize=9); axes[1].legend(fontsize=7)
fig.tight_layout(); fig.savefig(f"{O}/fig_trap.png", dpi=150)
json.dump(res, open(f"{O}/trap.json", "w"), indent=1)
