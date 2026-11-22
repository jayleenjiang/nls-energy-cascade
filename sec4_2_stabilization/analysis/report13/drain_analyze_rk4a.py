"""[rk4a copy 2026-10-27: dt-check block removed]
4.3.4a: drain of a large middle-mode excursion.  Runs: drain/n{3,4}_E{10,30,100}_dump.f32 (interval 0.002).
Usage: python3 drain_analyze.py HIDIR OUTDIR
Per stream: drain time tau_d = first t with I2 <= E0/2; angles (advisor-note convention) during t < tau_d;
mean of the endpoint growth rate d log I1/dt = 4 I2 sin th1 over t < tau_d (perfect locking at th1 = 2pi/3 gives
4 sin(2pi/3) I2 = 3.46 I2); time for M to fall below 2 x its NESS mean."""
import json, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
H, O = sys.argv[1:3]
DT_REC = 0.002


def load(path, n, nb=16):
    D = 2 * n - 1
    a = np.fromfile(path, dtype=np.float32).reshape(nb, -1, 16, D)
    return a.transpose(0, 2, 1, 3).reshape(nb * 16, -1, D).astype(np.float64)   # (stream, snap, D)


def ness_mean_M(n):
    a = load(f"{H}/ness/n{n}_T2_8_dump.f32", n)
    return float(a[:, :, :n].sum(-1).mean())


out = {}
fig, axes = plt.subplots(2, 3, figsize=(14, 7))
for ni, n in enumerate([3, 4]):
    Mness = ness_mean_M(n)
    out[f"n{n}_ness_mean_M"] = Mness
    for E0 in [10, 30, 100]:
        a = load(f"{H}/drain/n{n}_E{E0}_dump.f32", n)
        S, K, D = a.shape
        t = (np.arange(K) + 1) * DT_REC
        I = a[:, :, :n]; M = I.sum(-1)
        th2 = a[:, :, n]                         # long-chain theta_2 ; note angle th1 = -theta_2
        th1n = -th2
        tau = np.array([t[np.argmax(I[s, :, 1] <= E0 / 2)] if np.any(I[s, :, 1] <= E0 / 2) else np.nan for s in range(S)])
        mask = t[None, :] < tau[:, None]
        rate1 = 4 * I[:, :, 1] * np.sin(th1n)   # d log I1 / dt
        sel1 = th1n[mask]; selr = rate1[mask]; selI2 = I[:, :, 1][mask]
        tauM = np.array([t[np.argmax(M[s] <= 2 * Mness)] if np.any(M[s] <= 2 * Mness) else np.nan for s in range(S)])
        e = {"E0": E0, "n": n, "streams": S,
             "tau_drain_median": float(np.nanmedian(tau)), "tau_drain_q10_q90": [float(np.nanquantile(tau, .1)), float(np.nanquantile(tau, .9))],
             "undrained_streams": int(np.isnan(tau).sum()),
             "mean_sin_th1_during_drain": float(np.mean(np.sin(sel1))),
             "mean_rate_over_I2_during_drain": float(np.mean(selr / selI2)),
             "tau_M_median": float(np.nanmedian(tauM)), "tau_M_q90": float(np.nanquantile(tauM, .9)), "M_not_relaxed": int(np.isnan(tauM).sum())}
        if n == 3:
            th3n = a[:, :, n + 1]
            e["mean_sin_th3_during_drain"] = float(np.mean(np.sin(th3n[mask])))
        out[f"n{n}_E{E0}"] = e
        print(e, flush=True)
        ax = axes[ni, 0]
        ax.plot(t, I[:, :, 1].mean(0) / E0, label=f"E0={E0}")
        ax.set_xscale("log"); ax.set_xlabel("t"); ax.set_ylabel(r"$\langle I_2\rangle/E_0$"); ax.set_title(f"n={n}: middle mode drains", fontsize=9)
        ax = axes[ni, 1]
        h, edges = np.histogram(sel1, bins=72, range=(-np.pi, np.pi), density=True)
        ax.plot(0.5 * (edges[1:] + edges[:-1]), h, label=f"E0={E0}")
        ax.axvline(2 * np.pi / 3, color="grey", ls=":"); ax.axvline(-2 * np.pi / 3, color="grey", ls=":")
        ax.set_xlabel(r"$\theta_1$ (notes convention) while $I_2>E_0/2$"); ax.set_ylabel("density"); ax.set_title(f"n={n}: phase during the drain", fontsize=9)
        ax = axes[ni, 2]
        ax.plot(t, M.mean(0), label=f"E0={E0}")
        ax.axhline(Mness, color="grey", ls=":"); ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("t"); ax.set_ylabel(r"$\langle M\rangle$"); ax.set_title(f"n={n}: total action returns to the NESS (dotted)", fontsize=9)
    for ax in axes[ni]:
        ax.legend(fontsize=7)
# (time-step check removed: rk4a n3 drains have one dt per E0; the n4 dt check is in trap_scan)
fig.tight_layout(); fig.savefig(f"{O}/fig_drain.png", dpi=150)
json.dump(out, open(f"{O}/drain.json", "w"), indent=1)
