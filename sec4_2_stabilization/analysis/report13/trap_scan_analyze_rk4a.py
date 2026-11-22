"""4.3 extension, A1-A3: n = 4 middle-pair trap versus E0, gamma and bath temperatures, and tests of the linear theory
(theory_trap.py).  Runs: ext/e0scan, ext/gscan, ext/tscan (128 streams, start I2 = E0, others 0.01, random phases).
Definitions (per stream):
  pair state  : (I2 + I3)/M > 0.9 and cos th3 < -0.9            (th3 = 2(phi3 - phi2))
  formed      : pair state at t = 5
  tau_M       : first t with M <= 2 <M>_ref  (<M>_ref: NESS mean at (2,8) and its mirror (8,2), Gibbs mean for T1 = Tn)
Leak (primary): finite difference (M(t+dt) - M(t))/dt over consecutive pair-state snapshots with t >= 5.
Leak (cross-check): exact identity dM/dt = sum_ends 2 gamma (2 T_e - I_e dH/dI_e), dH/dI_1 = 2M - I1 + 2 I2 cos th2,
      dH/dI_4 = 2M - I4 + 2 I3 cos th4, averaged over pair-state snapshots with t >= 5.
      Prediction (frozen pair, gamma -> 0): dM/dt = -2 gamma (T1 f(I2, I3, th3) + T4 f(I3, I2, th3)),
      f = kappa^2 (2M - Delta)/(M (Delta^2 - kappa^2)), Delta = I_a - 2 I_b cos th3, kappa = 2 I_a; f = 0.4 for I2 = I3.
      Finite gamma, equal pair: f(gamma) from the stationary Lyapunov equation (theory_trap.f_exact): f(0.1) = 0.391.
      Endpoint action: <I_end> = 0.89 T_end / I (equal pair, gamma = 0.1).
Usage: python3 trap_scan_analyze_rk4a.py EXTDIR THEORYJSON OUTDIR TASKS_TSV NESS_N4_DUMP
[rk4a copy 2026-10-27: interval from tasks.tsv; M_ref for (2,8),(8,2) from the rk4a n = 4 NESS run]"""
import json, sys, glob, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
X, TJ, O, TASKS, NESS4 = sys.argv[1:6]
IV = {l.split("\t")[0] + "/" + l.split("\t")[1]: float(l.split("\t")[2].split("--interval ")[1].split()[0]) for l in open(TASKS) if l.strip()}


def gibbs_mean_M(T, n=4, chains=4096, burn=3000, steps=12000, thin=6, seed=3):
    from scipy.special import i0e
    rng = np.random.default_rng(seed)
    def logw(I):
        M = I.sum(-1); lw = -(M ** 2 - 0.5 * (I ** 2).sum(-1)) / (2 * T)
        for j in range(1, n):
            x = I[..., j - 1] * I[..., j] / T; lw += np.log(i0e(x)) + x
        return lw
    x = np.log(np.full((chains, n), np.sqrt(T) * 0.5)); lp = logw(np.exp(x)) + x.sum(-1); acc_M = []
    for k in range(burn + steps):
        y = x + 0.35 * rng.standard_normal(x.shape); lq = logw(np.exp(y)) + y.sum(-1)
        a = np.log(rng.random(chains)) < lq - lp; x[a] = y[a]; lp[a] = lq[a]
        if k >= burn and k % thin == 0:
            acc_M.append(np.exp(x).sum(-1).mean())
    return float(np.mean(acc_M))


_m4 = float(np.fromfile(NESS4, dtype=np.float32).reshape(-1, 7)[:, :4].astype(float).sum(-1).mean())
MREF = {(2.0, 8.0): _m4, (8.0, 2.0): _m4, (5.0, 5.0): 5.1246}
MREF[(1.0, 1.0)] = gibbs_mean_M(1.0)


def f_frozen(Ia, Ib, c):
    M = Ia + Ib; D = Ia - 2 * Ib * c; k = 2 * Ia
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(np.abs(D) > k, k ** 2 * (2 * M - D) / (M * (D ** 2 - k ** 2)), np.nan)


def summ(path):
    d = {}
    for line in open(path):
        k, *v = line.split()
        d[k] = v[0] if v else ""
    return d


def analyze(name, path):
    s = summ(path.replace("_dump.f32", "_summary.txt"))
    g, T1, T4, iv = float(s["gamma"]), float(s["T1"]), float(s["Tn"]), IV[os.path.basename(os.path.dirname(path)) + "/" + os.path.basename(path).replace("_dump.f32", "")]
    a = np.fromfile(path, dtype=np.float32).reshape(8, -1, 16, 7).transpose(0, 2, 1, 3).reshape(128, -1, 7).astype(float)
    t = (np.arange(a.shape[1]) + 1) * iv
    I = a[:, :, :4]; th2, th3, th4 = a[:, :, 4], a[:, :, 5], a[:, :, 6]
    M = I.sum(-1)
    pair = ((I[:, :, 1] + I[:, :, 2]) / M > 0.9) & (np.cos(th3) < -0.9)
    k5 = int(round(5 / iv)) - 1
    formed = pair[:, k5]
    ref = MREF[(T1, T4)]
    rel = M <= 2 * ref
    tau = np.array([t[np.argmax(r)] if r.any() else np.inf for r in rel])
    late = (t >= 5)[None, :] & pair
    dM = 2 * g * (2 * T1 - I[:, :, 0] * (2 * M - I[:, :, 0] + 2 * I[:, :, 1] * np.cos(th2))) \
        + 2 * g * (2 * T4 - I[:, :, 3] * (2 * M - I[:, :, 3] + 2 * I[:, :, 2] * np.cos(th4)))
    f1 = f_frozen(I[:, :, 1], I[:, :, 2], np.cos(th3)); f4 = f_frozen(I[:, :, 2], I[:, :, 1], np.cos(th3))
    pred = -2 * g * (T1 * f1 + T4 * f4)
    Ipair = 0.5 * (I[:, :, 1] + I[:, :, 2])
    e = {"gamma": g, "T1": T1, "Tn": T4, "E0": float(s.get("E0", 0)) if "E0" in s else None, "streams": 128,
         "M_ref": ref, "formed_fraction_t5": float(formed.mean()),
         "tau_M_median_all": float(np.median(tau)), "not_relaxed": int(np.isinf(tau).sum()), "T_run": float(t[-1]),
         "tau_M_median_formed": float(np.median(tau[formed])) if formed.any() else None,
         "tau_M_quartiles_formed": [float(np.quantile(tau[formed], q)) for q in (.25, .75)] if formed.sum() > 4 else None,
         "tau_M_median_not_formed": float(np.median(tau[~formed])) if (~formed).any() else None,
         "pair_snapshots_t>=5": int(late.sum())}
    fd = (M[:, 1:] - M[:, :-1]) / iv
    selfd = pair[:, :-1] & pair[:, 1:] & (t[None, :-1] >= 5)
    # M when the stream leaves the pair state for the last time before tau_M (formed streams)
    Mend = []
    for k in np.where(formed)[0]:
        kk = np.where(pair[k] & (t < tau[k]))[0]
        if len(kk):
            Mend.append(M[k, kk[-1]])
    e["M_at_t5_formed_mean"] = float(M[formed, k5].mean()) if formed.any() else None
    e["M_at_pair_end_median"] = float(np.median(Mend)) if Mend else None
    e["M_at_pair_end_quartiles"] = [float(np.quantile(Mend, .25)), float(np.quantile(Mend, .75))] if len(Mend) > 4 else None
    if selfd.sum() > 100:
        nb = 16; grp_ = np.array_split(np.arange(128), nb)
        gm = [fd[gg][selfd[gg]].mean() for gg in grp_ if selfd[gg].sum() > 50]
        from scipy.linalg import solve_continuous_lyapunov
        A = np.array([[-g * 6, -1.0], [5.0, -g * 2]]); S = solve_continuous_lyapunov(A, -2 * g * np.eye(2))
        fg = 6 * S[0, 0] + 2 * S[1, 1] - 2
        e.update({"dMdt_fd_pair": float(fd[selfd].mean()), "dMdt_fd_pair_se": float(np.std(gm, ddof=1) / np.sqrt(len(gm))),
                  "f_equal_pair_gamma": float(fg), "dMdt_pred_equal_pair_gamma": float(-2 * g * fg * (T1 + T4))})
    if late.sum() > 100:
        fin = np.isfinite(pred[late])
        e.update({"dMdt_measured_pair": float(dM[late].mean()), "dMdt_measured_se": float(dM[late].std() / np.sqrt(late.sum() / 50)),
                  "dMdt_pred_frozen_median": float(np.median(pred[late][fin])), "dMdt_pred_frozen_mean": float(pred[late][fin].mean()),
                  "frac_pair_snapshots_frozen_unstable": float((~fin).mean()),
                  "dMdt_pred_equal_pair": float(-2 * g * 0.4 * (T1 + T4)),
                  "endpoint_I1_times_I_over_T1": float((I[:, :, 0] * Ipair / T1)[late].mean()),
                  "endpoint_I4_times_I_over_T4": float((I[:, :, 3] * Ipair / T4)[late].mean()),
                  "pair_cos_th3_mean": float(np.cos(th3)[late].mean()), "pair_log_ratio_rms": float(np.sqrt((np.log(I[:, :, 2] / I[:, :, 1]) ** 2)[late].mean())),
                  "pair_state_points": {"logr": np.log(I[:, :, 2] / I[:, :, 1])[late][::max(1, late.sum() // 20000)].tolist(),
                                        "th3": np.mod(th3[late], 2 * np.pi)[::max(1, late.sum() // 20000)].tolist()}})
    return e, t, M, tau, formed


res = {}
fig, ax = plt.subplots(2, 3, figsize=(16, 9))
for grp in ["e0scan", "gscan", "tscan"]:
    for p in sorted(glob.glob(f"{X}/{grp}/n4_E*_dump.f32")):
        name = os.path.basename(p).replace("_dump.f32", "")
        if "dt1.25e-5" in name:
            continue
        E0 = float(name.split("_")[1][1:])
        e, t, M, tau, formed = analyze(name, p); e["E0"] = E0
        pts = e.pop("pair_state_points", None)
        res[f"{grp}/{name}"] = e
        if pts is not None:
            res[f"{grp}/{name}"]["_pts"] = pts
        print(grp, name, {k: v for k, v in e.items() if not k.startswith("_")}, flush=True)
        if grp == "e0scan":
            ax[0, 0].plot(t, M.mean(0), label=f"E0={E0:g}")
            s = np.sort(tau); fin = np.isfinite(s)
            ax[0, 1].step(np.append(0, s[fin]), 1 - np.arange(fin.sum() + 1) / len(s), where="post", label=f"E0={E0:g}")
# dt check
a = np.fromfile(f"{X}/e0scan/n4_E400_dt1.25e-5_dump.f32", dtype=np.float32).reshape(8, -1, 16, 7).transpose(0, 2, 1, 3).reshape(128, -1, 7)
b = np.fromfile(f"{X}/e0scan/n4_E400_dump.f32", dtype=np.float32).reshape(8, -1, 16, 7).transpose(0, 2, 1, 3).reshape(128, -1, 7)
chk = {}
for tt in (10, 50, 100):
    k = int(tt / 0.05) - 1
    ma, mb = a[:, k, :4].sum(-1), b[:, k, :4].sum(-1)
    chk[f"t{tt}"] = {"mean_M_dt2.5e-5": float(mb.mean()), "mean_M_dt1.25e-5": float(ma.mean()), "se": float(mb.std(ddof=1) / np.sqrt(128))}
res["dt_check_E400"] = chk
print(chk)
# lifetime vs E0 and gamma
E = sorted([v["E0"] for k, v in res.items() if k.startswith("e0scan")])
med = [res[f"e0scan/n4_E{e:g}"]["tau_M_median_formed"] for e in E]
q = [res[f"e0scan/n4_E{e:g}"]["tau_M_quartiles_formed"] for e in E]
ax[0, 2].errorbar(E, med, yerr=np.array([[m - qq[0], qq[1] - m] for m, qq in zip(med, q)]).T, fmt="o-", label="gamma=0.1, (2,8)")
for gg, mk in [(0.05, "s"), (0.02, "^")]:
    pts = [(v["E0"], v["tau_M_median_formed"]) for k, v in res.items() if k.startswith("gscan") and abs(v["gamma"] - gg) < 1e-9 and v["tau_M_median_formed"]]
    if pts:
        ax[0, 2].plot(*zip(*pts), mk, label=f"gamma={gg}")
ax[0, 2].set_xscale("log"); ax[0, 2].set_yscale("log"); ax[0, 2].set_xlabel("E0"); ax[0, 2].set_ylabel("median tau_M (streams in the pair at t=5)"); ax[0, 2].legend(fontsize=7)
for a_, ttl in [(ax[0, 0], "n=4, (2,8), gamma=0.1: <M>(t)"), (ax[0, 1], "survival P(tau_M > t)")]:
    a_.set_xscale("log"); a_.set_title(ttl, fontsize=9); a_.legend(fontsize=6); a_.set_xlabel("t")
ax[0, 0].set_yscale("log")
# leak: measured vs predicted
keys = [k for k, v in res.items() if not k.startswith("dt") and "dMdt_fd_pair" in v and v["E0"] >= 100]
mx = [res[k]["dMdt_fd_pair"] for k in keys]; se = [res[k]["dMdt_fd_pair_se"] for k in keys]; pe = [res[k]["dMdt_pred_equal_pair_gamma"] for k in keys]
ax[1, 0].errorbar(pe, mx, yerr=se, fmt="o", label="E0 >= 100: e0scan, gscan, tscan")
for k, x_, y_ in zip(keys, pe, mx):
    ax[1, 0].annotate(k.split("/")[1].replace("n4_", ""), (x_, y_), fontsize=5)
lim = [min(mx + pe) * 1.3, 0]; ax[1, 0].plot(lim, lim, "k:"); ax[1, 0].set_xlabel(r"predicted $-2\gamma f(\gamma)(T_1+T_n)$ (equal anti-phase pair)"); ax[1, 0].set_ylabel("measured dM/dt in the pair state")
ax[1, 0].legend(fontsize=7); ax[1, 0].set_title("leak of the trapped pair", fontsize=9)
# pair states vs damped basin
th = json.load(open(TJ)); td = json.load(open(TJ.replace("theory_trap.json", "theory_trap_damped.json")))
rr, tt_ = np.array(th["grid_r"]), np.array(th["grid_th"]); LG = np.array(td["LmaxG"]); gl = td["gammas"]
for a_, key, gg in [(ax[1, 1], "e0scan/n4_E200", 0.1), (ax[1, 2], "gscan/n4_E200_g0.02", 0.02)]:
    if key in res and "_pts" in res[key]:
        P = res[key]["_pts"]
        a_.hist2d(P["th3"], P["logr"], bins=[90, 60], range=[[0, 2 * np.pi], [np.log(rr[0]), np.log(rr[-1])]], cmap="Greys", cmin=1)
    for i, g_ in enumerate(gl):
        if g_ > 0:
            a_.contour(tt_, np.log(rr), (LG[i] < 0).astype(float), levels=[0.5], colors={0.02: "c", 0.05: "g", 0.1: "r"}[g_], linewidths=1)
    a_.set_xlim(np.pi - 1.2, np.pi + 1.2); a_.set_ylim(-1.2, 1.2)
    a_.set_xlabel(r"$\theta_3$"); a_.set_ylabel(r"$\log(I_3/I_2)$"); a_.set_title(f"{key}: pair states (t>=5)\ncontours: trapped basin, gamma 0.02 cyan / 0.05 green / 0.1 red", fontsize=8)
fig.tight_layout(); fig.savefig(f"{O}/fig_trap_scan.png", dpi=150)
json.dump({k: ({kk: vv for kk, vv in v.items() if kk != "_pts"} if isinstance(v, dict) else v) for k, v in res.items()}, open(f"{O}/trap_scan.json", "w"), indent=1)
