"""Action-angle runs: (i) boundary-angle densities (same layout as fig43_boundary_angles), (ii) n = 3 sweep of Lambda.
Usage: python3 aa_angles_sweep.py AARUNS OUTDIR"""
import json, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
A, O = sys.argv[1:3]


def summ(p):
    d = {}
    for line in open(p):
        k, *v = line.split()
        try:
            d[k] = [float(x) for x in v]
        except ValueError:
            d[k] = v
    return d


temps = ["2_8", "5_5", "8_2"]
rows = [("n3_hold2", 3, "n=3, I2 held"), ("n4_hold2", 4, "n=4, I2 held (case 1)"), ("n4_hold23", 4, "n=4, I2 and I3 held (case 2)")]
fig, ax = plt.subplots(3, 3, figsize=(13, 10), sharex=True)
for i, (case, n, title) in enumerate(rows):
    for j, T in enumerate(temps):
        D = 2 * n - 1
        a = np.fromfile(f"{A}/angles/{case}_T{T}_dump.f32", dtype=np.float32).reshape(4, -1, 16, D).reshape(-1, D)
        curves = [(-a[:, n], r"$\theta_1=-\theta_2$"), (a[:, n + 1], r"$\theta_3$")]
        if n == 4:
            curves.append((a[:, n + 2], r"$\theta_4$"))
        for x, lab in curves:
            h, e = np.histogram(np.remainder(x + np.pi, 2 * np.pi) - np.pi, bins=72, range=(-np.pi, np.pi), density=True)
            ax[i, j].plot(0.5 * (e[1:] + e[:-1]), h, label=lab)
        for v in (2 * np.pi / 3, -2 * np.pi / 3):
            ax[i, j].axvline(v, color="grey", ls=":", lw=1)
        ax[i, j].set_title(f"{title}\n(T1,Tn)=({T.replace('_', ',')})", fontsize=9)
        ax[i, j].legend(fontsize=7)
    ax[i, 0].set_ylabel("density")
for j in range(3):
    ax[2, j].set_xlabel(r"angle (from the held mode); dotted: $\pm2\pi/3$")
fig.suptitle("Boundary dynamics in action-angle variables (rk4a_chain: schemes4 RK4_A + Cartesian bath), gamma = 0.1", fontsize=10)
fig.tight_layout(); fig.savefig(f"{O}/fig43_boundary_angles_rk4a.png", dpi=150)

Ts = ["0.5", "1", "2", "5", "10"]
tab = {}
for a in Ts:
    for b in Ts:
        s = summ(f"{A}/sweep/n3_hold2_T{a}_{b}_summary.txt")
        lam = s["rate_2"]; I1 = s["meanI_1"][0]; I3 = s["meanI_3"][0]
        tab[f"{a}_{b}"] = {"Lambda": lam, "E": [-lam[0] / 4, lam[1] / 4], "Lambda_over_I1_plus_I3": lam[0] / (I1 + I3)}
L = np.array([[tab[f"{a}_{b}"]["Lambda"][0] for b in Ts] for a in Ts])
E = np.array([[tab[f"{a}_{b}"]["E"][0] for b in Ts] for a in Ts])
R = np.array([[tab[f"{a}_{b}"]["Lambda_over_I1_plus_I3"] for b in Ts] for a in Ts])
out = {"table": tab, "Lambda_min": float(L.min()), "Lambda_max": float(L.max()), "E_min": float(E.min()), "E_max": float(E.max()),
       "ratio_min": float(R.min()), "ratio_max": float(R.max()), "all_negative_E": bool((E < 0).all())}
print(json.dumps({k: v for k, v in out.items() if k != "table"}, indent=1))
json.dump(out, open(f"{O}/sweep_rk4a.json", "w"), indent=1)
