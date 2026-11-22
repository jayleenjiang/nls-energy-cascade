"""Sec. 3.3 figures from batch_09 (runs on the Mac after download, not on Unity).
Usage: python3 fig33_reproduce.py RESULTS_DIR OUT_DIR
  RESULTS_DIR = batch_09 results/ (subfolders resp/ and len/), OUT_DIR = new folder for figures and tables.
Current: bond at the centre of the chain, j = floor(n/2) + 1 (bond between modes n/2 and n/2 + 1), as in the old Sec. 3.3.
Sign: T1 < Tn in every run, so action flows from the right (hot) end to the left; the figures show -E[J] > 0.
Error: standard error over streams from <label>_current_streams.txt; band = 95% confidence interval (1.96 SE).
Streams with a non-finite value, or with time-averaged M above 3 x the median of the run (runaway, Sec. 3.1), are
excluded and counted in the tables.
Outputs: fig_current_response.png (3 panels BC1-BC3, layout of the old Fig. 6), fig_current_length.png (old Fig. 7),
the same with BC3b as a fourth panel (*_4bc.png), current_response.csv, current_length.csv, fits_length.txt."""
import os, sys, glob, csv
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
R, O = sys.argv[1:3]; os.makedirs(O, exist_ok=True)
plt.rcParams.update({"font.size": 12})
MATLAB = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"]


def read_run(prefix):
    d = {l.split()[0]: l.split()[1:] for l in open(prefix + "_summary.txt") if l.strip()}
    n = int(d["n"][0]); L = np.loadtxt(prefix + "_current_streams.txt", ndmin=2)
    M, J = L[:, 0], L[:, 1:]
    ok = np.isfinite(L).all(1)
    ok &= M <= 3 * np.median(M[ok]) if ok.any() else ok
    jc = n // 2                               # column of bond floor(n/2)+1 (columns are bonds 2..n)
    x = -J[ok, jc - 1]
    return {"n": n, "bc": d["bc"][0], "T1": float(d["T1"][0]), "Tn": float(d["Tn"][0]), "J": x.mean(),
            "se": x.std(ddof=1) / np.sqrt(ok.sum()), "used": int(ok.sum()), "excluded": int((~ok).sum()),
            "nonfinite_flag": int(d["nonfinite"][0]), "midpoint_failures": int(d["midpoint_failures"][0])}


def load(group):
    return [read_run(f[:-12]) for f in sorted(glob.glob(f"{R}/{group}/*_summary.txt"))]


resp, length = load("resp"), load("len")
for name, rows, keys in (("current_response.csv", resp, ("bc", "n", "T1", "Tn")), ("current_length.csv", length, ("bc", "n", "T1", "Tn"))):
    with open(f"{O}/{name}", "w", newline="") as f:
        w = csv.writer(f); w.writerow(list(keys) + ["dT", "minus_J_centre", "se", "streams_used", "streams_excluded", "nonfinite_flag", "midpoint_failures"])
        for r in sorted(rows, key=lambda r: (r["bc"], r["n"], r["Tn"] - r["T1"])):
            w.writerow([r[k] for k in keys] + [f"{r['Tn'] - r['T1']:g}", f"{r['J']:.6e}", f"{r['se']:.3e}", r["used"], r["excluded"], r["nonfinite_flag"], r["midpoint_failures"]])

for bcs, tag in ((("BC1", "BC2", "BC3"), ""), (("BC1", "BC2", "BC3", "BC3b"), "_4bc")):
    if resp:
        fig, ax = plt.subplots(1, len(bcs), figsize=(5.2 * len(bcs), 4.2))
        for i, bc in enumerate(bcs):
            ns = sorted({r["n"] for r in resp if r["bc"] == bc})
            for k, n in enumerate(ns):
                rr = sorted([r for r in resp if r["bc"] == bc and r["n"] == n], key=lambda r: r["Tn"] - r["T1"])
                ax[i].plot([r["Tn"] - r["T1"] for r in rr], [r["J"] for r in rr], color=MATLAB[k % 4], lw=1.8, label=f"n = {n}")
            ax[i].set_title(f"Boundary condition {bc[2:]}", fontweight="bold"); ax[i].set_xlabel(r"$\Delta T$"); ax[i].set_ylabel("Action current")
            ax[i].set_xlim(0, 4); ax[i].set_ylim(bottom=0); ax[i].legend(loc="upper left")
        fig.tight_layout(); fig.savefig(f"{O}/fig_current_response{tag}.png", dpi=200, facecolor="white"); plt.close(fig)
    if length:
        fig, ax = plt.subplots(1, len(bcs), figsize=(5.2 * len(bcs), 4.2))
        for i, bc in enumerate(bcs):
            rr = sorted([r for r in length if r["bc"] == bc], key=lambda r: r["n"])
            if not rr: continue
            n = np.array([r["n"] for r in rr], float); J = np.array([r["J"] for r in rr]); ci = 1.96 * np.array([r["se"] for r in rr])
            sl = -2 if bc in ("BC1", "BC2") else -1
            ax[i].loglog(n, J, color="blue", lw=1.8, label="Action current")
            ax[i].loglog(n, J + ci, color="#ff8080", lw=1.2, label="Confidence interval (+)")
            ax[i].loglog(n, np.clip(J - ci, 1e-12, None), color="#ff8080", lw=1.2, label="Confidence interval (-)")
            ax[i].loglog(n, 0.5 * J[0] * (n / n[0]) ** sl, color="black", lw=1.5, label=f"slope = {sl} ref")
            ax[i].set_xlim(n.min(), n.max()); ax[i].set_title(f"Boundary condition {bc[2:]}", fontweight="bold")
            ax[i].set_xlabel("Length"); ax[i].set_ylabel("Action current"); ax[i].legend(loc="upper right", fontsize=10)
        fig.tight_layout(); fig.savefig(f"{O}/fig_current_length{tag}.png", dpi=200, facecolor="white"); plt.close(fig)

with open(f"{O}/fits_length.txt", "w") as f:
    for bc in ("BC1", "BC2", "BC3", "BC3b"):
        rr = sorted([r for r in length if r["bc"] == bc], key=lambda r: r["n"])
        if len(rr) < 3: continue
        n = np.array([r["n"] for r in rr], float); J = np.array([r["J"] for r in rr]); se = np.array([r["se"] for r in rr])
        for lo in (20, 60, 100):
            m = n >= lo
            if m.sum() < 3 or (J[m] <= 0).any(): continue
            w = (J[m] / se[m]) ** 2                     # weights 1/var(log J)
            A = np.vstack([np.ones(m.sum()), np.log(n[m])]).T
            cov = np.linalg.inv(A.T @ (A * w[:, None])); b = cov @ (A.T @ (w * np.log(J[m])))
            f.write(f"{bc}: n >= {lo}: E[J] ~ n^({b[1]:.3f} +- {np.sqrt(cov[1, 1]):.3f})  ({m.sum()} lengths, weighted least squares in log-log)\n")
print(open(f"{O}/fits_length.txt").read())
