"""Figs. 6-7 of Sec. 3.3: action current against the temperature difference (linear response) and against n.

Usage: python3 make_figures.py [OUT]      (OUT: folder for the PDFs; default: this folder)
Reads ../results/tables/current_response.csv and current_length.csv (made from the runs by
../analysis/fig33_reproduce.py) and writes fig33_current_response.pdf and fig33_current_length.pdf.
Response: -E[J] at the centre bond with standard errors, and the weighted fit a dT + b dT^2 for each n.
Length: -E[J] against n on log-log axes with the 95% band, and a reference line of slope -2 (BC1, BC2) or -1.
"""
import sys, os, csv, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "..", "common"))
import paperstyle as S
import matplotlib.pyplot as plt
U = os.path.join(HERE, "..", "results", "tables"); OUT = sys.argv[1] if len(sys.argv) > 1 else HERE
def rows(f): return list(csv.DictReader(open(f)))
BCS = ["BC1", "BC2", "BC3", "BC3b"]; L = S.letters()
# ---- response
R = rows(f"{U}/current_response.csv")
fig, ax = plt.subplots(2, 2, figsize=(S.TEXTWIDTH, 4.4)); ax = ax.ravel()
for i, bc in enumerate(BCS):
    a = ax[i]
    for j, n in enumerate((15, 25, 35, 45)):
        r = [x for x in R if x["bc"] == bc and int(x["n"]) == n]
        dT = np.array([float(x["dT"]) for x in r]); y = np.array([float(x["minus_J_centre"]) for x in r]); se = np.array([float(x["se"]) for x in r])
        o = np.argsort(dT); dT, y, se = dT[o], y[o], se[o]
        A = np.column_stack([dT, dT**2]) / se[:, None]; coef = np.linalg.lstsq(A, y / se, rcond=None)[0]
        xx = np.linspace(0, dT.max(), 100); c = S.MATLAB[j]
        a.errorbar(dT, y, yerr=se, fmt="o", ms=2.5, mfc="none", mec=c, ecolor=c, elinewidth=0.6, capsize=1.5, mew=0.6)
        a.plot(xx, coef[0]*xx + coef[1]*xx**2, color=c, label=f"$n={n}$")
    a.set_xlim(0, 4.1); a.set_ylim(bottom=0)
    S.panel_title(a, L[i], bc); a.set_xlabel(r"$\Delta T$"); a.set_ylabel(r"$-\mathbb{E}[J]$")
    a.legend(loc="upper left")
fig.tight_layout(w_pad=0.8, h_pad=0.8); fig.savefig(f"{OUT}/fig33_current_response.pdf"); plt.close(fig)
# ---- length
R = rows(f"{U}/current_length.csv")
fig, ax = plt.subplots(2, 2, figsize=(S.TEXTWIDTH, 4.4)); ax = ax.ravel()
for i, bc in enumerate(BCS):
    a = ax[i]; r = [x for x in R if x["bc"] == bc]
    n = np.array([int(x["n"]) for x in r]); y = np.array([float(x["minus_J_centre"]) for x in r]); se = np.array([float(x["se"]) for x in r])
    o = np.argsort(n); n, y, se = n[o], y[o], se[o]
    a.fill_between(n, y - 1.96*se, y + 1.96*se, color=S.MATLAB[0], alpha=0.25, lw=0)
    a.loglog(n, y, color=S.MATLAB[0], label=r"$-\mathbb{E}[J]$")
    p = -2 if bc in ("BC1", "BC2") else -1
    ref = y[0] * (n / n[0]) ** p * 0.8
    a.loglog(n, ref, "k--", lw=0.8, label=f"slope ${p}$")
    S.panel_title(a, L[i], bc); a.set_xlabel("$n$"); a.set_ylabel(r"$-\mathbb{E}[J]$")
    import matplotlib.ticker as mt
    a.set_xticks([20, 50, 100, 200]); a.xaxis.set_major_formatter(mt.ScalarFormatter()); a.xaxis.set_minor_formatter(mt.NullFormatter())
    a.yaxis.set_minor_formatter(mt.NullFormatter())
    if bc in ("BC3", "BC3b"):
        a.set_yticks([0.05, 0.1, 0.2]); a.yaxis.set_major_formatter(mt.FormatStrFormatter("%g"))
    a.legend(loc="upper right")
fig.tight_layout(w_pad=0.8, h_pad=0.8); fig.savefig(f"{OUT}/fig33_current_length.pdf"); plt.close(fig)
print("ok")
