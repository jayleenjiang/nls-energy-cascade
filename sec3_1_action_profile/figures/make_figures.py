"""Figures 2-4 of Sec. 3.1: action profiles of the open chain (atlas), their collapse, and the closed chain.

Usage: python3 make_figures.py [OUT]      (OUT: folder for the PDFs; default: this folder)
Reads the run summaries in ../results (through ../analysis/load31.py) and writes
open_action_profile_atlas.pdf, open_action_profile_collapse.pdf, closed_action_profile.pdf.
Replicates of the same case are merged (stream-weighted means); bands are 95% intervals (1.96 standard errors).
"""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "common"))          # paperstyle.py
os.environ["S31"] = os.path.join(HERE, "..", "results"); sys.path.insert(0, os.path.join(HERE, "..", "analysis"))
OUT = sys.argv[1] if len(sys.argv) > 1 else HERE
import paperstyle as S
import matplotlib.pyplot as plt
from load31 import load, merge
BCS = ("BC1", "BC2", "BC3", "BC3b"); TEMPS = ("T10_T5", "T4_T6", "T6_T6"); NS = (25, 50, 100)
TT = {"T10_T5": "(10,5)", "T4_T6": "(4,6)", "T6_T6": "(6,6)"}
A, C, Dn = load("atlas"), load("closed"), load("dense"); L = S.letters()
def prof(group, key):
    m = merge(group[key]); n = len(m["mean_I"])
    return np.arange(n) / (n - 1), m["mean_I"], m["se_I"], m
def curve(ax, x, y, se, color, f=1.0, label=None, band=True):
    if band: ax.fill_between(x, f * (y - 1.96 * se), f * (y + 1.96 * se), color=color, alpha=0.25, lw=0)
    ax.plot(x, f * y, color=color, label=label)
XL = r"$x=(j-1)/(n-1)$"
# atlas
fig, ax = plt.subplots(4, 3, figsize=(S.TEXTWIDTH, 7.6), sharey="row"); k = 0
for r, bc in enumerate(BCS):
    for c, temp in enumerate(TEMPS):
        a = ax[r, c]
        for i, n in enumerate(NS):
            x, y, se, _ = prof(A, f"atlas_{bc}_{temp}_n{n}"); curve(a, x, y, se, S.MATLAB[i], label=f"$n={n}$")
        a.set_xlim(0, 1); S.panel_title(a, L[k], f"{bc}, $(T_1,T_n)={TT[temp]}$"); k += 1
        a.set_xlabel(XL); a.set_ylabel(r"$\langle I_j\rangle$"); a.tick_params(labelleft=True)
fig.tight_layout(w_pad=0.5, h_pad=0.6); S.top_legend(fig, ax[0, 0], 3); fig.savefig(f"{OUT}/open_action_profile_atlas.pdf"); plt.close(fig)
# collapse
fig, ax = plt.subplots(1, 2, figsize=(S.TEXTWIDTH, 2.4))
for c, bc in enumerate(("BC1", "BC2")):
    a = ax[c]
    for i, n in enumerate((20, 50, 100, 150)):
        x, y, se, _ = prof(Dn, f"dense_{bc}_T4_T6_n{n}"); curve(a, x, y, se, S.MATLAB[i], np.sqrt(n), label=f"$n={n}$", band=False)
    a.set_xlim(0, 1); S.panel_title(a, L[c], bc); a.set_xlabel(XL); a.set_ylabel(r"$\sqrt{n}\,\langle I_j\rangle$")
ax[0].legend(loc="lower right")
fig.tight_layout(w_pad=0.8); fig.savefig(f"{OUT}/open_action_profile_collapse.pdf"); plt.close(fig)
# closed
fig, ax = plt.subplots(2, 2, figsize=(S.TEXTWIDTH, 4.4)); k = 0
for r, T in enumerate((6, 10)):
    for c, resc in enumerate((False, True)):
        a = ax[r, c]
        for i, n in enumerate(NS):
            x, y, se, m = prof(C, f"closed_BC1_T{T}_T0_n{n}"); f = np.sqrt(n) if resc else 1.0
            curve(a, x, y, se, S.MATLAB[i], f, label=f"$n={n}$", band=False)
            if T == 6:
                xo, yo, seo, mo = prof(A, f"atlas_BC1_T6_T6_n{n}")
                a.plot(xo, f * yo, color="k", lw=0.5, label="open $(6,6)$" if n == 25 else None)
        if resc: a.axhline(np.sqrt(T), color="k", ls=":", lw=0.7, label=r"$\sqrt{T_{\rm left}}$")
        a.set_xlim(0, 1); a.set_xlabel(XL)
        a.set_ylabel(r"$\sqrt{n}\,\langle I_j\rangle$" if resc else r"$\langle I_j\rangle$")
        S.panel_title(a, L[k], f"$T_{{\\rm left}}={T}$"); k += 1
fig.tight_layout(w_pad=0.8, h_pad=0.8)
h0, l0 = ax[0, 0].get_legend_handles_labels(); h1, l1 = ax[0, 1].get_legend_handles_labels()
H = dict(zip(l0, h0)); H.update(dict(zip(l1, h1)))
keys = ["$n=25$", "$n=50$", "$n=100$", "open $(6,6)$", r"$\sqrt{T_{\rm left}}$"]
fig.legend([H[k] for k in keys], keys, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=5)
fig.savefig(f"{OUT}/closed_action_profile.pdf"); plt.close(fig)
print("ok")
