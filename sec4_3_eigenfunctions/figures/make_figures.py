"""Figs. 17-19 of Section 4.3, drawn from ../results/figure_data/ (made by ../analysis/reduce_figure_data.py).

  fig_psi_MH.pdf               psi_1, psi_2, psi_3 averaged over states with given (M, H); NESS and equilibrium
  fig_psi1_theta.pdf           psi_1 in the angles at I_1 = I_2 = I_3 = I
  fig_loc_backward_gscan.pdf   E_x[g(X_t)] from two initial states, five values of gamma
Usage: python3 make_figures.py [output folder]   (default: this folder)
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "common"))
import paperstyle as S
import matplotlib.pyplot as plt

DATA = os.path.join(HERE, "..", "results", "figure_data")
OUT = sys.argv[1] if len(sys.argv) > 1 else HERE
CASES = [("neq_dt1e-3", "NESS"), ("eq_dt1e-3", "equilibrium")]


def psi_MH():
    d = np.load(os.path.join(DATA, "psi_MH.npz")); L = S.letters(); k = 0
    fig, ax = plt.subplots(2, 3, figsize=(S.TEXTWIDTH, 4.0))
    for r, (case, lab) in enumerate(CASES):
        mb, hb, cnt = d[f"{case}_Medges"], d[f"{case}_Hedges"], d[f"{case}_count"]
        for c, j in enumerate((1, 2, 3)):
            with np.errstate(invalid="ignore"):
                mean = np.where(cnt >= 30, d[f"{case}_sum_psi{j}"] / cnt, np.nan)   # bins with < 30 samples blank
            v = np.nanmax(np.abs(mean)); a = ax[r, c]
            im = a.pcolormesh(mb, hb, mean.T, cmap=S.CMAP, vmin=-v, vmax=v, shading="flat", rasterized=False,
                              edgecolors="face", linewidth=0.1)
            S.colorbar(im, a)
            S.panel_title(a, L[k], f"{lab}, $\\psi_{j}$"); k += 1
            a.set_xlabel("$M$"); a.set_ylabel("$H$")
    fig.tight_layout(w_pad=0.6, h_pad=0.8)
    fig.savefig(os.path.join(OUT, "fig_psi_MH.pdf")); plt.close(fig)


def psi1_theta():
    d = np.load(os.path.join(DATA, "psi1_theta.npz")); th = d["grid"]; L = S.letters(); k = 0
    fig, ax = plt.subplots(2, 3, figsize=(S.TEXTWIDTH, 4.3))
    for r, (case, lab) in enumerate(CASES):
        for c, I in enumerate((1.0, 1.5, 2.0)):
            psi = d[f"{case}_I{I:g}"]; v = np.abs(psi).max(); a = ax[r, c]
            im = a.pcolormesh(th, th, psi, cmap=S.CMAP, vmin=-v, vmax=v, shading="gouraud", rasterized=True)
            S.colorbar(im, a); a.set_aspect("equal"); S.angle_axis(a)
            S.panel_title(a, L[k], f"{lab}, $I={I:g}$"); k += 1
            a.set_xlabel(r"$\theta_1$"); a.set_ylabel(r"$\theta_3$")
    fig.tight_layout(w_pad=0.6, h_pad=0.8)
    fig.savefig(os.path.join(OUT, "fig_psi1_theta.pdf")); plt.close(fig)


def loc_backward():
    d = np.load(os.path.join(DATA, "return_probabilities.npz")); NF = 12   # columns: 12 ball indicators, then test functions
    G = [("0.1", r"$\gamma=0.1$"), ("0.05", r"$\gamma=0.05$"), ("0.02", r"$\gamma=0.02$"), ("0.01", r"$\gamma=0.01$"),
         ("0.001", r"$\gamma=0.001$")]
    fig, ax = plt.subplots(2, 3, figsize=(S.TEXTWIDTH, 4.2)); L = S.letters(); kk = 0
    # c0: x = (1,1,1,0,0) (center 0), c2: x = (2,1,1,0,0) (center 2); dotted lines: return times of the Hamiltonian orbit
    for i, (c, ci, x0, trs) in enumerate((("c0", 0, "(1,1,1,0,0)", [0.65 * k for k in range(1, 8)]),
                                          ("c2", 2, "(2,1,1,0,0)", [0.50, 1.03, 1.55, 2.06, 2.56, 3.08]))):
        for j, (k, ttl, log) in enumerate(((2 * ci, "ball, $r=0.5$", True), (2 * ci + 1, "ball, $r=1$", True),
                                           (NF + 3, r"$g=\cos\theta_1$", False))):
            a = ax[i, j]
            for q, (g, lb) in enumerate(G):
                n = f"ex_{c}_neq_g{g}"; t, m, se = d[n + "_t"], d[n + "_mean"], d[n + "_se"]
                a.fill_between(t, m[:, k] - 2 * se[:, k], m[:, k] + 2 * se[:, k], color=S.MATLAB[q], alpha=.25, lw=0)
                a.plot(t, m[:, k], color=S.MATLAB[q], lw=0.8, label=lb)
            for tr in trs:
                a.axvline(tr, color="k", ls=":", lw=0.5)
            if log:
                a.set_yscale("log"); a.set_ylim(1e-4, 1.2)
            a.set_xlim(0, 5); a.set_xlabel("$t$"); a.set_ylabel(r"$\mathbb{E}_x[g(X_t)]$")
            S.panel_title(a, L[kk], f"$x={x0}$\n{ttl}"); kk += 1
    fig.tight_layout(w_pad=0.5, h_pad=0.8); S.top_legend(fig, ax[0, 0], 5)
    fig.savefig(os.path.join(OUT, "fig_loc_backward_gscan.pdf")); plt.close(fig)


psi_MH(); psi1_theta(); loc_backward()
print("written to", OUT)
