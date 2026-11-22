"""Figs. 13-16 of Section 4.2, drawn from ../results/figure_data/ (made by ../analysis/reduce_figure_data.py).

  fig43_boundary_angles.pdf       end-angle densities with the middle mode(s) held at 0
  fig43_worstcase_M.pdf           mean M(t) from the three initial phases, M0 = 200
  fig43_smallgamma_endpoints.pdf  median I_1 + I_n for gamma = 1e-3 ... 1e-6
  fig43_smallgamma_M.pdf          (E M(t) - M(0)) / gamma with a band of 2 standard errors
Usage: python3 make_figures.py [output folder]   (default: this folder)
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "common"))
import paperstyle as S
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

DATA = os.path.join(HERE, "..", "results", "figure_data")
OUT = sys.argv[1] if len(sys.argv) > 1 else HERE
load = lambda name: np.load(os.path.join(DATA, name))


def boundary_angles():
    d = load("boundary_angles.npz"); e = d["edges"]; x = 0.5 * (e[1:] + e[:-1])
    rows = [("n3_hold2", 3, r"$n=3$, $I_2=0$"), ("n4_hold2", 4, r"$n=4$, $I_2=0$"), ("n4_hold23", 4, r"$n=4$, $I_2=I_3=0$")]
    labels = {"theta1": r"$\theta_1$", "theta3": r"$\theta_3$", "theta4": r"$\theta_4$"}
    fig, ax = plt.subplots(3, 3, figsize=(S.TEXTWIDTH, 5.4)); L = S.letters(); k = 0
    for i, (case, n, title) in enumerate(rows):
        for j, T in enumerate(("2_8", "5_5", "8_2")):
            a = ax[i, j]
            for c in ("theta1", "theta3") + (("theta4",) if n == 4 else ()):
                a.plot(x, d[f"{case}_T{T}_{c}"], label=labels[c])
            for v in (2 * np.pi / 3, -2 * np.pi / 3):
                a.axvline(v, color="k", ls=":", lw=0.6)
            S.angle_axis(a, "x"); a.set_ylim(bottom=0)
            S.panel_title(a, L[k], f"{title}, $({T.replace('_', ',')})$"); k += 1
            a.legend(loc="upper right"); a.set_xlabel(r"$\theta$"); a.set_ylabel("density")
    fig.tight_layout(w_pad=0.5, h_pad=0.6); fig.savefig(os.path.join(OUT, "fig43_boundary_angles.pdf")); plt.close(fig)


def worstcase():
    Z = load("worstcase_M200.npz")
    LAB = {"worst": r"all angles $4\pi/3$", "anti": r"interior angles $\pi$", "random": "random phases"}
    fig, axs = plt.subplots(2, 2, figsize=(S.TEXTWIDTH, 4.2)); ax = axs.ravel(); L = S.letters()
    for i, n in enumerate((3, 4, 5, 6)):
        for j, ph in enumerate(("worst", "anti", "random")):
            k = f"n{n}_M200_{ph}"
            if f"{k}__t" in Z:                    # the opposite-phase start exists only for n >= 4
                ax[i].plot(np.r_[0, Z[f"{k}__t"]], np.r_[200, Z[f"{k}__mean"]], color=S.MATLAB[j], label=LAB[ph])
        ax[i].axhline(30, color="k", ls=":", lw=0.6)
        ax[i].set_xlim(0, 1100); ax[i].set_ylim(0, 205)
        S.panel_title(ax[i], L[i], f"$n={n}$"); ax[i].set_xlabel("$t$"); ax[i].set_ylabel(r"mean $M(t)$")
    ax[1].legend(loc="upper right")
    fig.tight_layout(w_pad=0.8, h_pad=0.8); fig.savefig(os.path.join(OUT, "fig43_worstcase_M.pdf")); plt.close(fig)


def smallgamma_endpoints():
    d = load("smallgamma_endpoints.npz")
    G = ["1e-3", "1e-4", "1e-5", "1e-6"]
    GL = [r"$\gamma=10^{-3}$", r"$\gamma=10^{-4}$", r"$\gamma=10^{-5}$", r"$\gamma=10^{-6}$"]
    fig, axs = plt.subplots(2, 2, figsize=(S.TEXTWIDTH, 4.2)); ax = axs.ravel(); L = S.letters()
    for i, n in enumerate((3, 4, 5, 6)):
        for k, g in enumerate(G):
            for ph, ls in (("worst", "-"), ("random", ":")):
                e = d[f"n{n}_g{g}_{ph}"]; t = np.arange(len(e)) * 1e-4
                ax[i].semilogy(t, e, ls, color=S.MATLAB[k], label=GL[k] if ph == "worst" else None)
        ax[i].axhline(0.02, color="k", lw=0.5)           # initial value I_1 + I_n = 0.02
        ax[i].set_xlim(0, 0.06); ax[i].set_ylim(2e-5, 200)
        S.panel_title(ax[i], L[i], f"$n={n}$"); ax[i].set_xlabel("$t$"); ax[i].set_ylabel(r"median $I_1+I_n$")
    h, l = ax[0].get_legend_handles_labels()
    h += [Line2D([], [], color="k", ls="-"), Line2D([], [], color="k", ls=":")]
    l += [r"all angles $4\pi/3$", "random phases"]
    fig.tight_layout(w_pad=0.8, h_pad=0.8)
    fig.legend(h, l, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, columnspacing=1.5)
    fig.savefig(os.path.join(OUT, "fig43_smallgamma_endpoints.pdf")); plt.close(fig)


def smallgamma_M():
    d = load("smallgamma_Mrise.npz"); L = S.letters(); k = 0; res = {}
    GL = {"1e-3": r"$\gamma=10^{-3}$", "1e-4": r"$\gamma=10^{-4}$"}
    fig, ax = plt.subplots(2, 4, figsize=(S.TEXTWIDTH, 3.3))
    for r, g in enumerate(("1e-3", "1e-4")):
        for c, n in enumerate((3, 4, 5, 6)):
            for j, (ph, lab) in enumerate((("worst", r"all angles $4\pi/3$"), ("random", "random phases"))):
                m, se = d[f"n{n}_g{g}_{ph}_mean"], d[f"n{n}_g{g}_{ph}_se"]
                t = (np.arange(len(m)) + 1) * 1e-3          # sampled every 1e-3, first sample at t = 1e-3
                ax[r, c].fill_between(t, m - 2 * se, m + 2 * se, color=S.MATLAB[j], alpha=0.25, lw=0)
                ax[r, c].plot(t, m, color=S.MATLAB[j], label=lab)
                q = int(np.argmax(m)); res[f"n{n}_g{g}_{ph}"] = (round(float(m[q]), 2), round(float(t[q]), 3))
            ax[r, c].axhline(0, color="k", lw=0.5); ax[r, c].set_xlim(0, 0.06); ax[r, c].set_ylim(-6, 3)
            S.panel_title(ax[r, c], L[k], f"$n={n}$, {GL[g]}"); k += 1
            ax[r, c].set_xlabel("$t$")
            if c == 0: ax[r, c].set_ylabel(r"$(\mathbb{E}M(t)-M(0))/\gamma$")
    fig.tight_layout(w_pad=0.3, h_pad=0.6); S.top_legend(fig, ax[0, 0], 2)
    fig.savefig(os.path.join(OUT, "fig43_smallgamma_M.pdf")); plt.close(fig)
    print("largest rise / gamma and its time:", res)


boundary_angles(); worstcase(); smallgamma_endpoints(); smallgamma_M()
print("written to", OUT)
