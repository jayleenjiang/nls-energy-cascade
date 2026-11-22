"""Figs. 9-12 of Section 4.1, drawn from the small arrays in ../results/figure_data/ (made by ../analysis/reduce_figure_data.py).

  fig_marginals_n3.pdf, fig_marginals_n4.pdf         two-dimensional marginals: simulation, network, log ratio, noise level
  fig_theta_slices_n3.pdf, fig_theta_slices_n4.pdf   NESS (network) and Gibbs densities in the angles at equal actions
Usage: python3 make_figures.py [output folder]   (default: this folder).  No TensorFlow needed.
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


def marginals(n):
    """3 rows (pairs of coordinates) x 4 columns: simulation, network, log(network/simulation), noise level"""
    d = np.load(os.path.join(DATA, f"marginals_n{n}.npz"))
    if n == 3:
        lab = ["$I_1$", "$I_2$", "$I_3$", r"$\theta_1$", r"$\theta_3$"]
    else:
        lab = ["$I_1$", "$I_2$", "$I_3$", "$I_4$", r"$\theta_2$", r"$\theta_3$", r"$\theta_4$"]
    heads = ["simulation", "network", "log(network/simulation)", "noise level"]
    fig, ax = plt.subplots(3, 4, figsize=(S.TEXTWIDTH, 4.9)); L = S.letters(); k = 0
    for r, (a, b) in enumerate(d["rows"]):
        ex, ey = d[f"r{r}_edges_x"], d[f"r{r}_edges_y"]
        A_, B_ = d[f"r{r}_simA"], d[f"r{r}_simB"]
        cnt = A_ + B_; T = cnt / cnt.sum(); mask = cnt < 2000           # bins with fewer than 2000 counts are blank
        E = d[f"r{r}_net"]; E = E / E.sum(); IA, IB = d[f"r{r}_netA"], d[f"r{r}_netB"]
        with np.errstate(all="ignore"):
            Z = np.where(mask, np.nan, np.log(E / T))
            # noise level: half the log ratio between the two validation sets plus half that between the two network halves
            N = np.where(mask, np.nan, np.log((A_ / A_.sum()) / (B_ / B_.sum())) / 2 + np.log((IA / IA.sum()) / (IB / IB.sum())) / 2)
            lt = np.log10(np.where(cnt > 0, T, np.nan)); le = np.log10(np.where(E > 0, E, np.nan))

        def wrms(Q):
            ok = ~np.isnan(Q); return np.sqrt(np.sum(T[ok] * Q[ok] ** 2) / np.sum(T[ok]))
        print(n, lab[a], lab[b], "rms", round(wrms(Z), 3), round(wrms(N), 3))
        vmin = np.nanpercentile(lt, 5); vmax = np.nanmax(lt)
        for c, (Q, lo, hi) in enumerate([(lt, vmin, vmax), (le, vmin, vmax), (Z, -0.2, 0.2), (N, -0.2, 0.2)]):
            x = ax[r, c]
            im = x.pcolormesh(ex, ey, Q.T, cmap=S.CMAP, vmin=lo, vmax=hi, shading="flat", rasterized=False,
                              edgecolors="face", linewidth=0.1)
            S.colorbar(im, x); S.panel_title(x, L[k], heads[c]); k += 1
            x.set_xlabel(lab[a]); x.set_ylabel(lab[b])
            if a >= n: S.angle_axis(x, "x")
            if b >= n: S.angle_axis(x, "y")
    fig.tight_layout(w_pad=1.0, h_pad=0.6)
    fig.savefig(os.path.join(OUT, f"fig_marginals_n{n}.pdf")); plt.close(fig)


def slices(n):
    """row 1: NESS network density, row 2: Gibbs density at T = 5; columns I = 0.3, 0.5, 1"""
    d = np.load(os.path.join(DATA, "theta_slices.npz")); g = d["grid"]
    fig, ax = plt.subplots(2, 3, figsize=(S.TEXTWIDTH, 4.3)); L = S.letters(); k = 0
    for r in range(2):
        for c, I in enumerate((0.3, 0.5, 1.0)):
            p = d[f"n{n}_{'nn' if r == 0 else 'gibbs'}_I{I}"]; a = ax[r, c]
            im = a.pcolormesh(g, g, p.T, cmap=S.CMAP, shading="gouraud", rasterized=True)
            S.colorbar(im, a); a.set_aspect("equal"); S.angle_axis(a)
            S.panel_title(a, L[k], f"{'NESS' if r == 0 else 'equilibrium'}, $I={I:g}$"); k += 1
            a.set_xlabel(r"$\theta_1$" if n == 3 else r"$\theta_2$"); a.set_ylabel(r"$\theta_3$")
    fig.tight_layout(w_pad=0.6, h_pad=0.8)
    fig.savefig(os.path.join(OUT, f"fig_theta_slices_n{n}.pdf")); plt.close(fig)


for n in (3, 4):
    marginals(n)
    slices(n)
print("written to", OUT)
