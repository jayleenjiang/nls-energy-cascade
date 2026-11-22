#!/usr/bin/env python3
"""Section 3.1 figures and tables from the v2 (unified Cartesian integrator) runs.
Same layout and style as paper_figures/scripts/plot_action_profiles.py and plot_collapse_bc1_bc2.py.
Inputs: results/{atlas,closed,dense,dtcheck}/*_summary.txt (S31 env var); replicates merged stream-weighted.
Outputs (OUT env var): Figures/*.pdf|png, data/*.csv
"""
import csv
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from load31 import load, merge

OUT = Path(os.environ.get("OUT", "out31"))
FIG, DATA = OUT / "Figures", OUT / "data"
FIG.mkdir(parents=True, exist_ok=True); DATA.mkdir(parents=True, exist_ok=True)
BCS = ("BC1", "BC2", "BC3", "BC3b"); TEMPS = ("T10_T5", "T4_T6", "T6_T6"); NS = (25, 50, 100)
N_STYLE = {25: ("#0072B2", "o", "-"), 50: ("#E69F00", "s", "--"), 100: ("#009E73", "^", "-.")}
DENSE_STYLE = {20: ("#3333CC", "o", "-"), 50: ("#1F85C6", "s", "--"), 100: ("#2CA25F", "^", "-."), 150: ("#E6550D", "D", ":")}
TEMP_TITLES = {"T10_T5": r"$(T_1,T_n)=(10,5)$", "T4_T6": r"$(T_1,T_n)=(4,6)$", "T6_T6": r"$(T_1,T_n)=(6,6)$"}
TV = {"T10_T5": (10, 5), "T4_T6": (4, 6), "T6_T6": (6, 6)}

mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
                     "mathtext.fontset": "dejavusans", "font.size": 8.0, "axes.labelsize": 8.5, "axes.titlesize": 9.0,
                     "legend.fontsize": 7.4, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.linewidth": 0.75,
                     "lines.linewidth": 1.35, "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600})

A, C, Dn, Dt = load("atlas"), load("closed"), load("dense"), load("dtcheck")


def prof(group, key):
    m = merge(group[key])
    if m is None:
        raise RuntimeError(key)
    n = len(m["mean_I"])
    return (np.arange(n) / (n - 1), m["mean_I"], m["se_I"], m)


def clean(ax):
    ax.grid(axis="y", color="0.90", linewidth=0.45); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)


def label(ax, s):
    ax.text(0.018, 0.955, s, transform=ax.transAxes, ha="left", va="top", fontweight="bold")


def save(fig, stem):
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight"); fig.savefig(FIG / f"{stem}.png", bbox_inches="tight", dpi=600); plt.close(fig)


def curve(ax, x, y, se, n, style, factor=1.0, every=None):
    color, marker, ls = style
    ax.fill_between(x, factor * (y - 1.96 * se), factor * (y + 1.96 * se), color=color, alpha=0.13, linewidth=0)
    ax.plot(x, factor * y, color=color, marker=marker, linestyle=ls, markevery=every or max(1, n // 9), markersize=2.4,
            markerfacecolor="white", markeredgewidth=0.55, label=fr"$n={n}$")


# ---------------------------------------------------------------- atlas (raw and rescaled)
def atlas(rescaled):
    fig, axes = plt.subplots(4, 3, figsize=(6.9, 7.7), sharex=True, sharey="row")
    for r, bc in enumerate(BCS):
        for c, temp in enumerate(TEMPS):
            ax = axes[r, c]
            for n in NS:
                x, y, se, _ = prof(A, f"atlas_{bc}_{temp}_n{n}")
                curve(ax, x, y, se, n, N_STYLE[n], np.sqrt(n) if rescaled else 1.0)
            if rescaled and temp == "T6_T6" and bc in ("BC1", "BC2"):
                ax.axhline(np.sqrt(6), color="0.25", linewidth=0.9, linestyle=":")
            if r == 0: ax.set_title(TEMP_TITLES[temp])
            if c == 0: ax.set_ylabel(f"{bc}\n" + (r"$\sqrt{n}\,\langle I_j\rangle$" if rescaled else r"$\langle I_j\rangle$"))
            if r == 3: ax.set_xlabel(r"normalized site $x=(j-1)/(n-1)$")
            label(ax, f"({chr(ord('a') + 3 * r + c)})"); ax.set_xlim(0, 1); clean(ax)
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.54, 0.999), ncol=3, frameon=False)
    fig.subplots_adjust(left=0.115, right=0.995, bottom=0.065, top=0.935, hspace=0.25, wspace=0.22)
    save(fig, "open_action_profile_collapse_atlas" if rescaled else "open_action_profile_atlas")


# ---------------------------------------------------------------- dense collapse BC1/BC2 at (4,6)
def collapse_dense():
    fig, axes = plt.subplots(1, 2, figsize=(6.9, 3.0), sharey=True)
    for c, bc in enumerate(("BC1", "BC2")):
        ax = axes[c]
        for n in (20, 50, 100, 150):
            x, y, se, _ = prof(Dn, f"dense_{bc}_T4_T6_n{n}")
            curve(ax, x, y, se, n, DENSE_STYLE[n], np.sqrt(n), every=max(1, n // 19))
        ax.set_title(bc); ax.set_xlabel(r"normalized site $x=(j-1)/(n-1)$"); ax.set_xlim(0, 1); clean(ax)
        label(ax, f"({chr(ord('a') + c)})")
    axes[0].set_ylabel(r"rescaled action $\sqrt{n}\,\langle I_j\rangle$")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.54, 1.02), ncol=4, frameon=False)
    fig.suptitle(r"Profile collapse at $(T_1,T_n)=(4,6)$", y=1.1)
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.17, top=0.86, wspace=0.12)
    save(fig, "open_action_profile_collapse_BC1_BC2_T4_T6_dense")


# ---------------------------------------------------------------- dense spatial means and fits
def scaling_dense():
    rows, fits = [], []
    fig, ax = plt.subplots(figsize=(3.6, 3.0))
    sty = {"BC1": ("#0072B2", "o"), "BC2": ("#E69F00", "s"), "BC3": ("#009E73", "^"), "BC3b": ("#CC79A7", "D")}
    for bc in BCS:
        ns, I, Ise = [], [], []
        for n in range(20, 151, 10):
            m = merge(Dn[f"dense_{bc}_T4_T6_n{n}"])
            ns.append(n); I.append(m["mean_I"].mean()); Ise.append(np.sqrt(np.sum(m["se_I"] ** 2)) / n)
            rows.append({"bc_label": bc, "T1": 4, "Tn": 6, "n": n, "streams": m["streams"], "spatial_mean_I": I[-1]})
        ns, I = np.array(ns, float), np.array(I)
        a, lc = np.polyfit(np.log(ns), np.log(I), 1)
        fit = lc + a * np.log(ns); r2 = 1 - np.sum((np.log(I) - fit) ** 2) / np.sum((np.log(I) - np.log(I).mean()) ** 2)
        fits.append({"bc_label": bc, "C": float(np.exp(lc)), "alpha": float(a), "loglog_R2": float(r2), "n_points": len(ns)})
        col, mk = sty[bc]
        ax.plot(ns, I, marker=mk, color=col, linestyle="none", markersize=3.2, markerfacecolor="white", markeredgewidth=0.7, label=bc)
        ax.plot(ns, np.exp(fit), color=col, linewidth=0.9, alpha=0.8)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(r"chain length $n$"); ax.set_ylabel(r"$\overline{I}_n$")
    ax.set_title(r"$(T_1,T_n)=(4,6)$, $n=20,\ldots,150$"); ax.legend(frameon=False, ncol=2); clean(ax)
    fig.subplots_adjust(left=0.17, right=0.98, bottom=0.17, top=0.90)
    save(fig, "spatial_mean_scaling_T4_T6_n20_150")
    with open(DATA / "spatial_mean_T4_T6_n20_150.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    with open(DATA / "scaling_fit_T4_T6_n20_150.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(fits[0])); w.writeheader(); w.writerows(fits)
    return fits


# ---------------------------------------------------------------- closed chain vs open (T,T)
def closed():
    fig, axes = plt.subplots(2, 2, figsize=(6.9, 4.8), sharex=True)
    rows = []
    for r, T in enumerate((6, 10)):
        for c, rescaled in enumerate((False, True)):
            ax = axes[r, c]
            for n in NS:
                x, y, se, m = prof(C, f"closed_BC1_T{T}_T0_n{n}")
                f = np.sqrt(n) if rescaled else 1.0
                curve(ax, x, y, se, n, N_STYLE[n], f)
                if T == 6:
                    xo, yo, seo, mo = prof(A, f"atlas_BC1_T6_T6_n{n}")
                    ax.plot(xo, f * yo, color="0.15", linewidth=0.6, linestyle="-", alpha=0.7,
                            label="open $(6,6)$" if n == 25 else None)
                    z = (y - yo) / np.hypot(se, seo)
                    if not rescaled:
                        rows.append({"T_left": T, "n": n, "closed_sqrt_n_spatial_mean": np.sqrt(n) * y.mean(),
                                     "open_TT_sqrt_n_spatial_mean": np.sqrt(n) * yo.mean(), "sqrt_T": np.sqrt(T),
                                     "pointwise_z_rms": float(np.sqrt(np.mean(z * z))), "pointwise_max_abs_z": float(np.abs(z).max()),
                                     "closed_mean_current": float(m["mean_current"].mean())})
                elif not rescaled:
                    rows.append({"T_left": T, "n": n, "closed_sqrt_n_spatial_mean": np.sqrt(n) * y.mean(),
                                 "open_TT_sqrt_n_spatial_mean": np.nan, "sqrt_T": np.sqrt(T), "pointwise_z_rms": np.nan,
                                 "pointwise_max_abs_z": np.nan, "closed_mean_current": float(m["mean_current"].mean())})
            if rescaled:
                ax.axhline(np.sqrt(T), color="0.25", linewidth=0.9, linestyle=":", label=r"$\sqrt{T_{\rm left}}$")
            ax.set_xlim(0, 1)
            ax.set_title(fr"$T_{{\rm left}}={T}$: " + (r"$\sqrt{n}$ rescaled" if rescaled else "raw"))
            ax.set_ylabel(r"$\sqrt{n}\,\langle I_j\rangle$" if rescaled else r"$\langle I_j\rangle$")
            if r == 1: ax.set_xlabel(r"normalized site $x=(j-1)/(n-1)$")
            ax.legend(frameon=False, loc="best"); label(ax, f"({chr(ord('a') + 2 * r + c)})"); clean(ax)
    fig.subplots_adjust(left=0.10, right=0.99, bottom=0.10, top=0.93, hspace=0.35, wspace=0.27)
    save(fig, "closed_action_profile_stationarity")
    with open(DATA / "closed_vs_open_gibbs.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


# ---------------------------------------------------------------- tables: atlas spatial means / slopes, dt check, merged profiles
def tables():
    rows, sc = [], []
    for bc in BCS:
        for temp in TEMPS:
            means = []
            for n in NS:
                x, y, se, m = prof(A, f"atlas_{bc}_{temp}_n{n}")
                means.append(y.mean())
                rows.append({"bc_label": bc, "temp_label": temp, "T1": TV[temp][0], "Tn": TV[temp][1], "n": n,
                             "streams": m["streams"], "mean_I_spatial": y.mean(), "mean_I_left": y[0], "mean_I_right": y[-1],
                             "se_I_left": se[0], "se_I_right": se[-1], "mean_current": float(m["mean_current"].mean())})
                with open(DATA / f"profile_{bc}_{temp}_n{n}.csv", "w", newline="") as f:
                    w = csv.writer(f); w.writerow(["j", "mean_I", "se_mean_I", "mean_current_bond_j"])
                    J = np.concatenate([[np.nan], m["mean_current"]])
                    for j in range(n): w.writerow([j + 1, y[j], se[j], J[j]])
            a = np.polyfit(np.log(NS), np.log(means), 1)[0]
            sc.append({"bc_label": bc, "temp_label": temp, "spatial_mean_I_values": ";".join(f"{v:.6f}" for v in means),
                       "loglog_slope_spatial_mean_I": float(a), "reference_slope": -0.5 if bc in ("BC1", "BC2") else 0.0})
    for name, data in (("logical_run_summary.csv", rows), ("spatial_mean_scaling.csv", sc)):
        with open(DATA / name, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0])); w.writeheader(); w.writerows(data)
    dt = []
    for bc in ("BC1", "BC3", "BC3b"):
        d = merge(Dt[f"dtcheck_{bc}_T10_T5_n100_dt5e-4"]); a = merge(A[f"atlas_{bc}_T10_T5_n100"])
        z = (d["mean_I"] - a["mean_I"]) / np.hypot(d["se_I"], a["se_I"])
        dt.append({"bc_label": bc, "n": 100, "T1": 10, "Tn": 5, "spatial_mean_dt5e-4": d["mean_I"].mean(), "spatial_mean_dt1e-3": a["mean_I"].mean(),
                   "relative_difference": d["mean_I"].mean() / a["mean_I"].mean() - 1, "pointwise_z_rms": float(np.sqrt(np.mean(z * z))),
                   "pointwise_max_abs_z": float(np.abs(z).max())})
    with open(DATA / "dt_check_T10_T5_n100.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(dt[0])); w.writeheader(); w.writerows(dt)
    for group, name in ((C, "closed"), (Dn, "dense")):
        for key in group:
            m = merge(group[key]); n = len(m["mean_I"])
            with open(DATA / f"profile_{key}.csv", "w", newline="") as f:
                w = csv.writer(f); w.writerow(["j", "mean_I", "se_mean_I", "mean_current_bond_j"])
                J = np.concatenate([[np.nan], m["mean_current"]])
                for j in range(n): w.writerow([j + 1, m["mean_I"][j], m["se_I"][j], J[j]])
    return sc, dt


if __name__ == "__main__":
    atlas(False); atlas(True); collapse_dense(); fits = scaling_dense(); closed(); sc, dt = tables()
    for r in sc: print(r["bc_label"], r["temp_label"], "slope %+.4f" % r["loglog_slope_spatial_mean_I"], r["spatial_mean_I_values"])
    for r in fits: print(r)
    for r in dt: print(r)
