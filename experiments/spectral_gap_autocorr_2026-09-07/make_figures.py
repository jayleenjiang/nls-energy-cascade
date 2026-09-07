#!/usr/bin/env python3
"""Publication figures for the frozen n=3 stationary-autocorrelation analysis."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
ANALYSIS = ROOT / "analysis"
FIGURES = ROOT / "figures"
FIGURES.mkdir(exist_ok=True)

OBS = [
    "cos_theta1", "sin_theta1", "cos_theta3", "cos_theta1_minus_theta3",
    "I2", "I1_plus_I3", "I1_minus_I3",
]
LABELS = {
    "cos_theta1": r"$\cos\theta_1$",
    "sin_theta1": r"$\sin\theta_1$",
    "cos_theta3": r"$\cos\theta_3$",
    "cos_theta1_minus_theta3": r"$\cos(\theta_1-\theta_3)$",
    "I2": r"$I_2$",
    "I1_plus_I3": r"$I_1+I_3$",
    "I1_minus_I3": r"$I_1-I_3$",
}
CASES = [
    ("driven_T2_T8", r"Driven: $(T_1,T_3)=(2,8)$"),
    ("equilibrium_T5_T5", r"Equal temperature: $(T_1,T_3)=(5,5)$"),
]
COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]

mpl.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "stix",
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "legend.fontsize": 7.4,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.bbox": "tight",
})


def read_rows(path):
    groups = defaultdict(list)
    with path.open() as f:
        for row in csv.DictReader(f):
            clean = {"observable": row["observable"]}
            for k, v in row.items():
                if k == "observable":
                    continue
                try: clean[k] = float(v)
                except (TypeError, ValueError): clean[k] = np.nan
            groups[row["observable"]].append(clean)
    return groups


def read_plateaus(path):
    out = {}
    with path.open() as f:
        for row in csv.DictReader(f):
            try:
                out[row["observable"]] = (
                    float(row["plateau_start"]), float(row["plateau_end"]), float(row["lambda_R"])
                )
            except ValueError:
                out[row["observable"]] = None
    return out


fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.7), sharex="col")
for col, (case, title) in enumerate(CASES):
    groups = read_rows(ANALYSIS / f"{case}_autocorrelations.csv")
    plateaus = read_plateaus(ANALYSIS / f"{case}_plateaus.csv")
    axc, axr = axes[0, col], axes[1, col]
    for color, obs in zip(COLORS, OBS):
        rows = groups[obs]
        t = np.array([r["time"] for r in rows])
        c = np.array([r["C"] for r in rows])
        se = np.array([r["SE_C"] for r in rows])
        lam = np.array([r["lambda_eff"] for r in rows])
        lo = np.array([r["lambda_ci_low"] for r in rows])
        hi = np.array([r["lambda_ci_high"] for r in rows])
        c0 = c[0]
        pos = (t > 0) & (c > 0)
        neg = (t > 0) & (c < 0)
        axc.plot(t[pos], c[pos] / c0, color=color, lw=1.4, label=LABELS[obs])
        if np.any(neg):
            axc.plot(t[neg], -c[neg] / c0, color=color, lw=1.0, ls=":", alpha=0.9)
        good = (t > 0) & np.isfinite(lam) & np.isfinite(lo) & np.isfinite(hi)
        axr.plot(t[good], lam[good], color=color, lw=1.15)
        axr.fill_between(t[good], lo[good], hi[good], color=color, alpha=0.08, lw=0)
        p = plateaus[obs]
        if p is not None:
            start, end, rate = p
            axr.plot([start, end], [rate, rate], color=color, lw=4.0, solid_capstyle="butt")
    axc.set_title(title)
    axc.set_xscale("log")
    axc.set_yscale("log")
    axc.set_xlim(0.01, 10)
    axc.set_ylim(1e-5, 2)
    axc.grid(True, which="major", alpha=0.2)
    axr.set_xscale("log")
    axr.set_xlim(0.01, 10)
    axr.set_ylim(-10, 2)
    axr.axhline(0, color="0.5", lw=0.8)
    axr.axvline(1.0, color="0.55", lw=0.8, ls="--")
    axr.grid(True, which="major", alpha=0.2)
    axr.set_xlabel(r"lag $t$")

axes[0, 0].set_ylabel(r"$|C_g(t)|/C_g(0)$")
axes[1, 0].set_ylabel(r"$\lambda_{\rm eff}(t)$")
axes[0, 0].legend(ncol=2, frameon=False, loc="lower left")
axes[0, 1].text(
    0.03, 0.04, "solid: $C_g>0$\ndotted: $C_g<0$",
    transform=axes[0, 1].transAxes, va="bottom", fontsize=7.5,
)
axes[1, 1].text(
    0.03, 0.04, r"thick segments: accepted plateaus" "\n" r"dashed vertical: $t_{\min}=1$",
    transform=axes[1, 1].transAxes, va="bottom", fontsize=7.5,
)
fig.subplots_adjust(hspace=0.12, wspace=0.20)
for suffix in ("pdf", "png"):
    fig.savefig(FIGURES / f"stationary_autocorrelation_rates.{suffix}", dpi=500)
plt.close(fig)


# Burn-in/stationarity check, normalized by C(0).
fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.3), sharey=True)
for ax, (case, title) in zip(axes, CASES):
    groups = read_rows(ANALYSIS / f"{case}_autocorrelations.csv")
    for color, obs in zip(COLORS, OBS):
        rows = groups[obs]
        t = np.array([r["time"] for r in rows])
        c = np.array([r["C"] for r in rows])
        ct = np.array([r["C_discard_first_quarter"] for r in rows])
        use = t <= 5
        ax.plot(t[use], (ct[use] - c[use]) / abs(c[0]), color=color, lw=1.1, label=LABELS[obs])
    ax.axhline(0, color="0.5", lw=0.8)
    ax.axhspan(-0.05, 0.05, color="0.8", alpha=0.25)
    ax.set_title(title)
    ax.set_xlabel(r"lag $t$")
    ax.grid(True, alpha=0.2)
axes[0].set_ylabel(r"$[C_g^{(3/4)}(t)-C_g^{(\rm full)}(t)]/|C_g(0)|$")
axes[0].legend(ncol=2, frameon=False, loc="lower left")
fig.subplots_adjust(wspace=0.15)
for suffix in ("pdf", "png"):
    fig.savefig(FIGURES / f"stationarity_discard_check.{suffix}", dpi=500)
plt.close(fig)
