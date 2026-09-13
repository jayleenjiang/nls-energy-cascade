#!/usr/bin/env python3
"""Generate the frozen equilibrium-audit report figures."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
ANALYSIS = ROOT / "analysis"
OUT = ROOT / "report" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 11,
    "legend.fontsize": 9,
    "figure.dpi": 180,
})


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


means = read_rows(ANALYSIS / "means_raw.csv")
lookup = {(row["label"], row["observable"]): row for row in means}
colors = {"T6": "#2166ac", "T10": "#b2182b"}

fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6))
observables = ["q_left_rate", "q_right_rate", "energy_current"]
labels = [r"$Q_L/t$", r"$Q_R/t$", r"$J_E$"]
x = np.arange(len(observables), dtype=float)
for offset, case in [(-0.10, "T6"), (0.10, "T10")]:
    rows = [lookup[(case, observable)] for observable in observables]
    values = np.array([float(row["estimate"]) for row in rows])
    lows = np.array([float(row["bootstrap_ci_low"]) for row in rows])
    highs = np.array([float(row["bootstrap_ci_high"]) for row in rows])
    axes[0].errorbar(
        x + offset, values, yerr=np.vstack([values - lows, highs - values]),
        fmt="o", capsize=4, color=colors[case], label=case.replace("T", r"$T=") + "$")
axes[0].axhline(0.0, color="black", linewidth=0.9)
axes[0].set_xticks(x, labels)
axes[0].set_ylabel("mean rate / current")
axes[0].legend(frameon=False)
axes[0].grid(axis="y", alpha=0.25)

for offset, case in [(-0.06, "T6"), (0.06, "T10")]:
    row = lookup[(case, "entropy_rate")]
    value = float(row["estimate"])
    low = float(row["bootstrap_ci_low"])
    high = float(row["bootstrap_ci_high"])
    axes[1].errorbar(
        [offset], [value], yerr=[[value - low], [high - value]], fmt="o",
        capsize=5, color=colors[case], label=case.replace("T", r"$T=") + "$")
    axes[1].annotate(
        rf"${float(row['sigma_from_zero']):+.2f}\,\mathrm{{SE}}$",
        (offset, value), xytext=(5, 8), textcoords="offset points",
        color=colors[case])
axes[1].axhline(0.0, color="black", linewidth=0.9)
axes[1].set_xlim(-0.25, 0.25)
axes[1].set_xticks([])
axes[1].set_ylabel(r"$\langle\Sigma^m/t\rangle$")
axes[1].ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
axes[1].legend(frameon=False)
axes[1].grid(axis="y", alpha=0.25)
fig.tight_layout()
for suffix in ("png", "pdf"):
    fig.savefig(OUT / f"equilibrium_means.{suffix}", bbox_inches="tight")
plt.close(fig)

bins = read_rows(ANALYSIS / "symmetry_bins.csv")
fits = {row["label"]: row for row in read_rows(ANALYSIS / "symmetry_fit.csv")}
fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6), sharey=True)
for axis, case in zip(axes, ("T6", "T10")):
    rows = [row for row in bins if row["label"] == case and row["supported"] == "1"]
    xx = np.array([float(row["abs_rate_center"]) for row in rows])
    yy = np.array([float(row["log_ratio_over_t"]) for row in rows])
    fit = fits[case]
    slope = float(fit["slope"])
    intercept = float(fit["intercept"])
    axis.scatter(xx, yy, s=13, alpha=0.72, color=colors[case], edgecolor="none")
    grid = np.linspace(0.0, xx.max(), 200)
    axis.plot(grid, intercept + slope * grid, color="black", linewidth=1.3,
              label="frozen WLS fit")
    axis.axhline(0.0, color="0.45", linestyle="--", linewidth=0.9,
                 label="equilibrium reference")
    axis.set_xlabel(r"$a=|\Sigma^m|/t$")
    axis.set_title(case.replace("T", r"$T=") + "$")
    axis.grid(alpha=0.2)
    axis.text(
        0.04, 0.94,
        rf"slope $={slope:+.3g}$" + "\n" +
        rf"95\% CI $=[{float(fit['slope_ci_low']):+.3g},"
        rf"{float(fit['slope_ci_high']):+.3g}]$",
        transform=axis.transAxes, va="top")
axes[0].set_ylabel(r"$t^{-1}\log[p(a)/p(-a)]$")
axes[0].legend(frameon=False, loc="lower left")
fig.tight_layout()
for suffix in ("png", "pdf"):
    fig.savefig(OUT / f"equilibrium_symmetry.{suffix}", bbox_inches="tight")
plt.close(fig)
