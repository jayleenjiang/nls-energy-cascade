#!/usr/bin/env python3
"""Publication-style figures for the controlled relaxation experiment."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
ANALYSIS = ROOT / "analysis"
FIGURES = ROOT / "figures"
CASES = [
    ("driven_dt1e-3", r"Driven, $\Delta t=10^{-3}$"),
    ("driven_dt2p5e-4", r"Driven, $\Delta t=2.5\times10^{-4}$"),
    ("equilibrium_dt1e-3", r"Equal $T$, $\Delta t=10^{-3}$"),
    ("equilibrium_dt2p5e-4", r"Equal $T$, $\Delta t=2.5\times10^{-4}$"),
]
RATE_OBS = [
    ("cos_theta3", r"$\cos\theta_3$"),
    ("I2", r"$I_2$"),
    ("I1_plus_I3", r"$I_1+I_3$"),
    ("cos_even_sum", r"$\cos\theta_1+\cos\theta_3$"),
]
ODD_OBS = [
    ("I1_minus_I3", r"$I_1-I_3$"),
    ("cos_odd_diff", r"$\cos\theta_1-\cos\theta_3$"),
    ("sin_odd_diff", r"$\sin\theta_1-\sin\theta_3$"),
]


def rows(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def main():
    FIGURES.mkdir(exist_ok=True)
    plt.rcParams.update({
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "legend.fontsize": 7.5,
        "figure.dpi": 160,
        "savefig.dpi": 300,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.5), sharex=True, sharey=True)
    for ax, (case, title) in zip(axes.flat, CASES):
        data = rows(ANALYSIS / f"{case}_autocorrelations.csv")
        plateau = {
            x["observable"]: x
            for x in rows(ANALYSIS / f"{case}_plateaus.csv")
        }
        for obs, label in RATE_OBS:
            r = [x for x in data if x["observable"] == obs]
            t = np.asarray([float(x["time"]) for x in r])
            y = np.asarray([float(x["lambda_eff"]) for x in r])
            lo = np.asarray([float(x["lambda_ci_low"]) for x in r])
            hi = np.asarray([float(x["lambda_ci_high"]) for x in r])
            mask = (t >= 0.8) & (t <= 4.5) & np.isfinite(y)
            line, = ax.plot(t[mask], y[mask], lw=1.20, label=label)
            ax.fill_between(t[mask], lo[mask], hi[mask], alpha=0.10)
            accepted = plateau[obs]
            if accepted["plateau_start"]:
                ax.plot(
                    [float(accepted["plateau_start"]), float(accepted["plateau_end"])],
                    [float(accepted["lambda_R"]), float(accepted["lambda_R"])],
                    color=line.get_color(), lw=3.2, solid_capstyle="butt",
                )
        ax.axvline(1.0, color="0.55", lw=0.8, ls=":")
        ax.axhline(0.0, color="0.7", lw=0.7)
        ax.set_title(title)
        ax.set_xlim(0.8, 3.6)
        ax.set_ylim(-2.0, -0.25)
        ax.grid(alpha=0.18)
    axes[1, 0].set_xlabel("lag $t$")
    axes[1, 1].set_xlabel("lag $t$")
    axes[0, 0].set_ylabel(r"$\lambda_{\rm eff}(t)$")
    axes[1, 0].set_ylabel(r"$\lambda_{\rm eff}(t)$")
    axes[0, 0].legend(ncol=2, frameon=False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGURES / f"controlled_local_rates.{ext}", bbox_inches="tight")
    plt.close(fig)

    curve = rows(ANALYSIS / "early_odd_fit_curves.csv")
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.5), sharex=True)
    for ax, (case, title) in zip(axes.flat, CASES):
        for obs, label in ODD_OBS:
            r = [x for x in curve if x["case"] == case and x["observable"] == obs]
            t = np.asarray([float(x["time"]) for x in r])
            y = np.asarray([float(x["C"]) for x in r])
            se = np.asarray([float(x["SE_C"]) for x in r])
            fit = np.asarray([float(x["damped_fit"]) for x in r])
            line, = ax.plot(t, fit, lw=1.4, label=label)
            ax.fill_between(t, y - se, y + se, color=line.get_color(), alpha=0.12)
            ax.plot(t[::5], y[::5], ".", ms=2.2, color=line.get_color())
        ax.axhline(0.0, color="0.55", lw=0.8)
        ax.set_title(title)
        ax.grid(alpha=0.18)
    axes[1, 0].set_xlabel("lag $t$")
    axes[1, 1].set_xlabel("lag $t$")
    axes[0, 0].set_ylabel("autocovariance")
    axes[1, 0].set_ylabel("autocovariance")
    axes[0, 0].legend(frameon=False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGURES / f"controlled_odd_damped_fits.{ext}", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
