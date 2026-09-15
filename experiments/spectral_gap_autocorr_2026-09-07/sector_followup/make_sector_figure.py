#!/usr/bin/env python3
"""Plot the predeclared early-window damped fits."""

from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
with (RESULTS / "sign_changing_fit_curves.csv").open(newline="") as f:
    curves = list(csv.DictReader(f))
with (RESULTS / "sign_changing_early_damped_fits.csv").open(newline="") as f:
    fits = list(csv.DictReader(f))

case_labels = {
    "driven_T2_T8": r"Driven: $(T_1,T_3)=(2,8)$",
    "equilibrium_T5_T5": r"Equal temperature: $(T_1,T_3)=(5,5)$",
}
obs_labels = {
    "sin_theta1": r"$\sin\theta_1$",
    "I1_minus_I3": r"$I_1-I_3$",
}

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 160,
})

fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.2), sharex=True)
for row, case in enumerate(case_labels):
    for col, obs in enumerate(obs_labels):
        ax = axes[row, col]
        d = [r for r in curves if r["case"] == case and r["observable"] == obs]
        f = next(r for r in fits if r["case"] == case and r["observable"] == obs)
        t = np.asarray([float(r["time"]) for r in d])
        c = np.asarray([float(r["C"]) for r in d])
        se = np.asarray([float(r["SE_C"]) for r in d])
        damped = np.asarray([float(r["damped_fit"]) for r in d])
        null = np.asarray([float(r["real_exponential_fit"]) for r in d])
        ax.fill_between(t, c - 1.96 * se, c + 1.96 * se,
                        color="#4C78A8", alpha=0.18, linewidth=0)
        ax.plot(t, c, color="#1f4e79", lw=1.5, label="saved mean correlation")
        ax.plot(t, damped, color="#D55E00", lw=1.6,
                label="damped fit")
        ax.plot(t, null, color="0.35", lw=1.1, ls="--",
                label="real-exponential null")
        ax.axhline(0.0, color="0.5", lw=0.7)
        ax.axvline(float(f["mean_correlation_first_zero"]), color="#009E73", lw=0.9, ls=":")
        ax.set_title(f"{case_labels[case]} — {obs_labels[obs]}")
        ax.set_ylabel(r"$C_g(t)$")
        ax.grid(alpha=0.2)
        ax.text(0.98, 0.95,
                rf"$\lambda_R={float(f['lambda_R']):.2f}$" + "\n" +
                rf"$\lambda_I={float(f['lambda_I']):.2f}$" + "\n" +
                rf"$\Delta\mathrm{{AIC}}={float(f['delta_AIC_damped_minus_exponential']):.0f}$",
                transform=ax.transAxes, ha="right", va="top", fontsize=8)
for ax in axes[-1]:
    ax.set_xlabel(r"lag $t$")
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False,
           bbox_to_anchor=(0.5, -0.01))
fig.tight_layout(rect=(0, 0.06, 1, 1))
fig.savefig(ROOT / "early_damped_sector_fits.pdf", bbox_inches="tight")
fig.savefig(ROOT / "early_damped_sector_fits.png", bbox_inches="tight", dpi=240)
