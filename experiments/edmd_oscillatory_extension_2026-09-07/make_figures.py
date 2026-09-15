#!/usr/bin/env python3
"""Publication figures for the targeted EDMD completeness test."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
AN = HERE / "analysis"
OUT = HERE / "figures"
OUT.mkdir(exist_ok=True)
CASES = ("driven_dt1e-3", "driven_dt2p5e-4",
         "equilibrium_dt1e-3", "equilibrium_dt2p5e-4")
LABEL = {"driven_dt1e-3": r"driven, $\Delta t=10^{-3}$",
         "driven_dt2p5e-4": r"driven, $\Delta t=2.5\times10^{-4}$",
         "equilibrium_dt1e-3": r"equilibrium, $\Delta t=10^{-3}$",
         "equilibrium_dt2p5e-4": r"equilibrium, $\Delta t=2.5\times10^{-4}$"}
AUTO = dict(zip(CASES, (5.4009123, 5.4996905, 5.0809736, 5.2636681)))


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 160,
                     "savefig.dpi": 300})

# Short-lag convergence.
fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.2), sharex=True)
for case, ax in zip(CASES, axes.flat):
    rr = [r for r in rows(AN / "oscillatory_convergence.csv")
          if r["case"] == case and r["dimension"] == "short_lag"]
    t = np.array([float(r["setting"]) for r in rr])
    re = np.array([float(r["lambda_real"]) for r in rr])
    im = np.array([float(r["lambda_imag"]) for r in rr])
    ax.plot(t, im, "o-", color="#1565c0", label=r"EDMD $\Im\lambda$")
    ax.axhline(AUTO[case], color="#d32f2f", ls="--", lw=1.2,
               label="autocorrelation")
    ax2 = ax.twinx()
    ax2.plot(t, re, "s--", color="#616161", ms=4, label=r"$\Re\lambda$")
    ax2.set_ylabel(r"$\Re\lambda$", color="#616161")
    ax.set_title(LABEL[case])
    ax.set_ylabel(r"$\Im\lambda$")
    ax.grid(alpha=.2)
axes[-1, 0].set_xlabel(r"EDMD lag $\tau$")
axes[-1, 1].set_xlabel(r"EDMD lag $\tau$")
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False)
fig.tight_layout(rect=(0, 0, 1, .94))
fig.savefig(OUT / "oscillatory_short_lag_convergence.pdf", bbox_inches="tight")
fig.savefig(OUT / "oscillatory_short_lag_convergence.png", bbox_inches="tight")
plt.close(fig)

# Dictionary convergence at tau=.05.
fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8))
for case in CASES:
    rr = [r for r in rows(AN / "oscillatory_convergence.csv")
          if r["case"] == case and r["dimension"] == "dictionary"]
    x = np.arange(3)
    axes[0].plot(x, [float(r["lambda_real"]) for r in rr], "o-", label=LABEL[case])
    axes[1].plot(x, [float(r["lambda_imag"]) for r in rr], "o-", label=LABEL[case])
for ax, ylabel in zip(axes, (r"$\Re\lambda$", r"$\Im\lambda$")):
    ax.set_xticks(range(3), ("E1", "E2", "E3"))
    ax.set_ylabel(ylabel)
    ax.set_xlabel("targeted dictionary")
    ax.grid(alpha=.2)
axes[1].legend(fontsize=7, frameon=False, loc="best")
fig.tight_layout()
fig.savefig(OUT / "oscillatory_dictionary_convergence.pdf", bbox_inches="tight")
fig.savefig(OUT / "oscillatory_dictionary_convergence.png", bbox_inches="tight")
plt.close(fig)

# Old/new leading real rates with CIs.
old_point = np.array((-0.9401736, -0.9631910, -0.9551489, -0.9407584))
old_lo = np.array((-0.954628, -0.980007, -0.972388, -0.954830))
old_hi = np.array((-0.920824, -0.941942, -0.932803, -0.922965))
rr = rows(AN / "leading_real_summary.csv")
new = np.array([float(r["bootstrap_mean_real"]) for r in rr])
new_lo = np.array([float(r["bootstrap_ci_low"]) for r in rr])
new_hi = np.array([float(r["bootstrap_ci_high"]) for r in rr])
y = np.arange(4)
fig, ax = plt.subplots(figsize=(7.2, 3.0))
ax.errorbar(old_point, y + .10, xerr=np.vstack((old_point-old_lo, old_hi-old_point)),
            fmt="o", color="#757575", capsize=2, label="original dictionary")
ax.errorbar(new, y - .10, xerr=np.vstack((new-new_lo, new_hi-new)),
            fmt="o", color="#1565c0", capsize=2, label="targeted extension")
ax.set_yticks(y, [LABEL[c] for c in CASES])
ax.set_xlabel(r"leading visible real rate $\Re\lambda$")
ax.grid(axis="x", alpha=.2)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT / "leading_real_old_vs_extended.pdf", bbox_inches="tight")
fig.savefig(OUT / "leading_real_old_vs_extended.png", bbox_inches="tight")
plt.close(fig)

# Modal weights.
obs = ("I2", "I1_plus_I3", "cos_even_sum", "cos_theta3")
pretty = (r"$I_2$", r"$I_1+I_3$", r"$\cos\theta_1+\cos\theta_3$", r"$\cos\theta_3$")
wrows = rows(AN / "leading_modal_weights.csv")
W = np.array([[float(next(r["weight_abs"] for r in wrows
                          if r["case"] == c and r["observable"] == o))
               for o in obs] for c in CASES])
fig, ax = plt.subplots(figsize=(6.8, 3.2))
im = ax.imshow(W, cmap="Blues", vmin=0, vmax=.55, aspect="auto")
ax.set_xticks(range(4), pretty)
ax.set_yticks(range(4), [LABEL[c] for c in CASES])
for i in range(4):
    for j in range(4):
        ax.text(j, i, f"{W[i,j]:.3f}", ha="center", va="center",
                color="white" if W[i,j] > .32 else "black")
fig.colorbar(im, ax=ax, label="absolute modal weight")
fig.tight_layout()
fig.savefig(OUT / "leading_modal_weights.pdf", bbox_inches="tight")
fig.savefig(OUT / "leading_modal_weights.png", bbox_inches="tight")
plt.close(fig)

# Eigenfunction slices: one 2x3 panel per case, real part normalized per mode.
for case in CASES:
    data = rows(AN / f"{case}_eigenfunction_slices.csv")
    fig, axes = plt.subplots(2, 3, figsize=(8.0, 5.0), sharex=True, sharey=True)
    for i, mode in enumerate(("oscillatory", "leading_real")):
        mode_rows = [r for r in data if r["mode"] == mode]
        maxabs = max(abs(float(r["phi_real"])) for r in mode_rows) or 1.0
        for j, action in enumerate((1.0, 2.0, 4.0)):
            q = [r for r in mode_rows if abs(float(r["I_equal"])-action) < 1e-9]
            th1 = np.unique([float(r["theta1"]) for r in q])
            th3 = np.unique([float(r["theta3"]) for r in q])
            Z = np.array([float(r["phi_real"]) / maxabs for r in q]).reshape(len(th1), len(th3))
            m = axes[i, j].pcolormesh(th3, th1, Z, shading="auto", cmap="RdBu_r", vmin=-1, vmax=1)
            axes[i, j].set_title(fr"{mode.replace('_',' ')}, $I={action:g}$", fontsize=8)
            if i == 1:
                axes[i, j].set_xlabel(r"$\theta_3$")
            if j == 0:
                axes[i, j].set_ylabel(r"$\theta_1$")
    cax = fig.add_axes((.91, .16, .018, .68))
    fig.colorbar(m, cax=cax, label="normalized real part")
    fig.suptitle(LABEL[case], y=.99)
    fig.subplots_adjust(left=.08, right=.88, bottom=.09, top=.91, wspace=.12, hspace=.28)
    fig.savefig(OUT / f"{case}_eigenfunction_slices.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{case}_eigenfunction_slices.png", bbox_inches="tight")
    plt.close(fig)

print(f"wrote figures to {OUT}")
