#!/usr/bin/env python3
"""Create fixed-format diagnostic figures from merged profile CSV files."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
PROFILES = ROOT / "curated_results" / "profiles"
FIGURES = ROOT / "figures"
BCS = ["BC1", "BC2", "BC3", "BC3b"]
TEMPS = ["T10_T2", "T4_T6", "T6_T6"]
NS = [25, 50, 100]
COLORS = {25: "#0072B2", 50: "#D55E00", 100: "#009E73"}


def load(path):
    with path.open() as handle:
        rows = list(csv.DictReader(line for line in handle if not line.startswith("#")))
    return rows


def save(fig, stem):
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / f"{stem}.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def profile_figures():
    for temp in TEMPS:
        for quantity in ("action", "sine"):
            fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), sharex=True)
            for ax, bc in zip(axes.flat, BCS):
                for n in NS:
                    rows = load(PROFILES / f"{bc}_{temp}_n{n}_profile.csv")
                    if quantity == "action":
                        y = np.array([float(row["mean_I"]) for row in rows])
                        se = np.array([float(row["se_mean_I"]) for row in rows])
                        x = np.arange(1, n + 1) / n
                    else:
                        y = np.array([float(row["mean_sin_theta"]) for row in rows[:-1]])
                        se = np.array([float(row["se_mean_sin_theta"]) for row in rows[:-1]])
                        x = (np.arange(1, n) + 0.5) / n
                    ax.plot(x, y, lw=1.25, color=COLORS[n], label=f"n={n}")
                    ax.fill_between(x, y - 1.96 * se, y + 1.96 * se,
                                    color=COLORS[n], alpha=0.13, linewidth=0)
                ax.axhline(0, color="0.65", lw=0.7)
                ax.set_title(bc)
                ax.grid(alpha=0.18)
            axes[1, 0].set_xlabel("normalized position")
            axes[1, 1].set_xlabel("normalized position")
            axes[0, 0].set_ylabel(r"$\langle I_j\rangle$" if quantity == "action"
                                  else r"$\langle\sin\theta_j\rangle$")
            axes[1, 0].set_ylabel(r"$\langle I_j\rangle$" if quantity == "action"
                                  else r"$\langle\sin\theta_j\rangle$")
            axes[0, 0].legend(frameon=False, ncol=3, fontsize=8)
            fig.tight_layout()
            save(fig, f"{quantity}_profiles_{temp}")


def midpoint_figure():
    with (ROOT / "analysis" / "logical_run_summary.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.8), sharey=False)
    markers = {"BC1": "o", "BC2": "s", "BC3": "^", "BC3b": "D"}
    for ax, temp in zip(axes, TEMPS):
        for bc in BCS:
            subset = sorted((r for r in rows if r["temp_label"] == temp and r["bc_label"] == bc),
                            key=lambda r: int(r["n"]))
            n = np.array([int(r["n"]) for r in subset])
            y = np.array([float(r["mean_I_mid"]) for r in subset])
            se = np.array([float(r["se_I_mid_conservative"]) for r in subset])
            ax.errorbar(n, y, yerr=1.96 * se, marker=markers[bc], lw=1.2,
                        capsize=2, label=bc)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks(NS, labels=[str(n) for n in NS])
        ax.set_title(temp.replace("_", ", "))
        ax.set_xlabel("chain length n")
        ax.grid(alpha=0.2, which="both")
    axes[0].set_ylabel(r"central $\langle I\rangle$")
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    save(fig, "midpoint_scaling")


if __name__ == "__main__":
    profile_figures()
    midpoint_figure()
