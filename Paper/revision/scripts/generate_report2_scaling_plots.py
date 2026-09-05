#!/usr/bin/env python3
"""Generate separate scaling-line plots for report (2).tex."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path("/Users/jayleenjiang/Documents/NLS")
OUT = ROOT / "Paper/revision_2026-06-19/report_assets"
SIMD_DIR = ROOT / "flux/flux_data/fixed_simd_2026-06-19"
CANON_DIR = (
    ROOT
    / "Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4"
)
NS = np.array([10, 20, 30, 40], dtype=float)


plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "legend.fontsize": 8.5,
        "legend.frameon": False,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.22,
        "grid.linestyle": ":",
        "lines.linewidth": 1.8,
        "lines.markersize": 5,
    }
)


def fit_power_law(means: np.ndarray) -> tuple[float, float, float]:
    log_n = np.log(NS)
    log_y = np.log(means)
    slope, intercept = np.polyfit(log_n, log_y, 1)
    pred = slope * log_n + intercept
    r2 = 1.0 - np.sum((log_y - pred) ** 2) / np.sum((log_y - log_y.mean()) ** 2)
    return math.exp(intercept), -slope, float(r2)


def load_simd() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    means, stds, ses = [], [], []
    for n in NS.astype(int):
        samples = np.loadtxt(SIMD_DIR / f"flux_n{n}.txt", comments="#")
        means.append(float(samples.mean()))
        stds.append(float(samples.std(ddof=1)))
        ses.append(float(samples.std(ddof=1) / math.sqrt(samples.size)))
    return np.array(means), np.array(stds), np.array(ses)


def load_canonical() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    means, stds, ses, burnins = [], [], [], []
    for n in NS.astype(int):
        with (CANON_DIR / f"n{n}_summary.csv").open(newline="") as f:
            row = next(csv.DictReader(f))
        means.append(float(row["mean_action_current"]))
        stds.append(float(row["sample_sd"]))
        ses.append(float(row["standard_error"]))
        burnins.append(float(row["burnin"]))
    return np.array(means), np.array(stds), np.array(ses), np.array(burnins)


def make_scaling_plot(
    means: np.ndarray,
    ses: np.ndarray,
    out_stem: str,
    title: str,
    color: str,
) -> dict[str, float]:
    coefficient, alpha, r2 = fit_power_law(means)
    grid = np.linspace(NS.min(), NS.max(), 200)
    fit = coefficient * grid ** (-alpha)

    fig, ax = plt.subplots(figsize=(4.6, 3.15))
    ax.errorbar(
        NS,
        means,
        yerr=ses,
        fmt="o",
        color=color,
        ecolor=color,
        elinewidth=1.0,
        capsize=3,
        label="simulation mean",
        zorder=3,
    )
    ax.plot(
        grid,
        fit,
        "--",
        color="#333333",
        label=rf"${coefficient:.2f}\,n^{{-{alpha:.3f}}}$",
        zorder=2,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("chain length $n$")
    ax.set_ylabel("mean current $\\langle J\\rangle$")
    ax.set_title(title)
    ax.legend(loc="upper right")
    ax.text(
        0.04,
        0.06,
        rf"$R^2={r2:.4f}$",
        transform=ax.transAxes,
        fontsize=9,
        color="#333333",
    )
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"{out_stem}.{suffix}")
    plt.close(fig)
    return {"coefficient": coefficient, "alpha": alpha, "r2": r2}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    simd_means, simd_stds, simd_ses = load_simd()
    canon_means, canon_stds, canon_ses, burnins = load_canonical()

    metadata = {
        "simd": {
            "means": simd_means.tolist(),
            "stds": simd_stds.tolist(),
            "ses": simd_ses.tolist(),
            **make_scaling_plot(
                simd_means,
                simd_ses,
                "simd_fixed_scaling_line",
                "SIMD fixed-current scaling",
                "#0072B2",
            ),
        },
        "canonical": {
            "means": canon_means.tolist(),
            "stds": canon_stds.tolist(),
            "ses": canon_ses.tolist(),
            "burnins": burnins.tolist(),
            **make_scaling_plot(
                canon_means,
                canon_ses,
                "canonical_scaling_line",
                "Canonical current scaling",
                "#D55E00",
            ),
        },
    }
    (OUT / "report2_scaling_plot_metrics.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
