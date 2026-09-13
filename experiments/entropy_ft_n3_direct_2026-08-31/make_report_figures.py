#!/usr/bin/env python3
"""Generate fixed descriptive figures for the n=3 direct-sampling report."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw" / "n3_blocks.csv"
ANALYSIS = ROOT / "analysis"
OUT = ROOT / "report" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 160,
    "savefig.dpi": 300,
})

entropy, residual = np.loadtxt(
    RAW, delimiter=",", skiprows=1, usecols=(5, 8), unpack=True)
counts = np.genfromtxt(
    ANALYSIS / "negative_tail_counts.csv", delimiter=",", names=True)
ift = np.genfromtxt(
    ANALYSIS / "medium_entropy_ift.csv", delimiter=",", names=True)

width = float(counts["fd_bin_width"][0])
lower = np.floor(entropy.min() / width) * width
upper = np.ceil(entropy.max() / width) * width
edges = np.arange(lower, upper + 1.01 * width, width)

fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.8))
axes[0].hist(entropy, bins=edges, histtype="step", color="#2166ac",
             linewidth=1.2)
axes[0].axvline(0.0, color="black", linewidth=0.8, linestyle="--")
axes[0].set_yscale("log")
axes[0].set_xlabel(r"medium entropy $\Sigma^{\rm m}_{20}$")
axes[0].set_ylabel("raw block count")
axes[0].text(0.03, 0.94, "41 negative blocks / 1,000,064",
             transform=axes[0].transAxes, va="top")

times = counts["time"]
negative = counts["negative_count"]
axes[1].plot(times, negative, "o-", color="#b2182b", linewidth=1.2)
axes[1].set_yscale("symlog", linthresh=1.0)
axes[1].set_xlabel(r"aggregation time $t$")
axes[1].set_ylabel(r"raw count with $\Sigma^{\rm m}_t<0$")
axes[1].set_xticks(times)
axes[1].tick_params(axis="x", rotation=45)
axes[1].grid(axis="y", alpha=0.2)
fig.tight_layout()
for suffix in ("png", "pdf"):
    fig.savefig(OUT / f"negative_tail_support.{suffix}", bbox_inches="tight")
plt.close(fig)
fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.8))
axes[0].hist(residual, bins=180, histtype="stepfilled", alpha=0.35,
             color="#4d9221", edgecolor="#276419")
axes[0].axvline(0.0, color="black", linewidth=0.8, linestyle="--")
axes[0].set_xlabel(r"$Q_L+Q_R-\Delta E$ at $t=20$")
axes[0].set_ylabel("raw block count")
axes[0].ticklabel_format(axis="x", style="sci", scilimits=(0, 0))

axes[1].plot(ift["time"], ift["exponential_weight_ess"], "o-",
             color="#762a83", label="weight ESS")
axes[1].set_yscale("log")
axes[1].set_xlabel(r"aggregation time $t$")
axes[1].set_ylabel(r"ESS of $\exp(-\Sigma^{\rm m}_t)$")
axes[1].set_xticks(ift["time"])
axes[1].tick_params(axis="x", rotation=45)
axes[1].grid(axis="y", alpha=0.2)
fig.tight_layout()
for suffix in ("png", "pdf"):
    fig.savefig(OUT / f"balance_and_ift_resolution.{suffix}",
                bbox_inches="tight")
plt.close(fig)
