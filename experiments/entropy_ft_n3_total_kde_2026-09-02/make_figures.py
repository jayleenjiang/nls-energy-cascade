#!/usr/bin/env python3
"""Figures for the frozen n=3 KDE total-entropy audit."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
ANALYSIS = ROOT / "analysis"
OUT = ROOT / "report" / "figures"
OUT.mkdir(parents=True, exist_ok=True)


def rows(name: str) -> list[dict[str, str]]:
    with (ANALYSIS / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

tail = rows("equilibrium_kde_tail_accuracy.csv")
regions = ["0-80", "80-95", "95-99", "99-100"]
figure, axes = plt.subplots(1, 2, figsize=(9.2, 3.5))
for dataset, color, marker in (("T6", "#1f77b4", "o"),
                               ("T10", "#d95f02", "s")):
    selected = {row["region"]: row for row in tail if row["dataset"] == dataset}
    values = [float(selected[region]["rmse"]) for region in regions]
    axes[0].plot(regions, values, marker=marker, color=color,
                 linewidth=1.7, label=fr"$T={dataset[1:]}$")
axes[0].axhline(0.10, color="0.35", linestyle="--", linewidth=1,
                label="overall target 0.10")
axes[0].axhline(0.25, color="0.55", linestyle=":", linewidth=1,
                label="tail target 0.25")
axes[0].set_yscale("log")
axes[0].set_ylabel(r"RMSE of $\Delta s_{\rm sys}^{\rm KDE}-\Delta s_{\rm sys}^{\rm exact}$")
axes[0].set_xlabel("maximum-endpoint energy percentile")
axes[0].legend(frameon=False, fontsize=8)
axes[0].grid(alpha=0.2, axis="y")

stability = rows("driven_kde_stability.csv")
all_increment = next(row for row in stability
                     if row["quantity"] == "system_entropy_increment_disagreement"
                     and row["region"] == "all")
tail_increment = next(row for row in stability
                      if row["quantity"] == "system_entropy_increment_disagreement"
                      and row["region"] == "lowest_density_1_percent")
endpoint = next(row for row in stability
                if row["quantity"] == "centered_endpoint_log_density_disagreement")
names = ["endpoint\nlog density", "system\nincrement",
         "lowest-density 1%\nincrement"]
values = [float(endpoint["rmse"]), float(all_increment["rmse"]),
          float(tail_increment["rmse"])]
targets = [0.15, 0.10, 0.25]
x = np.arange(3)
axes[1].bar(x - 0.18, values, width=0.36, color="#8c564b", label="observed RMSE")
axes[1].bar(x + 0.18, targets, width=0.36, color="0.75", label="frozen maximum")
axes[1].set_xticks(x, names)
axes[1].set_yscale("log")
axes[1].set_ylabel("log-density units")
axes[1].legend(frameon=False, fontsize=8)
axes[1].grid(alpha=0.2, axis="y")
figure.tight_layout()
figure.savefig(OUT / "kde_accuracy_audit.pdf", bbox_inches="tight")
figure.savefig(OUT / "kde_accuracy_audit.png", dpi=180, bbox_inches="tight")
plt.close(figure)

bins = rows("symmetric_bin_counts.csv")
shown = [row for row in bins if int(row["bin_index"]) <= 15]
centers = np.array([float(row["a_center"]) for row in shown])
positive = np.array([int(row["positive_count"]) for row in shown])
negative = np.array([int(row["negative_count"]) for row in shown])
figure, axis = plt.subplots(figsize=(5.8, 3.5))
axis.step(centers, positive, where="mid", color="#1f77b4", label=r"$N(+a)$")
axis.step(centers, negative, where="mid", color="#d62728", label=r"$N(-a)$")
axis.axhline(20, color="0.4", linestyle="--", linewidth=1,
             label="frozen per-side minimum")
axis.set_ylim(bottom=0)
axis.set_xlabel(r"medium entropy magnitude $a$")
axis.set_ylabel("raw count")
axis.legend(frameon=False)
axis.grid(alpha=0.2, axis="y")
figure.tight_layout()
figure.savefig(OUT / "medium_entropy_raw_support.pdf", bbox_inches="tight")
figure.savefig(OUT / "medium_entropy_raw_support.png", dpi=180,
               bbox_inches="tight")
plt.close(figure)
