#!/usr/bin/env python3
"""Publication figures for the complete stationary action-profile dataset."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PROFILE_DIR = ROOT / "data" / "profiles"
FIGURE_DIR = ROOT / "Figures"
SUMMARY_FILE = ROOT / "data" / "logical_run_summary.csv"

BCS = ("BC1", "BC2", "BC3", "BC3b")
TEMPS = ("T10_T2", "T4_T6", "T6_T6")
NS = (25, 50, 100)

# Okabe-Ito colors; markers and line styles keep the plots legible in grayscale.
N_STYLE = {
    25: ("#0072B2", "o", "-"),
    50: ("#D55E00", "s", "--"),
    100: ("#009E73", "^", "-."),
}
BC_STYLE = {
    "BC1": ("#0072B2", "o", "-"),
    "BC2": ("#E69F00", "s", "--"),
    "BC3": ("#009E73", "^", "-"),
    "BC3b": ("#CC79A7", "D", "--"),
}
TEMP_TITLES = {
    "T10_T2": r"$(T_1,T_n)=(10,2)$",
    "T4_T6": r"$(T_1,T_n)=(4,6)$",
    "T6_T6": r"$(T_1,T_n)=(6,6)$",
}


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["STIX Two Text", "STIXGeneral", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "font.size": 8.5,
            "axes.labelsize": 9.5,
            "axes.titlesize": 10,
            "legend.fontsize": 8.5,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.45,
            "savefig.transparent": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def profile_path(bc: str, temp: str, n: int) -> Path:
    return PROFILE_DIR / f"{bc}_{temp}_n{n}_profile.csv"


def load_profile(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(line for line in handle if not line.startswith("#")))
    return {
        "j": np.array([int(row["j"]) for row in rows], dtype=int),
        "mean_I": np.array([float(row["mean_I"]) for row in rows]),
        "se_mean_I": np.array([float(row["se_mean_I"]) for row in rows]),
    }


def load_summary() -> list[dict[str, str]]:
    with SUMMARY_FILE.open(newline="") as handle:
        return list(csv.DictReader(handle))


def panel_label(index: int) -> str:
    return f"({chr(ord('a') + index)})"


def save_figure(fig: mpl.figure.Figure, stem: str) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_DIR / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIGURE_DIR / f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def plot_profile_atlas() -> None:
    fig, axes = plt.subplots(
        len(BCS),
        len(TEMPS),
        figsize=(7.15, 8.35),
        sharex=True,
        sharey="row",
        constrained_layout=False,
    )

    for row, bc in enumerate(BCS):
        for col, temp in enumerate(TEMPS):
            ax = axes[row, col]
            for n in NS:
                data = load_profile(profile_path(bc, temp, n))
                x = (data["j"] - 1) / (n - 1)
                mean = data["mean_I"]
                ci = 1.96 * data["se_mean_I"]
                color, marker, linestyle = N_STYLE[n]
                markevery = max(1, n // 10)
                ax.fill_between(x, mean - ci, mean + ci, color=color, alpha=0.12, linewidth=0)
                ax.plot(
                    x,
                    mean,
                    color=color,
                    linestyle=linestyle,
                    marker=marker,
                    markersize=2.8,
                    markevery=markevery,
                    markerfacecolor="white",
                    markeredgewidth=0.65,
                    label=fr"$n={n}$",
                )

            if row == 0:
                ax.set_title(TEMP_TITLES[temp], pad=7)
            if col == 0:
                ax.set_ylabel(fr"{bc}: $\langle I_j\rangle$")
            if row == len(BCS) - 1:
                ax.set_xlabel(r"normalized site $\xi=(j-1)/(n-1)$")
            ax.text(
                0.025,
                0.93,
                panel_label(row * len(TEMPS) + col),
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontweight="normal",
            )
            ax.set_xlim(0, 1)
            ax.grid(axis="y", color="0.88", linewidth=0.5)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.54, 0.995), ncol=3, frameon=False)
    fig.subplots_adjust(left=0.12, right=0.99, bottom=0.065, top=0.945, hspace=0.25, wspace=0.22)
    save_figure(fig, "action_profile_atlas")


def plot_midpoint_scaling() -> None:
    rows = load_summary()
    fig, axes = plt.subplots(1, 3, figsize=(7.15, 2.85), sharex=True, sharey=True)

    inventory_rows: list[dict[str, str | int | float]] = []
    for col, temp in enumerate(TEMPS):
        ax = axes[col]
        for bc in BCS:
            subset = sorted(
                (row for row in rows if row["temp_label"] == temp and row["bc_label"] == bc),
                key=lambda row: int(row["n"]),
            )
            n = np.array([int(row["n"]) for row in subset], dtype=float)
            y = np.array([float(row["mean_I_mid"]) for row in subset])
            se = np.array([float(row["se_I_mid_conservative"]) for row in subset])
            y0, se0 = y[0], se[0]
            yn = y / y0
            # Conservative propagation treats the n=25 normalizer as independent.
            sen = yn * np.sqrt((se / y) ** 2 + (se0 / y0) ** 2)
            slope, intercept = np.polyfit(np.log(n), np.log(y), 1)
            fit_n = np.geomspace(NS[0], NS[-1], 120)
            fit_norm = np.exp(intercept) * fit_n**slope / y0
            color, marker, linestyle = BC_STYLE[bc]
            ax.errorbar(
                n,
                yn,
                yerr=1.96 * sen,
                color=color,
                marker=marker,
                linestyle="none",
                markersize=4.0,
                markerfacecolor="white",
                markeredgewidth=0.9,
                capsize=2.0,
                label=bc,
            )
            ax.plot(fit_n, fit_norm, color=color, linestyle=linestyle, linewidth=1.35)
            inventory_rows.append(
                {
                    "bc_label": bc,
                    "temp_label": temp,
                    "n_values": "25;50;100",
                    "midpoint_values": ";".join(f"{value:.16g}" for value in y),
                    "normalized_midpoint_values": ";".join(f"{value:.16g}" for value in yn),
                    "loglog_slope": float(slope),
                }
            )

        ref_n = np.geomspace(NS[0], NS[-1], 120)
        ax.plot(ref_n, (ref_n / NS[0]) ** (-0.5), color="0.25", linewidth=1.0, linestyle=":",
                label=r"$n^{-1/2}$" if col == 0 else None)
        ax.plot(ref_n, np.ones_like(ref_n), color="0.55", linewidth=1.0, linestyle=(0, (5, 2)),
                label=r"$n^0$" if col == 0 else None)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.xaxis.set_major_locator(mticker.FixedLocator(NS))
        ax.xaxis.set_major_formatter(mticker.FixedFormatter([str(value) for value in NS]))
        ax.xaxis.set_minor_locator(mticker.NullLocator())
        ax.xaxis.set_minor_formatter(mticker.NullFormatter())
        ax.set_xlim(22.5, 112)
        ax.set_ylim(0.47, 1.23)
        ax.set_title(TEMP_TITLES[temp])
        ax.set_xlabel(r"chain length $n$")
        ax.text(0.04, 0.94, panel_label(col), transform=ax.transAxes, ha="left", va="top")
        ax.grid(which="major", axis="both", color="0.88", linewidth=0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[0].set_ylabel(r"normalized midpoint action $\langle I\rangle_{\rm mid}/\langle I\rangle_{\rm mid,25}$")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.52, 1.025), ncol=6, frameon=False)
    fig.subplots_adjust(left=0.105, right=0.99, bottom=0.19, top=0.78, wspace=0.17)
    save_figure(fig, "action_midpoint_scaling")

    with (ROOT / "data" / "figure_midpoint_source.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(inventory_rows[0]))
        writer.writeheader()
        writer.writerows(inventory_rows)


def write_figure_inventory() -> None:
    rows = []
    for bc in BCS:
        for temp in TEMPS:
            for n in NS:
                path = profile_path(bc, temp, n)
                data = load_profile(path)
                rows.append(
                    {
                        "bc_label": bc,
                        "temp_label": temp,
                        "n": n,
                        "profile_rows": len(data["j"]),
                        "source_file": str(path.relative_to(ROOT)),
                    }
                )
    with (ROOT / "data" / "figure_profile_inventory.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    configure_style()
    write_figure_inventory()
    plot_profile_atlas()
    plot_midpoint_scaling()


if __name__ == "__main__":
    main()
