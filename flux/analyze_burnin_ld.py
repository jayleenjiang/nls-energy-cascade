#!/usr/bin/env python3
"""Analyze burn-in relaxation and finite-tau action-current tails.

The C++ runner writes two kinds of outputs:

* transient: *_profile_timeseries.csv and *_flux_timeseries.csv
* tau:       *_blocks.csv and *_summary.csv

This script generates the plots requested for the burn-in/large-deviation
diagnostic study.  It is descriptive: the rate-proxy plots are only evidence
for or against tau-collapse, not a proof of a large-deviation principle.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-dir", type=Path, required=True)
    parser.add_argument(
        "--threshold-step",
        type=float,
        default=0.01,
        help="A-grid spacing for P[J_tau > A].",
    )
    parser.add_argument(
        "--threshold-max",
        type=float,
        default=0.70,
        help="Largest A threshold to tabulate.",
    )
    parser.add_argument(
        "--survival-fit-low",
        type=float,
        default=0.005,
        help="Lowest survival probability included in descriptive slope fit.",
    )
    parser.add_argument(
        "--survival-fit-high",
        type=float,
        default=0.20,
        help="Highest survival probability included in descriptive slope fit.",
    )
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def plot_transient_profiles(exp_dir: Path, out_dir: Path) -> None:
    paths = sorted(exp_dir.glob("transient_n*_profile_timeseries.csv"))
    if not paths:
        return

    for quantity, column, ylabel, suffix in [
        (
            "terminal",
            "mean_terminal_action",
            "mean terminal action",
            "terminal_profile",
        ),
        (
            "cumulative",
            "mean_cumulative_action",
            "mean cumulative action over [0,T]",
            "cumulative_profile",
        ),
    ]:
        fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.5), sharey=False)
        axes = axes.ravel()
        for axis, path in zip(axes, paths):
            rows = read_csv(path)
            n = int(rows[0]["n"])
            times = sorted({float(row["time"]) for row in rows})
            cmap = plt.get_cmap("viridis")
            for index, time in enumerate(times):
                subset = [row for row in rows if float(row["time"]) == time]
                modes = np.array([int(row["mode"]) for row in subset])
                values = np.array([float(row[column]) for row in subset])
                x = modes / max(1, n - 1)
                color = cmap(index / max(1, len(times) - 1))
                axis.plot(x, values, lw=1.2, color=color, label=f"T={time:g}")
            axis.set_title(f"n={n}")
            axis.set_xlabel(r"normalized mode $j/(n-1)$")
            axis.set_ylabel(ylabel)
            axis.grid(True, alpha=0.25)
        handles, labels = axes[-1].get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            loc="center right",
            frameon=False,
            fontsize=8,
            ncol=1,
        )
        fig.suptitle(f"No-burn-in {quantity} action profiles")
        fig.tight_layout(rect=(0, 0, 0.88, 0.95))
        fig.savefig(out_dir / f"burnin_{suffix}.png", dpi=200)
        fig.savefig(out_dir / f"burnin_{suffix}.pdf")
        plt.close(fig)


def plot_transient_flux(exp_dir: Path, out_dir: Path) -> None:
    paths = sorted(exp_dir.glob("transient_n*_flux_timeseries.csv"))
    if not paths:
        return

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4), sharex=True)
    for path in paths:
        rows = read_csv(path)
        n = int(rows[0]["n"])
        time = np.array([float(row["time"]) for row in rows])
        cumulative = np.array(
            [float(row["mean_cumulative_current"]) for row in rows]
        )
        cumulative_se = np.array(
            [float(row["se_cumulative_current"]) for row in rows]
        )
        interval = np.array(
            [float(row["mean_last_interval_current"]) for row in rows]
        )
        interval_se = np.array(
            [float(row["se_last_interval_current"]) for row in rows]
        )
        axes[0].errorbar(
            time,
            cumulative,
            yerr=1.96 * cumulative_se,
            marker="o",
            ms=3,
            lw=1.1,
            capsize=2,
            label=f"n={n}",
        )
        axes[1].errorbar(
            time,
            interval,
            yerr=1.96 * interval_se,
            marker="o",
            ms=3,
            lw=1.1,
            capsize=2,
            label=f"n={n}",
        )
    axes[0].set_title(r"cumulative mean current over $[0,T]$")
    axes[1].set_title(r"last-interval mean current")
    for axis in axes:
        axis.set_xlabel("T")
        axis.set_ylabel("mean action current")
        axis.grid(True, alpha=0.25)
        axis.legend(frameon=False)
    fig.suptitle("No-burn-in current relaxation diagnostic")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out_dir / "burnin_flux_timeseries.png", dpi=200)
    fig.savefig(out_dir / "burnin_flux_timeseries.pdf")
    plt.close(fig)


def tau_values_from_blocks(path: Path) -> tuple[int, np.ndarray]:
    summary_path = path.with_name(path.name.replace("_blocks.csv", "_summary.csv"))
    summary = read_csv(summary_path)[0]
    n = int(summary["n"])
    tau_block = float(summary["tau_block"])
    raw = np.genfromtxt(path, delimiter=",", names=True)
    names = raw.dtype.names or ()
    block_names = [name for name in names if name.startswith("block_")]
    blocks = np.column_stack([np.atleast_1d(raw[name]) for name in block_names])
    taus = tau_block * np.arange(1, blocks.shape[1] + 1)
    prefix_values = np.cumsum(blocks, axis=1) / np.arange(1, blocks.shape[1] + 1)
    return n, np.column_stack([taus, prefix_values.T]).T


def load_tau_prefix(path: Path) -> tuple[int, np.ndarray, list[float]]:
    summary_path = path.with_name(path.name.replace("_blocks.csv", "_summary.csv"))
    summary = read_csv(summary_path)[0]
    n = int(summary["n"])
    tau_block = float(summary["tau_block"])
    raw = np.genfromtxt(path, delimiter=",", names=True)
    names = raw.dtype.names or ()
    block_names = [name for name in names if name.startswith("block_")]
    blocks = np.column_stack([np.atleast_1d(raw[name]) for name in block_names])
    values_by_tau = []
    taus = []
    cumulative = np.zeros(blocks.shape[0])
    for k in range(blocks.shape[1]):
        cumulative += blocks[:, k]
        taus.append((k + 1) * tau_block)
        values_by_tau.append(cumulative / float(k + 1))
    return n, np.column_stack(values_by_tau), taus


def analyze_tau_tails(
    exp_dir: Path,
    out_dir: Path,
    threshold_step: float,
    threshold_max: float,
    fit_low: float,
    fit_high: float,
) -> None:
    paths = sorted(exp_dir.glob("tau_n*_blocks.csv"))
    if not paths:
        return

    thresholds = np.arange(
        threshold_step,
        threshold_max + 0.5 * threshold_step,
        threshold_step,
    )
    surface_rows: list[dict] = []
    slope_rows: list[dict] = []

    for path in paths:
        n, prefix_values, taus = load_tau_prefix(path)
        fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4), sharex=True)
        cmap = plt.get_cmap("plasma")

        for tau_index, tau in enumerate(taus):
            values = prefix_values[:, tau_index]
            probs = np.array([np.mean(values > A) for A in thresholds])
            finite = probs > 0
            color = cmap(tau_index / max(1, len(taus) - 1))

            for A, prob in zip(thresholds, probs):
                surface_rows.append(
                    {
                        "n": n,
                        "tau": tau,
                        "A": A,
                        "samples": values.size,
                        "probability": prob,
                        "log_probability": math.log(prob)
                        if prob > 0
                        else float("-inf"),
                        "minus_log_probability_over_tau": -math.log(prob) / tau
                        if prob > 0
                        else float("inf"),
                    }
                )

            if np.any(finite):
                axes[0].plot(
                    thresholds[finite],
                    np.log(probs[finite]),
                    lw=1.2,
                    color=color,
                    label=rf"$\tau={tau:g}$",
                )
                axes[1].plot(
                    thresholds[finite],
                    -np.log(probs[finite]) / tau,
                    lw=1.2,
                    color=color,
                    label=rf"$\tau={tau:g}$",
                )

            fit_mask = finite & (probs >= fit_low) & (probs <= fit_high)
            if int(np.sum(fit_mask)) >= 3:
                slope, intercept = np.polyfit(
                    thresholds[fit_mask], np.log(probs[fit_mask]), 1
                )
                pred = intercept + slope * thresholds[fit_mask]
                residual = float(np.sum((np.log(probs[fit_mask]) - pred) ** 2))
                total = float(
                    np.sum(
                        (np.log(probs[fit_mask]) - np.mean(np.log(probs[fit_mask])))
                        ** 2
                    )
                )
                r2 = 1.0 - residual / total if total > 0 else float("nan")
                slope_rows.append(
                    {
                        "n": n,
                        "tau": tau,
                        "samples": values.size,
                        "mean": float(np.mean(values)),
                        "std": float(np.std(values, ddof=1)),
                        "lambda": -float(slope),
                        "intercept": float(intercept),
                        "r_squared": r2,
                        "fit_points": int(np.sum(fit_mask)),
                        "fit_A_min": float(np.min(thresholds[fit_mask])),
                        "fit_A_max": float(np.max(thresholds[fit_mask])),
                    }
                )

        axes[0].set_title(rf"$n={n}$: $\log P[\bar J_\tau>A]$")
        axes[0].set_ylabel("log survival probability")
        axes[1].set_title(rf"$n={n}$: rate proxy")
        axes[1].set_ylabel(r"$-\tau^{-1}\log P[\bar J_\tau>A]$")
        for axis in axes:
            axis.set_xlabel("threshold A")
            axis.grid(True, alpha=0.25)
            axis.legend(frameon=False, fontsize=7, ncol=2)
        fig.tight_layout()
        fig.savefig(out_dir / f"tau_tail_survival_n{n}.png", dpi=200)
        fig.savefig(out_dir / f"tau_tail_survival_n{n}.pdf")
        plt.close(fig)

        selected_As = [0.02, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20]
        fig, axis = plt.subplots(figsize=(6.0, 4.2))
        for A in selected_As:
            rows = [
                row
                for row in surface_rows
                if row["n"] == n and math.isclose(row["A"], A, abs_tol=1e-12)
            ]
            rows = [row for row in rows if row["probability"] > 0]
            if len(rows) >= 2:
                axis.plot(
                    [row["tau"] for row in rows],
                    [row["log_probability"] for row in rows],
                    marker="o",
                    ms=3,
                    lw=1.1,
                    label=f"A={A:g}",
                )
        axis.set_xlabel(r"$\tau$")
        axis.set_ylabel(r"$\log P[\bar J_\tau>A]$")
        axis.set_title(rf"$n={n}$: fixed-threshold dependence on $\tau$")
        axis.grid(True, alpha=0.25)
        axis.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(out_dir / f"tau_logprob_vs_tau_n{n}.png", dpi=200)
        fig.savefig(out_dir / f"tau_logprob_vs_tau_n{n}.pdf")
        plt.close(fig)

    if surface_rows:
        write_csv(
            out_dir / "tau_survival_surface.csv",
            surface_rows,
            [
                "n",
                "tau",
                "A",
                "samples",
                "probability",
                "log_probability",
                "minus_log_probability_over_tau",
            ],
        )
    if slope_rows:
        write_csv(
            out_dir / "tau_survival_slope_fits.csv",
            slope_rows,
            [
                "n",
                "tau",
                "samples",
                "mean",
                "std",
                "lambda",
                "intercept",
                "r_squared",
                "fit_points",
                "fit_A_min",
                "fit_A_max",
            ],
        )


def main() -> None:
    args = parse_args()
    exp_dir = args.experiment_dir
    out_dir = exp_dir / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_transient_profiles(exp_dir, out_dir)
    plot_transient_flux(exp_dir, out_dir)
    analyze_tau_tails(
        exp_dir,
        out_dir,
        args.threshold_step,
        args.threshold_max,
        args.survival_fit_low,
        args.survival_fit_high,
    )


if __name__ == "__main__":
    main()
