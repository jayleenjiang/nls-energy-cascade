#!/usr/bin/env python3
"""Fit descriptive large-deviation time scaling of action-current tails.

For each chain length, threshold, and tail, fit

    log P_t = intercept - t I(A)

using only raw probabilities with enough observed events and lying in a
declared rare-tail window.  Different averaging times are aggregated from the
same underlying streams, so these OLS fits are descriptive rate-function
proxies rather than independent-sample hypothesis tests.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("survival_csv", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--minimum-raw-count", type=int, default=20)
    parser.add_argument("--maximum-probability", type=float, default=0.2)
    parser.add_argument("--minimum-time-points", type=int, default=4)
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("time-scaling analysis produced no qualified fits")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fit_line(time: np.ndarray, log_probability: np.ndarray) -> tuple[float, float, float]:
    slope, intercept = np.polyfit(time, log_probability, 1)
    fitted = slope * time + intercept
    residual = float(np.sum((log_probability - fitted) ** 2))
    total = float(np.sum((log_probability - np.mean(log_probability)) ** 2))
    r_squared = 1.0 - residual / total if total > 0.0 else float("nan")
    return float(-slope), float(intercept), r_squared


def plot_chain(output_dir: Path, n: int, rows: list[dict[str, object]]) -> None:
    selected = [row for row in rows if int(row["n"]) == n]
    if not selected:
        return
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.1))
    for tail, marker, label in [
        ("plus", "o", r"$P[J_t\geq A]$"),
        ("minus", "s", r"$P[J_t\leq-A]$"),
    ]:
        tail_rows = sorted(
            (row for row in selected if row["tail"] == tail),
            key=lambda row: float(row["A"]),
        )
        if not tail_rows:
            continue
        threshold = np.asarray([float(row["A"]) for row in tail_rows])
        rate = np.asarray([float(row["rate_proxy"]) for row in tail_rows])
        r_squared = np.asarray([float(row["r_squared"]) for row in tail_rows])
        axes[0].plot(threshold, rate, marker=marker, markersize=4, label=label)
        axes[1].plot(threshold, r_squared, marker=marker, markersize=4, label=label)
    axes[0].set_xlabel(r"threshold $A$")
    axes[0].set_ylabel(r"rate proxy $I_\pm(A)$")
    axes[0].grid(alpha=0.2)
    axes[0].legend(fontsize=8)
    axes[1].axhline(0.98, color="k", linestyle="--", linewidth=1.0, label=r"$R^2=0.98$")
    axes[1].set_xlabel(r"threshold $A$")
    axes[1].set_ylabel(r"linearity of $\log P$ versus $t$: $R^2$")
    axes[1].set_ylim(0.9, 1.005)
    axes[1].grid(alpha=0.2)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output_dir / f"action_tail_time_scaling_n{n}.png", dpi=220)
    fig.savefig(output_dir / f"action_tail_time_scaling_n{n}.pdf")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if args.minimum_raw_count <= 0:
        raise ValueError("minimum raw count must be positive")
    if not 0.0 < args.maximum_probability < 1.0:
        raise ValueError("maximum probability must lie in (0,1)")
    if args.minimum_time_points < 3:
        raise ValueError("at least three time points are required")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    grouped: dict[tuple[int, float], list[dict[str, str]]] = defaultdict(list)
    for row in read_rows(args.survival_csv):
        grouped[(int(row["n"]), float(row["A"]))].append(row)

    output: list[dict[str, object]] = []
    for (n, threshold), rows in sorted(grouped.items()):
        for tail, count_column, probability_column in [
            ("plus", "plus_count", "p_plus_raw"),
            ("minus", "minus_count", "p_minus_raw"),
        ]:
            qualified = sorted(
                (
                    row
                    for row in rows
                    if int(row[count_column]) >= args.minimum_raw_count
                    and 0.0 < float(row[probability_column]) <= args.maximum_probability
                ),
                key=lambda row: float(row["tau"]),
            )
            if len(qualified) < args.minimum_time_points:
                continue
            time = np.asarray([float(row["tau"]) for row in qualified])
            probability = np.asarray(
                [float(row[probability_column]) for row in qualified]
            )
            rate, intercept, r_squared = fit_line(time, np.log(probability))
            output.append(
                {
                    "n": n,
                    "A": threshold,
                    "tail": tail,
                    "time_points": len(qualified),
                    "t_min": float(np.min(time)),
                    "t_max": float(np.max(time)),
                    "minimum_raw_count": args.minimum_raw_count,
                    "maximum_probability": args.maximum_probability,
                    "rate_proxy": rate,
                    "intercept": intercept,
                    "r_squared": r_squared,
                    "minimum_observed_count": min(
                        int(row[count_column]) for row in qualified
                    ),
                    "maximum_included_probability": max(probability),
                }
            )

    write_rows(args.output_dir / "action_tail_time_scaling.csv", output)
    for n in sorted({int(row["n"]) for row in output}):
        plot_chain(args.output_dir, n, output)
    print(f"wrote action-tail time-scaling analysis to {args.output_dir}")


if __name__ == "__main__":
    main()
