#!/usr/bin/env python3
"""Analyze finite-time averaged action-current distributions.

Each canonical simulation stores four disjoint equal blocks per trajectory.
This script compares averaging windows of one, two, and four blocks while
preserving trajectory clusters in the bootstrap.  It intentionally reports
descriptive finite-time statistics rather than asserting an asymptotic
large-deviation law from one window.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("summaries", nargs="+", type=Path)
    parser.add_argument("--bootstrap", type=int, default=2_000)
    parser.add_argument("--seed", type=int, default=20260619)
    parser.add_argument("--output-prefix", type=Path, required=True)
    return parser.parse_args()


def moments(values: np.ndarray) -> dict[str, float]:
    mean = float(np.mean(values))
    centered = values - mean
    variance = float(np.mean(centered**2))
    sd = math.sqrt(variance)
    if sd > 0:
        skewness = float(np.mean((centered / sd) ** 3))
        excess_kurtosis = float(np.mean((centered / sd) ** 4) - 3.0)
    else:
        skewness = float("nan")
        excess_kurtosis = float("nan")
    return {
        "mean": mean,
        "sd": sd,
        "skewness": skewness,
        "excess_kurtosis": excess_kurtosis,
        "negative_fraction": float(np.mean(values < 0.0)),
    }


def read_run(summary_path: Path) -> dict:
    with summary_path.open(newline="") as stream:
        summary_rows = list(csv.DictReader(stream))
    if len(summary_rows) != 1:
        raise ValueError(f"{summary_path}: expected one row")
    summary = summary_rows[0]
    suffix = "_summary.csv"
    if not summary_path.name.endswith(suffix):
        raise ValueError(f"{summary_path}: expected filename ending {suffix}")
    samples_path = summary_path.with_name(
        summary_path.name[: -len(suffix)] + "_samples.csv"
    )
    raw = np.genfromtxt(samples_path, delimiter=",", names=True)
    block = np.column_stack(
        [
            np.atleast_1d(raw[f"block_{index}_action_current"])
            for index in range(4)
        ]
    ).astype(float)
    first_half = np.atleast_1d(raw["first_half_action_current"]).astype(float)
    second_half = np.atleast_1d(raw["second_half_action_current"]).astype(float)
    full = np.atleast_1d(raw["time_averaged_action_current"]).astype(float)
    if not np.allclose(first_half, np.mean(block[:, :2], axis=1)):
        raise ValueError(f"{samples_path}: first-half column mismatch")
    if not np.allclose(second_half, np.mean(block[:, 2:], axis=1)):
        raise ValueError(f"{samples_path}: second-half column mismatch")
    if not np.allclose(full, np.mean(block, axis=1)):
        raise ValueError(f"{samples_path}: full-window column mismatch")
    return {
        "n": int(summary["n"]),
        "dt": float(summary["dt"]),
        "block_duration": float(summary["current_block_duration"]),
        "block": block,
        "samples_path": str(samples_path),
        "summary_path": str(summary_path),
    }


def window_arrays(run: dict) -> list[tuple[float, np.ndarray]]:
    block = run["block"]
    duration = run["block_duration"]
    return [
        (duration, block.reshape(-1)),
        (
            2.0 * duration,
            np.column_stack(
                [np.mean(block[:, :2], axis=1), np.mean(block[:, 2:], axis=1)]
            ).reshape(-1),
        ),
        (4.0 * duration, np.mean(block, axis=1)),
    ]


def cluster_bootstrap(
    block: np.ndarray, duration_multiple: int, count: int, rng: np.random.Generator
) -> np.ndarray:
    trajectories = block.shape[0]
    output = np.empty((count, 4), dtype=float)
    for replicate in range(count):
        selected = block[rng.integers(0, trajectories, size=trajectories)]
        if duration_multiple == 1:
            values = selected.reshape(-1)
        elif duration_multiple == 2:
            values = np.column_stack(
                [
                    np.mean(selected[:, :2], axis=1),
                    np.mean(selected[:, 2:], axis=1),
                ]
            ).reshape(-1)
        elif duration_multiple == 4:
            values = np.mean(selected, axis=1)
        else:
            raise ValueError("duration_multiple must be 1, 2, or 4")
        stat = moments(values)
        output[replicate] = [
            stat["sd"],
            stat["skewness"],
            stat["excess_kurtosis"],
            stat["negative_fraction"],
        ]
    return output


def main() -> None:
    args = parse_args()
    runs = sorted(
        [read_run(path) for path in args.summaries],
        key=lambda run: run["n"],
    )
    rng = np.random.default_rng(args.seed)
    rows: list[dict] = []

    for run in runs:
        for multiple, (duration, values) in zip(
            (1, 2, 4), window_arrays(run)
        ):
            stat = moments(values)
            bootstrap = cluster_bootstrap(
                run["block"], multiple, args.bootstrap, rng
            )
            lower = np.quantile(bootstrap, 0.025, axis=0)
            upper = np.quantile(bootstrap, 0.975, axis=0)
            rows.append(
                {
                    "n": run["n"],
                    "dt": run["dt"],
                    "window": duration,
                    "observations": values.size,
                    "independent_trajectory_clusters": run["block"].shape[0],
                    "mean": stat["mean"],
                    "sd": stat["sd"],
                    "sd_ci_lower": lower[0],
                    "sd_ci_upper": upper[0],
                    "skewness": stat["skewness"],
                    "skewness_ci_lower": lower[1],
                    "skewness_ci_upper": upper[1],
                    "excess_kurtosis": stat["excess_kurtosis"],
                    "kurtosis_ci_lower": lower[2],
                    "kurtosis_ci_upper": upper[2],
                    "negative_fraction": stat["negative_fraction"],
                    "negative_fraction_ci_lower": lower[3],
                    "negative_fraction_ci_upper": upper[3],
                    "window_times_variance": duration * stat["sd"] ** 2,
                }
            )

    output_prefix = args.output_prefix
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    table_path = output_prefix.with_name(
        output_prefix.name + "_window_statistics.csv"
    )
    with table_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    colors = ["#1f4e79", "#6c8ebf", "#b03a2e"]
    for run in runs:
        fig, axis = plt.subplots(figsize=(6.0, 4.3))
        for color, (duration, values) in zip(colors, window_arrays(run)):
            standardized = (values - np.mean(values)) / np.std(values)
            axis.hist(
                standardized,
                bins=45,
                density=True,
                histtype="step",
                linewidth=1.5,
                color=color,
                label=rf"$\tau={duration:g}$",
            )
        grid = np.linspace(-4.0, 4.0, 400)
        normal = np.exp(-0.5 * grid**2) / math.sqrt(2.0 * math.pi)
        axis.plot(grid, normal, "k--", linewidth=1.2, label="standard normal")
        axis.set_xlabel(r"standardized finite-time current")
        axis.set_ylabel("density")
        axis.set_title(f"n={run['n']}")
        axis.legend(frameon=False)
        axis.grid(True, alpha=0.2)
        fig.tight_layout()
        fig.savefig(
            output_prefix.with_name(
                output_prefix.name + f"_n{run['n']}_windows.pdf"
            )
        )
        fig.savefig(
            output_prefix.with_name(
                output_prefix.name + f"_n{run['n']}_windows.png"
            ),
            dpi=240,
        )
        plt.close(fig)

    fig, axis = plt.subplots(figsize=(6.0, 4.3))
    for run in runs:
        selected = [row for row in rows if row["n"] == run["n"]]
        axis.plot(
            [row["window"] for row in selected],
            [row["window_times_variance"] for row in selected],
            "o-",
            label=f"n={run['n']}",
        )
    axis.set_xscale("log", base=2)
    axis.set_xlabel(r"averaging window $\tau$")
    axis.set_ylabel(r"$\tau\,\mathrm{Var}(\overline{J}_\tau)$")
    axis.grid(True, alpha=0.2)
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(
        output_prefix.with_name(output_prefix.name + "_variance_scaling.pdf")
    )
    fig.savefig(
        output_prefix.with_name(output_prefix.name + "_variance_scaling.png"),
        dpi=240,
    )

    summary = {
        "runs": [
            {
                "n": run["n"],
                "dt": run["dt"],
                "summary_path": run["summary_path"],
                "samples_path": run["samples_path"],
            }
            for run in runs
        ],
        "bootstrap_replicates": args.bootstrap,
        "bootstrap_seed": args.seed,
        "bootstrap_unit": "trajectory cluster",
        "claim_scope": (
            "Finite-time descriptive statistics at three disjoint averaging "
            "scales. No asymptotic exponential-tail or fluctuation-theorem "
            "claim is inferred from these data alone."
        ),
    }
    summary_path = output_prefix.with_name(
        output_prefix.name + "_window_analysis.json"
    )
    with summary_path.open("w") as stream:
        json.dump(summary, stream, indent=2)
        stream.write("\n")

    print(f"wrote {table_path}")
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()
