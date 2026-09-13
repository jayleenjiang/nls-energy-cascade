#!/usr/bin/env python3
"""Validate and analyze canonical NLS action-current runs.

The script treats each trajectory's finite-time average as one independent
replicate.  It verifies summary files against raw samples, fits

    E[J_n] = C n^alpha

on a selected timestep, and obtains a Monte Carlo confidence interval for the
exponent by resampling trajectories within each chain length.  A separate table
reports timestep sensitivity when matching runs are available.
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


MODEL_VERSION = "gibbs-canonical-v1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "summaries",
        nargs="+",
        type=Path,
        help="Canonical *_summary.csv files.",
    )
    parser.add_argument(
        "--primary-dt",
        type=float,
        required=True,
        help="Timestep used for the primary scaling fit.",
    )
    parser.add_argument("--bootstrap", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=20260619)
    parser.add_argument("--output-prefix", type=Path, required=True)
    return parser.parse_args()


def read_one_summary(path: Path) -> dict:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 1:
        raise ValueError(f"{path}: expected exactly one summary row")
    row = rows[0]
    if row["model_version"] != MODEL_VERSION:
        raise ValueError(
            f"{path}: model_version={row['model_version']!r}, "
            f"expected {MODEL_VERSION!r}"
        )

    numeric_int = {
        "n",
        "batches",
        "lanes",
        "n_trajectories",
        "bond",
        "seed",
        "threads",
        "projection_count",
    }
    numeric_float = {
        "T1",
        "Tn",
        "gamma",
        "dt",
        "burnin",
        "measure",
        "mean_action_current",
        "sample_sd",
        "standard_error",
        "normal95_ci_lower",
        "normal95_ci_upper",
        "current_block_duration",
        "mean_first_half_current",
        "mean_second_half_current",
        "mean_second_minus_first",
        "paired_difference_se",
        "action_floor",
        "projection_rate",
        "projection_events_per_trajectory_time",
        "elapsed_seconds",
    }
    for key in numeric_int:
        row[key] = int(row[key])
    for key in numeric_float:
        row[key] = float(row[key])

    suffix = "_summary.csv"
    if not path.name.endswith(suffix):
        raise ValueError(f"{path}: filename must end with {suffix}")
    sample_path = path.with_name(path.name[: -len(suffix)] + "_samples.csv")
    raw = np.genfromtxt(sample_path, delimiter=",", names=True)
    samples = np.atleast_1d(raw["time_averaged_action_current"]).astype(float)
    if samples.size != row["n_trajectories"]:
        raise ValueError(
            f"{sample_path}: {samples.size} rows, "
            f"summary claims {row['n_trajectories']}"
        )

    mean = float(np.mean(samples))
    sd = float(np.std(samples, ddof=1))
    se = sd / math.sqrt(samples.size)
    tolerances = {
        "mean_action_current": mean,
        "sample_sd": sd,
        "standard_error": se,
    }
    for key, recomputed in tolerances.items():
        if not math.isclose(row[key], recomputed, rel_tol=2e-12, abs_tol=2e-14):
            raise ValueError(
                f"{path}: {key}={row[key]:.17g}, "
                f"raw-sample value={recomputed:.17g}"
            )
    row["samples"] = samples
    row["summary_path"] = str(path)
    row["sample_path"] = str(sample_path)
    return row


def same_float(a: float, b: float) -> bool:
    return math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-15)


def select_primary(rows: list[dict], primary_dt: float) -> list[dict]:
    selected = [row for row in rows if same_float(row["dt"], primary_dt)]
    if not selected:
        raise ValueError(f"no runs found at primary dt={primary_dt}")
    selected.sort(key=lambda row: row["n"])
    seen: set[int] = set()
    for row in selected:
        if row["n"] in seen:
            raise ValueError(
                f"multiple primary-dt runs found for n={row['n']}; "
                "combine replicates explicitly before fitting"
            )
        seen.add(row["n"])
        if row["mean_action_current"] <= 0:
            raise ValueError(
                f"n={row['n']}: log scaling requires positive mean current"
            )
    if len(selected) < 3:
        raise ValueError("at least three chain lengths are required")
    return selected


def linear_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    slope, intercept = np.polyfit(x, y, 1)
    prediction = intercept + slope * x
    residual = float(np.sum((y - prediction) ** 2))
    total = float(np.sum((y - np.mean(y)) ** 2))
    r_squared = 1.0 - residual / total if total > 0 else float("nan")
    return float(slope), float(intercept), r_squared


def bootstrap_exponent(
    rows: list[dict], count: int, seed: int
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    log_n = np.log([row["n"] for row in rows])
    exponents = np.empty(count, dtype=float)
    for index in range(count):
        means = []
        for row in rows:
            samples = row["samples"]
            draw = rng.choice(samples, size=samples.size, replace=True)
            means.append(float(np.mean(draw)))
        if min(means) <= 0:
            exponents[index] = np.nan
        else:
            exponents[index] = np.polyfit(log_n, np.log(means), 1)[0]
    return exponents[np.isfinite(exponents)]


def timestep_comparisons(rows: list[dict]) -> list[dict]:
    comparisons: list[dict] = []
    groups: dict[tuple, list[dict]] = {}
    for row in rows:
        key = (
            row["n"],
            row["T1"],
            row["Tn"],
            row["burnin"],
            row["measure"],
        )
        groups.setdefault(key, []).append(row)
    for key, group in groups.items():
        group.sort(key=lambda row: row["dt"])
        for fine, coarse in zip(group, group[1:]):
            fine_mean = fine["mean_action_current"]
            coarse_mean = coarse["mean_action_current"]
            difference = coarse_mean - fine_mean
            pooled_se = math.sqrt(
                fine["standard_error"] ** 2 + coarse["standard_error"] ** 2
            )
            comparisons.append(
                {
                    "n": key[0],
                    "T1": key[1],
                    "Tn": key[2],
                    "burnin": key[3],
                    "measure": key[4],
                    "fine_dt": fine["dt"],
                    "coarse_dt": coarse["dt"],
                    "fine_mean": fine_mean,
                    "coarse_mean": coarse_mean,
                    "relative_difference": difference / fine_mean,
                    "difference_over_pooled_se": (
                        difference / pooled_se if pooled_se > 0 else float("nan")
                    ),
                }
            )
    return comparisons


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    rows = [read_one_summary(path) for path in args.summaries]
    primary = select_primary(rows, args.primary_dt)

    log_n = np.log([row["n"] for row in primary])
    log_mean = np.log([row["mean_action_current"] for row in primary])
    exponent, log_prefactor, r_squared = linear_fit(log_n, log_mean)
    bootstrap = bootstrap_exponent(primary, args.bootstrap, args.seed)
    if bootstrap.size < 0.99 * args.bootstrap:
        raise RuntimeError(
            "more than 1% of bootstrap fits had a non-positive resampled mean"
        )
    exponent_ci = np.quantile(bootstrap, [0.025, 0.975])

    output_prefix = args.output_prefix
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    run_table: list[dict] = []
    for row in sorted(rows, key=lambda item: (item["dt"], item["n"])):
        stationarity_z = (
            row["mean_second_minus_first"] / row["paired_difference_se"]
            if row["paired_difference_se"] > 0
            else float("nan")
        )
        run_record = {
            key: row[key]
            for key in (
                "n",
                "T1",
                "Tn",
                "dt",
                "burnin",
                "measure",
                "n_trajectories",
                "mean_action_current",
                "sample_sd",
                "standard_error",
                "normal95_ci_lower",
                "normal95_ci_upper",
                "current_block_duration",
                "mean_first_half_current",
                "mean_second_half_current",
                "mean_second_minus_first",
                "paired_difference_se",
                "action_floor",
                "projection_count",
                "projection_rate",
                "projection_events_per_trajectory_time",
                "seed",
            )
        }
        run_record["stationarity_z"] = stationarity_z
        run_table.append(run_record)
    write_csv(
        output_prefix.with_name(output_prefix.name + "_runs.csv"),
        run_table,
        list(run_table[0]),
    )

    comparisons = timestep_comparisons(rows)
    comparison_fields = [
        "n",
        "T1",
        "Tn",
        "burnin",
        "measure",
        "fine_dt",
        "coarse_dt",
        "fine_mean",
        "coarse_mean",
        "relative_difference",
        "difference_over_pooled_se",
    ]
    write_csv(
        output_prefix.with_name(output_prefix.name + "_dt_sensitivity.csv"),
        comparisons,
        comparison_fields,
    )

    result = {
        "model_version": MODEL_VERSION,
        "primary_dt": args.primary_dt,
        "chain_lengths": [row["n"] for row in primary],
        "exponent": exponent,
        "prefactor": math.exp(log_prefactor),
        "r_squared_log_fit": r_squared,
        "bootstrap_replicates_requested": args.bootstrap,
        "bootstrap_replicates_valid": int(bootstrap.size),
        "bootstrap_seed": args.seed,
        "exponent_normalized_95_ci": [
            float(exponent_ci[0]),
            float(exponent_ci[1]),
        ],
        "stationarity_max_abs_z": max(
            abs(row["stationarity_z"]) for row in run_table
        ),
        "stationarity_flagged_runs_abs_z_ge_2": [
            {
                "n": row["n"],
                "dt": row["dt"],
                "stationarity_z": row["stationarity_z"],
            }
            for row in run_table
            if abs(row["stationarity_z"]) >= 2.0
        ],
        "interpretation_scope": (
            "Monte Carlo uncertainty conditional on the selected chain lengths, "
            "time window, timestep, and power-law model; it does not include "
            "finite-size or discretization-model uncertainty."
        ),
    }
    with output_prefix.with_name(
        output_prefix.name + "_scaling.json"
    ).open("w") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")

    n_values = np.array([row["n"] for row in primary], dtype=float)
    means = np.array(
        [row["mean_action_current"] for row in primary], dtype=float
    )
    errors = np.array(
        [1.959963984540054 * row["standard_error"] for row in primary],
        dtype=float,
    )
    grid = np.geomspace(n_values.min(), n_values.max(), 200)
    fit = math.exp(log_prefactor) * grid**exponent

    fig, axis = plt.subplots(figsize=(5.8, 4.3))
    axis.errorbar(
        n_values,
        means,
        yerr=errors,
        fmt="o",
        capsize=3,
        color="#1f4e79",
        label="trajectory means (normal 95% CI)",
    )
    axis.plot(
        grid,
        fit,
        color="#b03a2e",
        label=(
            rf"$\langle J\rangle={math.exp(log_prefactor):.3g}"
            rf"n^{{{exponent:.3f}}}$"
        ),
    )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("chain length $n$")
    axis.set_ylabel("mean action current")
    axis.legend(frameon=False)
    axis.grid(True, which="both", alpha=0.2)
    fig.tight_layout()
    figure_path = output_prefix.with_name(output_prefix.name + "_scaling.pdf")
    fig.savefig(figure_path)
    fig.savefig(
        output_prefix.with_name(output_prefix.name + "_scaling.png"), dpi=240
    )

    print(json.dumps(result, indent=2))
    print(f"wrote {figure_path}")


if __name__ == "__main__":
    main()
