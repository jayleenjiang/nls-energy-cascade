#!/usr/bin/env python3
"""Audit a bidirectional discrete path-ratio FT experiment.

The forward and reverse files are produced by NLS_discrete_path_ft.cpp.  This
script evaluates the two integral fluctuation relations and the bidirectional
Crooks histogram relation

    log p_F(Sigma=s) / p_R(Sigma_R=-s) = s.

Only raw bins with adequate counts on both sides are fitted.  IFT uncertainty
is obtained by bootstrap resampling independent trajectory groups, while the
exponential-weight ESS and maximum weight fraction diagnose rare-event
support.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("forward_csv", type=Path)
    parser.add_argument("reverse_csv", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bins", type=int, default=80)
    parser.add_argument("--min-count", type=int, default=50)
    parser.add_argument("--groups", type=int, default=128)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260828)
    return parser.parse_args()


def load_csv(path: Path) -> dict[str, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=None,
                         encoding="utf-8")
    if data.ndim == 0:
        data = np.asarray([data], dtype=data.dtype)
    return {name: np.asarray(data[name]) for name in data.dtype.names or ()}


def log_mean_exp(log_values: np.ndarray) -> float:
    maximum = float(np.max(log_values))
    return maximum + math.log(float(np.mean(np.exp(log_values - maximum))))


def ift_metrics(sigma: np.ndarray, groups: int, bootstrap: int,
                rng: np.random.Generator) -> dict[str, float | list[float]]:
    log_weights = -np.asarray(sigma, dtype=float)
    maximum = float(np.max(log_weights))
    scaled = np.exp(log_weights - maximum)
    scaled_sum = float(np.sum(scaled))
    ess = scaled_sum * scaled_sum / float(np.sum(scaled * scaled))
    max_fraction = float(np.max(scaled) / scaled_sum)
    log_ift = log_mean_exp(log_weights)

    group_count = min(groups, len(log_weights))
    index_groups = np.array_split(np.arange(len(log_weights)), group_count)
    group_sums = np.asarray([float(np.sum(np.exp(log_weights[g] - maximum)))
                             for g in index_groups])
    group_sizes = np.asarray([len(g) for g in index_groups], dtype=float)
    draws = rng.integers(0, group_count, size=(bootstrap, group_count))
    bootstrap_sums = np.sum(group_sums[draws], axis=1)
    bootstrap_sizes = np.sum(group_sizes[draws], axis=1)
    bootstrap_logs = maximum + np.log(bootstrap_sums / bootstrap_sizes)
    low, high = np.quantile(bootstrap_logs, [0.025, 0.975])

    return {
        "samples": int(len(sigma)),
        "mean_sigma": float(np.mean(sigma)),
        "std_sigma": float(np.std(sigma, ddof=1)),
        "log_mean_exp_minus_sigma": float(log_ift),
        "log_mean_exp_ci95": [float(low), float(high)],
        "exponential_weight_ess": float(ess),
        "maximum_weight_fraction": max_fraction,
        "ift_ci_includes_zero": bool(low <= 0.0 <= high),
    }


def weighted_line(x: np.ndarray, y: np.ndarray,
                  standard_error: np.ndarray) -> dict[str, float]:
    design = np.column_stack([np.ones_like(x), x])
    weight = 1.0 / np.square(standard_error)
    normal = design.T @ (weight[:, None] * design)
    covariance = np.linalg.inv(normal)
    coefficients = covariance @ (design.T @ (weight * y))
    fitted = design @ coefficients
    residual = y - fitted
    chi_square = float(np.sum(np.square(residual / standard_error)))
    degrees = max(1, len(x) - 2)
    return {
        "intercept": float(coefficients[0]),
        "slope": float(coefficients[1]),
        "intercept_se": float(math.sqrt(covariance[0, 0])),
        "slope_se": float(math.sqrt(covariance[1, 1])),
        "chi_square": chi_square,
        "degrees_of_freedom": int(degrees),
        "reduced_chi_square": chi_square / degrees,
        "identity_weighted_rmse": float(
            math.sqrt(np.average(np.square(y - x), weights=weight))
        ),
    }


def crooks_bins(forward: np.ndarray, reverse: np.ndarray, bins: int,
                 min_count: int) -> tuple[list[dict[str, float | int]],
                                          dict[str, float | int | bool]]:
    pooled = np.concatenate([np.abs(forward), np.abs(reverse)])
    limit = float(np.quantile(pooled, 0.9975))
    limit = max(limit, 1.0e-8)
    edges = np.linspace(-limit, limit, bins + 1)
    forward_counts, _ = np.histogram(forward, bins=edges)
    reverse_counts, _ = np.histogram(reverse, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    rows: list[dict[str, float | int]] = []
    accepted_x: list[float] = []
    accepted_y: list[float] = []
    accepted_se: list[float] = []
    for index, center in enumerate(centers):
        mirror = bins - 1 - index
        count_forward = int(forward_counts[index])
        count_reverse = int(reverse_counts[mirror])
        accepted = count_forward >= min_count and count_reverse >= min_count
        log_ratio = math.nan
        standard_error = math.nan
        if count_forward > 0 and count_reverse > 0:
            log_ratio = math.log(
                (count_forward / len(forward)) /
                (count_reverse / len(reverse))
            )
            standard_error = math.sqrt(1.0 / count_forward +
                                       1.0 / count_reverse)
        rows.append({
            "bin_left": float(edges[index]),
            "bin_right": float(edges[index + 1]),
            "bin_center": float(center),
            "forward_count": count_forward,
            "reverse_mirror_count": count_reverse,
            "log_probability_ratio": float(log_ratio),
            "standard_error": float(standard_error),
            "accepted": int(accepted),
        })
        if accepted:
            accepted_x.append(float(center))
            accepted_y.append(float(log_ratio))
            accepted_se.append(float(standard_error))

    if len(accepted_x) < 3:
        summary: dict[str, float | int | bool] = {
            "accepted_bins": len(accepted_x),
            "fit_available": False,
        }
    else:
        summary = {
            "accepted_bins": len(accepted_x),
            "fit_available": True,
            **weighted_line(np.asarray(accepted_x), np.asarray(accepted_y),
                            np.asarray(accepted_se)),
        }
    return rows, summary


def write_rows(path: Path, rows: list[dict[str, float | int]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def plot_results(output: Path, forward: np.ndarray, reverse: np.ndarray,
                 rows: list[dict[str, float | int]],
                 crooks: dict[str, float | int | bool]) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.3))
    pooled = np.concatenate([forward, -reverse])
    edges = np.linspace(float(np.quantile(pooled, 0.0025)),
                        float(np.quantile(pooled, 0.9975)), 90)
    axes[0].hist(forward, bins=edges, density=True, histtype="step",
                 linewidth=1.8, label=r"forward $\Sigma_F$")
    axes[0].hist(-reverse, bins=edges, density=True, histtype="step",
                 linewidth=1.8, label=r"reverse $-\Sigma_R$")
    axes[0].set_yscale("log")
    axes[0].set_xlabel(r"path entropy $s$")
    axes[0].set_ylabel("density")
    axes[0].legend(frameon=False)

    accepted = [row for row in rows if row["accepted"]]
    x = np.asarray([row["bin_center"] for row in accepted], dtype=float)
    y = np.asarray([row["log_probability_ratio"] for row in accepted],
                   dtype=float)
    error = np.asarray([row["standard_error"] for row in accepted],
                       dtype=float)
    axes[1].errorbar(x, y, yerr=error, fmt="o", markersize=3.5,
                     capsize=2, label="raw-bin ratio")
    if len(x):
        grid = np.linspace(float(np.min(x)), float(np.max(x)), 200)
        axes[1].plot(grid, grid, color="black", linestyle="--",
                     label=r"FT: $y=s$")
        if crooks.get("fit_available"):
            axes[1].plot(
                grid,
                float(crooks["intercept"]) + float(crooks["slope"]) * grid,
                label=(f"fit: {float(crooks['intercept']):.3f} + "
                       f"{float(crooks['slope']):.3f} s"),
            )
    axes[1].set_xlabel(r"entropy-bin center $s$")
    axes[1].set_ylabel(r"$\log[p_F(s)/p_R(-s)]$")
    axes[1].legend(frameon=False)
    figure.tight_layout()
    figure.savefig(output / "discrete_path_ft.png", dpi=220)
    figure.savefig(output / "discrete_path_ft.pdf")
    plt.close(figure)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    forward_data = load_csv(args.forward_csv)
    reverse_data = load_csv(args.reverse_csv)
    if not np.all(forward_data["finite"] == 1) or not np.all(
            reverse_data["finite"] == 1):
        raise RuntimeError("nonfinite trajectories present")

    forward = np.asarray(forward_data["sigma_total"], dtype=float)
    reverse = np.asarray(reverse_data["sigma_total"], dtype=float)
    rng = np.random.default_rng(args.seed)
    forward_ift = ift_metrics(forward, args.groups, args.bootstrap, rng)
    reverse_ift = ift_metrics(reverse, args.groups, args.bootstrap, rng)
    rows, crooks = crooks_bins(
        forward, reverse, args.bins, args.min_count)
    write_rows(args.output_dir / "crooks_bins.csv", rows)

    numerical = {}
    for label, data in (("forward", forward_data),
                        ("reverse", reverse_data)):
        numerical[label] = {
            "midpoint_failures": int(np.sum(data["midpoint_failures"])),
            "maximum_abs_energy_balance_error": float(
                np.max(np.abs(data["energy_balance_error"]))),
            "rms_energy_balance_error": float(
                np.sqrt(np.mean(np.square(data["energy_balance_error"])))),
            "mean_kernel_minus_heat": float(
                np.mean(data["kernel_minus_heat"])),
            "rms_kernel_minus_heat": float(
                np.sqrt(np.mean(np.square(data["kernel_minus_heat"])))),
        }

    gate = {
        "numerical_integrity": bool(
            numerical["forward"]["midpoint_failures"] == 0 and
            numerical["reverse"]["midpoint_failures"] == 0
        ),
        "forward_ift": bool(
            forward_ift["ift_ci_includes_zero"] and
            forward_ift["exponential_weight_ess"] >= 1000 and
            forward_ift["maximum_weight_fraction"] <= 0.01
        ),
        "reverse_ift": bool(
            reverse_ift["ift_ci_includes_zero"] and
            reverse_ift["exponential_weight_ess"] >= 1000 and
            reverse_ift["maximum_weight_fraction"] <= 0.01
        ),
        "crooks_histogram": bool(
            crooks.get("fit_available") and
            int(crooks["accepted_bins"]) >= 8 and
            abs(float(crooks["slope"]) - 1.0) <= max(
                0.1, 2.0 * float(crooks["slope_se"])) and
            abs(float(crooks["intercept"])) <= max(
                0.1, 2.0 * float(crooks["intercept_se"]))
        ),
    }
    gate["overall"] = bool(all(gate.values()))
    summary = {
        "forward_file": str(args.forward_csv),
        "reverse_file": str(args.reverse_csv),
        "analysis": {
            "bins": args.bins,
            "minimum_raw_count": args.min_count,
            "bootstrap_groups": args.groups,
            "bootstrap_replicates": args.bootstrap,
            "bootstrap_seed": args.seed,
        },
        "forward_ift": forward_ift,
        "reverse_ift": reverse_ift,
        "crooks": crooks,
        "numerical": numerical,
        "gates": gate,
    }
    (args.output_dir / "audit.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    plot_results(args.output_dir, forward, reverse, rows, crooks)

    with (args.output_dir / "audit.md").open("w") as stream:
        stream.write("# Discrete path-ratio FT audit\n\n")
        stream.write(f"Overall: **{'PASS' if gate['overall'] else 'BLOCKED'}**\n\n")
        stream.write("| Gate | Status |\n|---|---:|\n")
        for name, passed in gate.items():
            stream.write(f"| {name} | {'PASS' if passed else 'BLOCKED'} |\n")
        stream.write("\n## Integral fluctuation relation\n\n")
        stream.write("| ensemble | log mean exp(-Sigma) | 95% CI | ESS | max fraction |\n")
        stream.write("|---|---:|---:|---:|---:|\n")
        for label, values in (("forward", forward_ift),
                              ("reverse", reverse_ift)):
            ci = values["log_mean_exp_ci95"]
            stream.write(
                f"| {label} | {values['log_mean_exp_minus_sigma']:.6g} | "
                f"[{ci[0]:.6g}, {ci[1]:.6g}] | "
                f"{values['exponential_weight_ess']:.1f} | "
                f"{values['maximum_weight_fraction']:.4g} |\n"
            )
        stream.write("\n## Bidirectional Crooks histogram\n\n")
        if crooks.get("fit_available"):
            stream.write(
                f"Accepted bins: {crooks['accepted_bins']}; "
                f"slope = {crooks['slope']:.6g} +/- "
                f"{crooks['slope_se']:.3g}; intercept = "
                f"{crooks['intercept']:.6g} +/- "
                f"{crooks['intercept_se']:.3g}; identity weighted RMSE = "
                f"{crooks['identity_weighted_rmse']:.6g}.\n"
            )
        else:
            stream.write("Insufficient two-sided raw-bin support.\n")
    print(f"Wrote discrete path-ratio FT audit to {args.output_dir}")


if __name__ == "__main__":
    main()
