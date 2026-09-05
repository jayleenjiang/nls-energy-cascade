#!/usr/bin/env python3
"""Adaptive-range robustness check for fluctuation-symmetry slopes.

The production analyzer uses a predeclared fixed-width histogram.  This script
keeps that result intact and provides a second, explicitly labelled
robustness estimator.  It removes the primary estimator's discontinuity when
the overall first percentile crosses zero: the symmetric fit range is instead
set by a fixed quantile of the positive and negative-magnitude samples
separately.  Bins remain equal-width and symmetric, while their number adapts
to the rarer sign's effective sample size.
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

from analyze_entropy_ft import aggregate_blocks, integrated_autocorrelation_time, load_run


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("blocks", nargs="+", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--taus", default="20,40,60,80,100,120,140,160,180,200")
    parser.add_argument("--max-bins", type=int, default=60)
    parser.add_argument("--min-effective-count", type=float, default=50.0)
    parser.add_argument("--range-quantile", type=float, default=0.99)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--max-acf-lag", type=int, default=20)
    parser.add_argument("--seed", type=int, default=2026082801)
    return parser.parse_args()


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def weighted_fit(
    x: np.ndarray, y: np.ndarray, variance: np.ndarray
) -> tuple[float, float, float]:
    weights = 1.0 / np.sqrt(np.maximum(variance, np.finfo(float).tiny))
    design = np.column_stack([x, np.ones_like(x)])
    coefficients, *_ = np.linalg.lstsq(
        design * weights[:, None], y * weights, rcond=None
    )
    slope, intercept = coefficients
    fitted = slope * x + intercept
    weighted_mean = np.average(y, weights=weights * weights)
    residual = np.sum((weights * (y - fitted)) ** 2)
    total = np.sum((weights * (y - weighted_mean)) ** 2)
    r_squared = 1.0 - residual / total if total > 0.0 else float("nan")
    return float(slope), float(intercept), float(r_squared)


def interval(values: list[float]) -> tuple[float, float]:
    finite = np.asarray([value for value in values if math.isfinite(value)])
    if finite.size < 20:
        return float("nan"), float("nan")
    return tuple(float(value) for value in np.quantile(finite, [0.025, 0.975]))


def adaptive_symmetry(
    values: np.ndarray,
    tau: float,
    max_bins: int,
    min_effective_count: float,
    range_quantile: float,
    iat: float,
    bootstrap: int,
    rng: np.random.Generator,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    flattened = values.ravel()
    positive = flattened[flattened > 0.0]
    negative_magnitude = -flattened[flattened < 0.0]
    empty = {
        "adaptive_slope": float("nan"),
        "adaptive_intercept": float("nan"),
        "adaptive_r_squared": float("nan"),
        "adaptive_slope_ci_low": float("nan"),
        "adaptive_slope_ci_high": float("nan"),
        "adaptive_intercept_ci_low": float("nan"),
        "adaptive_intercept_ci_high": float("nan"),
        "adaptive_bins_used": 0,
        "adaptive_a_max": float("nan"),
    }
    if positive.size == 0 or negative_magnitude.size == 0:
        return [], empty

    overlap = min(float(np.max(positive)), float(np.max(negative_magnitude)))
    a_max = min(
        float(np.quantile(positive, range_quantile)),
        float(np.quantile(negative_magnitude, range_quantile)),
        overlap,
    )
    if not math.isfinite(a_max) or a_max <= 0.0:
        return [], empty
    positive_in_range = positive[positive <= a_max]
    negative_in_range = negative_magnitude[negative_magnitude <= a_max]
    minority_count = min(positive_in_range.size, negative_in_range.size)
    raw_count_target = max(1, int(math.ceil(min_effective_count * iat)))
    bin_count = min(max_bins, minority_count // raw_count_target)
    if bin_count < 3:
        return [], empty

    edges = np.linspace(0.0, a_max, bin_count + 1)

    plus_by_stream = np.zeros((values.shape[0], edges.size - 1), dtype=int)
    minus_by_stream = np.zeros_like(plus_by_stream)
    for stream in range(values.shape[0]):
        plus_by_stream[stream], _ = np.histogram(values[stream], bins=edges)
        minus_by_stream[stream], _ = np.histogram(-values[stream], bins=edges)
    plus = plus_by_stream.sum(axis=0)
    minus = minus_by_stream.sum(axis=0)
    effective_plus = plus / iat
    effective_minus = minus / iat
    fit_mask = (
        (effective_plus >= min_effective_count)
        & (effective_minus >= min_effective_count)
        & (plus > 0)
        & (minus > 0)
    )
    centers = 0.5 * (edges[:-1] + edges[1:])
    raw_symmetry = np.full(centers.size, np.nan)
    nonzero = (plus > 0) & (minus > 0)
    raw_symmetry[nonzero] = np.log(plus[nonzero] / minus[nonzero]) / tau
    variance = np.full(centers.size, np.nan)
    variance[nonzero] = iat * (
        1.0 / plus[nonzero] + 1.0 / minus[nonzero]
    ) / (tau * tau)

    rows = []
    for index in range(centers.size):
        rows.append(
            {
                "a_low": float(edges[index]),
                "a_high": float(edges[index + 1]),
                "a_center": float(centers[index]),
                "plus_count": int(plus[index]),
                "minus_count": int(minus[index]),
                "effective_plus_count": float(effective_plus[index]),
                "effective_minus_count": float(effective_minus[index]),
                "symmetry_raw": float(raw_symmetry[index]),
                "fit_used": int(fit_mask[index]),
            }
        )
    if np.count_nonzero(fit_mask) < 3:
        return rows, {
            **empty,
            "adaptive_bins_used": int(np.count_nonzero(fit_mask)),
            "adaptive_a_max": a_max,
        }

    slope, intercept, r_squared = weighted_fit(
        centers[fit_mask], raw_symmetry[fit_mask], variance[fit_mask]
    )
    boot_slopes: list[float] = []
    boot_intercepts: list[float] = []
    stream_count = values.shape[0]
    for _ in range(bootstrap):
        selected = rng.integers(0, stream_count, size=stream_count)
        boot_plus = plus_by_stream[selected].sum(axis=0)
        boot_minus = minus_by_stream[selected].sum(axis=0)
        valid = fit_mask & (boot_plus > 0) & (boot_minus > 0)
        if np.count_nonzero(valid) < 3:
            continue
        response = np.log(boot_plus[valid] / boot_minus[valid]) / tau
        boot_variance = iat * (
            1.0 / boot_plus[valid] + 1.0 / boot_minus[valid]
        ) / (tau * tau)
        boot_slope, boot_intercept, _ = weighted_fit(
            centers[valid], response, boot_variance
        )
        boot_slopes.append(boot_slope)
        boot_intercepts.append(boot_intercept)
    slope_low, slope_high = interval(boot_slopes)
    intercept_low, intercept_high = interval(boot_intercepts)
    return rows, {
        "adaptive_slope": slope,
        "adaptive_intercept": intercept,
        "adaptive_r_squared": r_squared,
        "adaptive_slope_ci_low": slope_low,
        "adaptive_slope_ci_high": slope_high,
        "adaptive_intercept_ci_low": intercept_low,
        "adaptive_intercept_ci_high": intercept_high,
        "adaptive_bins_used": int(np.count_nonzero(fit_mask)),
        "adaptive_a_max": a_max,
    }


def plot_slopes(
    output_dir: Path,
    n: int,
    summaries: list[dict[str, object]],
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0))
    for axis, observable, target, title in [
        (axes[0], "entropy_rate", 1.0, "medium entropy"),
        (axes[1], "heat_current", 0.4, "bath heat"),
        (axes[2], "action_current", None, "action current"),
    ]:
        selected = sorted(
            (
                row
                for row in summaries
                if row["n"] == n
                and row["observable"] == observable
                and math.isfinite(float(row["adaptive_slope"]))
            ),
            key=lambda row: float(row["tau"]),
        )
        if selected:
            tau = np.asarray([float(row["tau"]) for row in selected])
            slope = np.asarray([float(row["adaptive_slope"]) for row in selected])
            low = np.asarray([float(row["adaptive_slope_ci_low"]) for row in selected])
            high = np.asarray([float(row["adaptive_slope_ci_high"]) for row in selected])
            error = np.vstack([slope - low, high - slope])
            error[~np.isfinite(error)] = 0.0
            axis.errorbar(1.0 / tau, slope, yerr=error, fmt="o-", capsize=3)
        if target is not None:
            axis.axhline(target, color="k", linestyle="--", label=f"reference {target:g}")
            axis.legend(fontsize=8)
        axis.set_title(title)
        axis.set_xlabel(r"$1/t$")
        axis.set_ylabel("adaptive-bin symmetry slope")
        axis.set_xlim(0.0, 0.055)
        axis.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / f"adaptive_symmetry_slopes_n{n}.png", dpi=220)
    fig.savefig(output_dir / f"adaptive_symmetry_slopes_n{n}.pdf")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    taus = [float(item) for item in args.taus.split(",") if item]
    rng = np.random.default_rng(args.seed)
    all_bins: list[dict[str, object]] = []
    all_summaries: list[dict[str, object]] = []
    chain_lengths = []

    for path in args.blocks:
        metadata, raw = load_run(path)
        n = int(metadata["n"])
        chain_lengths.append(n)
        base_tau = float(metadata["block_time"])
        for tau in taus:
            factor_float = tau / base_tau
            factor = int(round(factor_float))
            if factor <= 0 or not math.isclose(factor_float, factor, abs_tol=1.0e-12):
                continue
            if factor > raw.shape[1]:
                continue
            aggregated = aggregate_blocks(raw, factor, base_tau)
            for observable in ["entropy_rate", "heat_current", "action_current"]:
                values = aggregated[observable]
                iat = integrated_autocorrelation_time(values, args.max_acf_lag)
                rows, fit = adaptive_symmetry(
                    values,
                    tau,
                    args.max_bins,
                    args.min_effective_count,
                    args.range_quantile,
                    iat,
                    args.bootstrap,
                    rng,
                )
                for row in rows:
                    all_bins.append(
                        {
                            "source": str(path),
                            "n": n,
                            "tau": tau,
                            "observable": observable,
                            **row,
                        }
                    )
                all_summaries.append(
                    {
                        "source": str(path),
                        "n": n,
                        "tau": tau,
                        "observable": observable,
                        "n_samples": values.size,
                        "negative_count": int(np.count_nonzero(values < 0.0)),
                        "autocorrelation_time": iat,
                        "max_bins": args.max_bins,
                        "min_effective_count": args.min_effective_count,
                        "range_quantile": args.range_quantile,
                        **fit,
                    }
                )

    write_rows(args.output_dir / "adaptive_symmetry_bins.csv", all_bins)
    write_rows(args.output_dir / "adaptive_symmetry_summary.csv", all_summaries)
    for n in chain_lengths:
        plot_slopes(args.output_dir, n, all_summaries)
    print(f"wrote adaptive-bin robustness analysis to {args.output_dir}")


if __name__ == "__main__":
    main()
