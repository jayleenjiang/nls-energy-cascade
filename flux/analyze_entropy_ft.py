#!/usr/bin/env python3
"""Analyze joint bath-heat, entropy-production, and action-current blocks.

The primary Gallavotti--Cohen diagnostic is

    R_t(a) = [log p_t(a) - log p_t(-a)] / t,

for the medium entropy-production rate a = Sigma_t^m/t.  Symmetric histogram
bins are used.  The plus-four estimate is written for visualization only;
only bins with nonzero raw counts on both sides and sufficient effective
counts enter the weighted fit.

Input files are ``*_blocks.csv`` outputs from NLS_entropy_ft.cpp.  Metadata are
read from the matching ``*_summary.csv`` files.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm


BLOCK_COLUMNS = {
    "stream_id": 0,
    "block_id": 1,
    "q_left": 2,
    "q_right": 3,
    "delta_energy": 4,
    "entropy_medium": 5,
    "entropy_rate": 6,
    "action_current": 7,
    "energy_balance_error": 8,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("blocks", nargs="+", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--taus",
        default="20,40,60,80,100,120,140,160,180,200",
        help="comma-separated aggregation windows",
    )
    parser.add_argument("--symmetric-bins", type=int, default=40)
    parser.add_argument("--min-effective-count", type=float, default=20.0)
    parser.add_argument("--bootstrap", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260826)
    parser.add_argument("--max-acf-lag", type=int, default=400)
    parser.add_argument("--target-tail-count", type=int, default=200)
    return parser.parse_args()


def read_single_row(path: Path) -> dict[str, str]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 1:
        raise ValueError(f"expected one metadata row in {path}, found {len(rows)}")
    return rows[0]


def load_run(path: Path) -> tuple[dict[str, str], np.ndarray]:
    if not path.name.endswith("_blocks.csv"):
        raise ValueError(f"input must end in _blocks.csv: {path}")
    summary = path.with_name(path.name.replace("_blocks.csv", "_summary.csv"))
    if not summary.exists():
        raise FileNotFoundError(f"missing matching summary: {summary}")
    metadata = read_single_row(summary)
    values = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    if values.shape[1] != len(BLOCK_COLUMNS):
        raise ValueError(
            f"{path}: expected {len(BLOCK_COLUMNS)} columns, got {values.shape[1]}"
        )
    order = np.lexsort(
        (values[:, BLOCK_COLUMNS["block_id"]], values[:, BLOCK_COLUMNS["stream_id"]])
    )
    values = values[order]
    stream_ids = values[:, BLOCK_COLUMNS["stream_id"]].astype(int)
    unique, counts = np.unique(stream_ids, return_counts=True)
    if unique.size == 0 or not np.all(counts == counts[0]):
        raise ValueError(f"{path}: streams do not contain equal block counts")
    blocks_per_stream = int(counts[0])
    reshaped = values.reshape(unique.size, blocks_per_stream, values.shape[1])
    expected_stream = np.arange(unique.size)
    if not np.array_equal(unique, expected_stream):
        raise ValueError(f"{path}: stream ids must be contiguous from zero")
    expected_blocks = np.arange(blocks_per_stream)
    if not np.all(
        reshaped[:, :, BLOCK_COLUMNS["block_id"]].astype(int) == expected_blocks
    ):
        raise ValueError(f"{path}: block ids must be contiguous within each stream")
    return metadata, reshaped


def aggregate_blocks(data: np.ndarray, factor: int, base_tau: float) -> dict[str, np.ndarray]:
    usable = (data.shape[1] // factor) * factor
    if usable == 0:
        raise ValueError(f"not enough blocks to aggregate by factor {factor}")
    grouped = data[:, :usable, :].reshape(data.shape[0], usable // factor, factor, -1)
    tau = base_tau * factor
    q_left = grouped[:, :, :, BLOCK_COLUMNS["q_left"]].sum(axis=2)
    q_right = grouped[:, :, :, BLOCK_COLUMNS["q_right"]].sum(axis=2)
    delta_energy = grouped[:, :, :, BLOCK_COLUMNS["delta_energy"]].sum(axis=2)
    entropy = grouped[:, :, :, BLOCK_COLUMNS["entropy_medium"]].sum(axis=2)
    action = grouped[:, :, :, BLOCK_COLUMNS["action_current"]].mean(axis=2)
    balance = grouped[:, :, :, BLOCK_COLUMNS["energy_balance_error"]].sum(axis=2)
    return {
        "tau": np.asarray(tau),
        "q_left": q_left,
        "q_right": q_right,
        "heat_current": (q_left - q_right) / (2.0 * tau),
        "delta_energy_rate": delta_energy / tau,
        "entropy_rate": entropy / tau,
        "action_current": action,
        "balance_rate": balance / tau,
    }


def integrated_autocorrelation_time(values: np.ndarray, max_lag: int) -> float:
    """Positive-sequence IAT pooled over equal-length independent streams."""
    if values.ndim != 2:
        raise ValueError("values must have shape (streams, blocks)")
    if values.shape[1] < 3:
        return 1.0
    centered = values - np.mean(values)
    denominator = float(np.sum(centered * centered))
    if not np.isfinite(denominator) or denominator <= 0.0:
        return 1.0
    upper = min(max_lag, values.shape[1] - 1)
    positive_sum = 0.0
    for lag in range(1, upper + 1):
        covariance = float(np.sum(centered[:, :-lag] * centered[:, lag:]))
        scale = values.shape[1] / (values.shape[1] - lag)
        rho = scale * covariance / denominator
        if not np.isfinite(rho) or rho <= 0.0:
            break
        positive_sum += rho
    return max(1.0, 1.0 + 2.0 * positive_sum)


def weighted_line_fit(x: np.ndarray, y: np.ndarray, variance: np.ndarray) -> tuple[float, float, float]:
    weights = 1.0 / np.sqrt(np.maximum(variance, np.finfo(float).tiny))
    design = np.column_stack([x, np.ones_like(x)])
    coefficients, *_ = np.linalg.lstsq(design * weights[:, None], y * weights, rcond=None)
    slope, intercept = coefficients
    fitted = slope * x + intercept
    weighted_mean = np.average(y, weights=weights * weights)
    residual = np.sum((weights * (y - fitted)) ** 2)
    total = np.sum((weights * (y - weighted_mean)) ** 2)
    r_squared = float(1.0 - residual / total) if total > 0.0 else float("nan")
    return float(slope), float(intercept), r_squared


def percentile_interval(values: list[float]) -> tuple[float, float]:
    finite = np.asarray([value for value in values if np.isfinite(value)])
    if finite.size < 20:
        return float("nan"), float("nan")
    low, high = np.quantile(finite, [0.025, 0.975])
    return float(low), float(high)


def symmetry_analysis(
    entropy_rate: np.ndarray,
    tau: float,
    symmetric_bins: int,
    min_effective_count: float,
    autocorrelation_time: float,
    bootstrap: int,
    rng: np.random.Generator,
) -> tuple[list[dict[str, float]], dict[str, float]]:
    flattened = entropy_rate.ravel()
    negative = flattened[flattened < 0.0]
    positive = flattened[flattened > 0.0]
    if negative.size == 0 or positive.size == 0:
        return [], {
            "ft_slope": float("nan"),
            "ft_intercept": float("nan"),
            "ft_r_squared": float("nan"),
            "ft_slope_ci_low": float("nan"),
            "ft_slope_ci_high": float("nan"),
            "ft_intercept_ci_low": float("nan"),
            "ft_intercept_ci_high": float("nan"),
            "ft_bins_used": 0,
        }

    q01, q99 = np.quantile(flattened, [0.01, 0.99])
    observed_overlap = min(float(np.max(positive)), float(-np.min(negative)))
    central_overlap = min(float(q99), float(-q01)) if q01 < 0.0 else 0.0
    a_max = central_overlap if central_overlap > 0.0 else observed_overlap
    a_max = min(a_max, observed_overlap)
    if not np.isfinite(a_max) or a_max <= 0.0:
        return [], {
            "ft_slope": float("nan"),
            "ft_intercept": float("nan"),
            "ft_r_squared": float("nan"),
            "ft_slope_ci_low": float("nan"),
            "ft_slope_ci_high": float("nan"),
            "ft_intercept_ci_low": float("nan"),
            "ft_intercept_ci_high": float("nan"),
            "ft_bins_used": 0,
        }

    edges = np.linspace(0.0, a_max, symmetric_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    plus_counts_by_stream = np.zeros((entropy_rate.shape[0], symmetric_bins), dtype=int)
    minus_counts_by_stream = np.zeros_like(plus_counts_by_stream)
    for stream in range(entropy_rate.shape[0]):
        plus_counts_by_stream[stream], _ = np.histogram(entropy_rate[stream], bins=edges)
        minus_counts_by_stream[stream], _ = np.histogram(-entropy_rate[stream], bins=edges)
    plus_counts = plus_counts_by_stream.sum(axis=0)
    minus_counts = minus_counts_by_stream.sum(axis=0)
    sample_count = flattened.size
    effective_plus = plus_counts / autocorrelation_time
    effective_minus = minus_counts / autocorrelation_time
    fit_mask = (
        (plus_counts > 0)
        & (minus_counts > 0)
        & (effective_plus >= min_effective_count)
        & (effective_minus >= min_effective_count)
    )

    plus_four = (plus_counts + 2.0) / (sample_count + 4.0)
    minus_four = (minus_counts + 2.0) / (sample_count + 4.0)
    symmetry_plot = (np.log(plus_four) - np.log(minus_four)) / tau
    symmetry_raw = np.full(symmetric_bins, np.nan)
    nonzero = (plus_counts > 0) & (minus_counts > 0)
    symmetry_raw[nonzero] = (
        np.log(plus_counts[nonzero] / sample_count)
        - np.log(minus_counts[nonzero] / sample_count)
    ) / tau
    variance = np.full(symmetric_bins, np.nan)
    variance[nonzero] = autocorrelation_time * (
        1.0 / plus_counts[nonzero] + 1.0 / minus_counts[nonzero]
    ) / (tau * tau)

    rows: list[dict[str, float]] = []
    for index in range(symmetric_bins):
        rows.append(
            {
                "a_low": float(edges[index]),
                "a_high": float(edges[index + 1]),
                "a_center": float(centers[index]),
                "plus_count": int(plus_counts[index]),
                "minus_count": int(minus_counts[index]),
                "effective_plus_count": float(effective_plus[index]),
                "effective_minus_count": float(effective_minus[index]),
                "symmetry_plus_four": float(symmetry_plot[index]),
                "symmetry_raw": float(symmetry_raw[index]),
                "fit_used": int(fit_mask[index]),
            }
        )

    if np.count_nonzero(fit_mask) < 3:
        return rows, {
            "ft_slope": float("nan"),
            "ft_intercept": float("nan"),
            "ft_r_squared": float("nan"),
            "ft_slope_ci_low": float("nan"),
            "ft_slope_ci_high": float("nan"),
            "ft_intercept_ci_low": float("nan"),
            "ft_intercept_ci_high": float("nan"),
            "ft_bins_used": int(np.count_nonzero(fit_mask)),
        }

    slope, intercept, r_squared = weighted_line_fit(
        centers[fit_mask], symmetry_raw[fit_mask], variance[fit_mask]
    )
    slope_bootstrap: list[float] = []
    intercept_bootstrap: list[float] = []
    stream_count = entropy_rate.shape[0]
    for _ in range(bootstrap):
        selected = rng.integers(0, stream_count, size=stream_count)
        plus = plus_counts_by_stream[selected].sum(axis=0)
        minus = minus_counts_by_stream[selected].sum(axis=0)
        valid = fit_mask & (plus > 0) & (minus > 0)
        if np.count_nonzero(valid) < 3:
            continue
        response = (np.log(plus[valid]) - np.log(minus[valid])) / tau
        boot_variance = autocorrelation_time * (
            1.0 / plus[valid] + 1.0 / minus[valid]
        ) / (tau * tau)
        boot_slope, boot_intercept, _ = weighted_line_fit(
            centers[valid], response, boot_variance
        )
        slope_bootstrap.append(boot_slope)
        intercept_bootstrap.append(boot_intercept)
    slope_low, slope_high = percentile_interval(slope_bootstrap)
    intercept_low, intercept_high = percentile_interval(intercept_bootstrap)
    return rows, {
        "ft_slope": slope,
        "ft_intercept": intercept,
        "ft_r_squared": r_squared,
        "ft_slope_ci_low": slope_low,
        "ft_slope_ci_high": slope_high,
        "ft_intercept_ci_low": intercept_low,
        "ft_intercept_ci_high": intercept_high,
        "ft_bins_used": int(np.count_nonzero(fit_mask)),
    }


def coupling_metrics(action: np.ndarray, heat: np.ndarray) -> dict[str, float]:
    x = action.ravel()
    y = heat.ravel()
    correlation = float(np.corrcoef(x, y)[0, 1])
    design = np.column_stack([x, np.ones_like(x)])
    slope, intercept = np.linalg.lstsq(design, y, rcond=None)[0]
    residual = y - slope * x - intercept
    variance_y = float(np.var(y, ddof=1))
    residual_variance = float(np.var(residual, ddof=1))
    return {
        "heat_on_action_slope": float(slope),
        "heat_on_action_intercept": float(intercept),
        "pearson_correlation": correlation,
        "residual_variance": residual_variance,
        "heat_variance": variance_y,
        "residual_variance_fraction": residual_variance / variance_y
        if variance_y > 0.0
        else float("nan"),
    }


def stationarity_metrics(values: np.ndarray) -> dict[str, float]:
    quarter = max(1, values.shape[1] // 4)
    first = np.mean(values[:, :quarter], axis=1)
    last = np.mean(values[:, -quarter:], axis=1)
    difference = last - first
    se = float(np.std(difference, ddof=1) / math.sqrt(difference.size))
    return {
        "first_quarter_mean": float(np.mean(first)),
        "last_quarter_mean": float(np.mean(last)),
        "last_minus_first": float(np.mean(difference)),
        "paired_se": se,
        "paired_z": float(np.mean(difference) / se) if se > 0.0 else float("nan"),
    }


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_run(
    output_dir: Path,
    n: int,
    aggregated_by_tau: dict[float, dict[str, np.ndarray]],
    ft_rows: list[dict[str, object]],
    ft_summary: list[dict[str, object]],
    rng: np.random.Generator,
) -> None:
    selected_taus = [tau for tau in [20, 40, 80, 120, 160, 200] if tau in aggregated_by_tau]
    if not selected_taus:
        selected_taus = sorted(aggregated_by_tau)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    for tau in selected_taus:
        values = aggregated_by_tau[tau]["entropy_rate"].ravel()
        q01, q99 = np.quantile(values, [0.01, 0.99])
        axes[0].hist(
            values,
            bins=80,
            range=(q01, q99),
            density=True,
            histtype="step",
            linewidth=1.2,
            label=rf"$t={tau:g}$",
        )
    axes[0].set_xlabel(r"entropy-production rate $a=\Sigma_t^{\rm m}/t$")
    axes[0].set_ylabel("PDF")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.2)

    for tau in selected_taus:
        rows = [row for row in ft_rows if row["n"] == n and row["tau"] == tau]
        if not rows:
            continue
        x = np.asarray([row["a_center"] for row in rows])
        y = np.asarray([row["symmetry_plus_four"] for row in rows])
        used = np.asarray([row["fit_used"] for row in rows], dtype=bool)
        axes[1].plot(x, y, color="0.75", linewidth=0.8)
        if np.any(used):
            axes[1].plot(x[used], y[used], marker="o", markersize=3, label=rf"$t={tau:g}$")
    limits = axes[1].get_xlim()
    upper = max(0.0, limits[1])
    axes[1].plot([0.0, upper], [0.0, upper], "k--", linewidth=1.2, label=r"$R_t(a)=a$")
    axes[1].set_xlabel(r"$a$")
    axes[1].set_ylabel(r"$R_t(a)=t^{-1}\log[p_t(a)/p_t(-a)]$")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / f"entropy_pdf_symmetry_n{n}.png", dpi=220)
    fig.savefig(output_dir / f"entropy_pdf_symmetry_n{n}.pdf")
    plt.close(fig)

    rows = sorted((row for row in ft_summary if row["n"] == n), key=lambda row: row["tau"])
    if rows:
        tau = np.asarray([row["tau"] for row in rows], dtype=float)
        slope = np.asarray([row["ft_slope"] for row in rows], dtype=float)
        low = np.asarray([row["ft_slope_ci_low"] for row in rows], dtype=float)
        high = np.asarray([row["ft_slope_ci_high"] for row in rows], dtype=float)
        finite = np.isfinite(slope)
        fig, ax = plt.subplots(figsize=(5.8, 4.2))
        if np.any(finite):
            error = np.vstack([slope[finite] - low[finite], high[finite] - slope[finite]])
            error[~np.isfinite(error)] = 0.0
            ax.errorbar(1.0 / tau[finite], slope[finite], yerr=error, fmt="o-", capsize=3)
        ax.axhline(1.0, color="k", linestyle="--", label="FT slope 1")
        ax.set_xlabel(r"$1/t$")
        ax.set_ylabel(r"fitted symmetry slope $m_t$")
        ax.grid(alpha=0.2)
        ax.legend()
        fig.tight_layout()
        fig.savefig(output_dir / f"ft_slope_vs_time_n{n}.png", dpi=220)
        fig.savefig(output_dir / f"ft_slope_vs_time_n{n}.pdf")
        plt.close(fig)

    first_tau = min(aggregated_by_tau)
    last_tau = max(aggregated_by_tau)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    for ax, tau in zip(axes, [first_tau, last_tau]):
        action = aggregated_by_tau[tau]["action_current"].ravel()
        heat = aggregated_by_tau[tau]["heat_current"].ravel()
        if action.size > 8000:
            indices = rng.choice(action.size, size=8000, replace=False)
            action = action[indices]
            heat = heat[indices]
        ax.hexbin(action, heat, gridsize=45, mincnt=1, bins="log", cmap="viridis")
        ax.set_title(rf"$t={tau:g}$")
        ax.set_xlabel(r"action current $\overline{J}_M$")
        ax.set_ylabel(r"bath heat current $\overline{J}_E$")
        ax.grid(alpha=0.15)
    fig.tight_layout()
    fig.savefig(output_dir / f"heat_action_coupling_n{n}.png", dpi=220)
    fig.savefig(output_dir / f"heat_action_coupling_n{n}.pdf")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    requested_taus = [float(value) for value in args.taus.split(",") if value.strip()]
    rng = np.random.default_rng(args.seed)

    ft_summary_rows: list[dict[str, object]] = []
    ft_bin_rows: list[dict[str, object]] = []
    coupling_rows: list[dict[str, object]] = []
    stationarity_rows: list[dict[str, object]] = []
    recommendation_rows: list[dict[str, object]] = []

    for path in args.blocks:
        metadata, raw = load_run(path)
        n = int(metadata["n"])
        base_tau = float(metadata["block_time"])
        aggregated_by_tau: dict[float, dict[str, np.ndarray]] = {}
        run_ft_rows: list[dict[str, object]] = []
        run_ft_summary: list[dict[str, object]] = []

        for tau in requested_taus:
            factor_float = tau / base_tau
            factor = int(round(factor_float))
            if factor <= 0 or not math.isclose(factor_float, factor, rel_tol=0.0, abs_tol=1e-12):
                continue
            if factor > raw.shape[1]:
                continue
            aggregated = aggregate_blocks(raw, factor, base_tau)
            aggregated_by_tau[tau] = aggregated
            entropy = aggregated["entropy_rate"]
            sample_count = entropy.size
            iat = integrated_autocorrelation_time(entropy, args.max_acf_lag)
            effective_sample_size = sample_count / iat
            flattened = entropy.ravel()
            mean = float(np.mean(flattened))
            std = float(np.std(flattened, ddof=1))
            negative_count = int(np.count_nonzero(flattened < 0.0))
            negative_probability = negative_count / sample_count
            normal_symmetry_slope = 2.0 * mean / (tau * std * std) if std > 0.0 else float("nan")

            bins, fit = symmetry_analysis(
                entropy,
                tau,
                args.symmetric_bins,
                args.min_effective_count,
                iat,
                args.bootstrap,
                rng,
            )
            for row in bins:
                full = {"source": str(path), "n": n, "tau": tau, **row}
                ft_bin_rows.append(full)
                run_ft_rows.append(full)
            summary = {
                "source": str(path),
                "n": n,
                "tau": tau,
                "n_samples": sample_count,
                "autocorrelation_time": iat,
                "effective_sample_size": effective_sample_size,
                "mean_entropy_rate": mean,
                "std_entropy_rate": std,
                "negative_count": negative_count,
                "negative_probability": negative_probability,
                "normal_symmetry_slope": normal_symmetry_slope,
                "mean_balance_error_rate": float(np.mean(aggregated["balance_rate"])),
                "rms_balance_error_rate": float(
                    np.sqrt(np.mean(aggregated["balance_rate"] ** 2))
                ),
                **fit,
            }
            ft_summary_rows.append(summary)
            run_ft_summary.append(summary)

            coupling = coupling_metrics(
                aggregated["action_current"], aggregated["heat_current"]
            )
            coupling_rows.append(
                {
                    "source": str(path),
                    "n": n,
                    "tau": tau,
                    "n_samples": sample_count,
                    "mean_action_current": float(np.mean(aggregated["action_current"])),
                    "mean_heat_current": float(np.mean(aggregated["heat_current"])),
                    **coupling,
                }
            )

            for observable in ["entropy_rate", "action_current", "heat_current"]:
                stationarity_rows.append(
                    {
                        "source": str(path),
                        "n": n,
                        "tau": tau,
                        "observable": observable,
                        **stationarity_metrics(aggregated[observable]),
                    }
                )

            correlation_penalty = iat
            if negative_probability > 0.0:
                required = math.ceil(
                    args.target_tail_count * correlation_penalty / negative_probability
                )
            else:
                required = float("inf")
            recommendation_rows.append(
                {
                    "source": str(path),
                    "n": n,
                    "tau": tau,
                    "negative_probability": negative_probability,
                    "autocorrelation_time": iat,
                    "target_effective_negative_count": args.target_tail_count,
                    "recommended_raw_sample_count": required,
                }
            )

        plot_run(
            args.output_dir,
            n,
            aggregated_by_tau,
            run_ft_rows,
            run_ft_summary,
            rng,
        )

    write_rows(args.output_dir / "ft_summary.csv", ft_summary_rows)
    write_rows(args.output_dir / "ft_symmetric_bins.csv", ft_bin_rows)
    write_rows(args.output_dir / "coupling_summary.csv", coupling_rows)
    write_rows(args.output_dir / "stationarity_summary.csv", stationarity_rows)
    write_rows(
        args.output_dir / "sample_size_recommendations.csv", recommendation_rows
    )

    print(f"wrote analysis to {args.output_dir}")
    for row in ft_summary_rows:
        print(
            "n={n} t={tau:g} N={n_samples} neg={negative_count} "
            "m={ft_slope:.6g} CI=[{ft_slope_ci_low:.6g},"
            "{ft_slope_ci_high:.6g}] bins={ft_bins_used}".format(**row)
        )


if __name__ == "__main__":
    main()
