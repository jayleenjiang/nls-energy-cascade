#!/usr/bin/env python3
"""Cross-fitted n=2 NESS-density pilot for the total-entropy endpoint term.

This is deliberately a validation-first pilot.  A smoothing scale is selected
only from an equal-temperature data set, where the exact stationary density
ratio is known.  The selected estimator is then applied unchanged to an
independent driven data set.  Failure of the equilibrium density-ratio check
blocks any total-entropy interpretation of the driven result.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter
from scipy.special import logsumexp


@dataclass
class Blocks:
    stream: np.ndarray
    block: np.ndarray
    start: np.ndarray
    end: np.ndarray
    energy_start: np.ndarray
    energy_end: np.ndarray
    entropy_medium: np.ndarray
    balance_error: np.ndarray


def load_blocks(path: Path) -> Blocks:
    values = np.genfromtxt(path, delimiter=",", names=True, dtype=None,
                           encoding=None)
    if values.size == 0:
        raise ValueError(f"empty endpoint file: {path}")
    return Blocks(
        stream=np.asarray(values["stream_id"], dtype=int),
        block=np.asarray(values["block_id"], dtype=int),
        start=np.column_stack([
            values["log_action_1_start"], values["log_action_2_start"],
            values["theta_start"],
        ]).astype(float),
        end=np.column_stack([
            values["log_action_1_end"], values["log_action_2_end"],
            values["theta_end"],
        ]).astype(float),
        energy_start=np.asarray(values["energy_start"], dtype=float),
        energy_end=np.asarray(values["energy_end"], dtype=float),
        entropy_medium=np.asarray(values["entropy_medium"], dtype=float),
        balance_error=np.asarray(values["energy_balance_error"], dtype=float),
    )


def wrap_angle(theta: np.ndarray) -> np.ndarray:
    return (theta + np.pi) % (2.0 * np.pi) - np.pi


def density_evaluator(train: np.ndarray, bins: tuple[int, int, int],
                      sigma: float):
    train = train.copy()
    train[:, 2] = wrap_angle(train[:, 2])
    mean = train[:, :2].mean(axis=0)
    scale = train[:, :2].std(axis=0, ddof=1)
    standardized = train.copy()
    standardized[:, :2] = (standardized[:, :2] - mean) / scale

    low = np.quantile(standardized[:, :2], 0.0001, axis=0) - 0.5
    high = np.quantile(standardized[:, :2], 0.9999, axis=0) + 0.5
    edges = [
        np.linspace(low[0], high[0], bins[0] + 1),
        np.linspace(low[1], high[1], bins[1] + 1),
        np.linspace(-np.pi, np.pi, bins[2] + 1),
    ]
    counts, _ = np.histogramdd(standardized, bins=edges)
    smooth = gaussian_filter(counts, sigma=(sigma, sigma, sigma),
                             mode=("nearest", "nearest", "wrap"))
    centers = [0.5 * (edge[1:] + edge[:-1]) for edge in edges]
    theta_centers = np.concatenate([
        centers[2][-1:] - 2.0 * np.pi,
        centers[2],
        centers[2][:1] + 2.0 * np.pi,
    ])
    extended = np.concatenate([smooth[:, :, -1:], smooth,
                               smooth[:, :, :1]], axis=2)
    floor = max(float(smooth.sum()) * 1.0e-15, 1.0e-300)
    log_extended = np.log(np.maximum(extended, floor))
    log_interpolator = RegularGridInterpolator(
        (centers[0], centers[1], theta_centers), log_extended,
        bounds_error=False, fill_value=np.nan,
    )
    count_interpolator = RegularGridInterpolator(
        (centers[0], centers[1], theta_centers), extended,
        bounds_error=False, fill_value=0.0,
    )

    def evaluate(points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        transformed = points.copy()
        transformed[:, :2] = (transformed[:, :2] - mean) / scale
        transformed[:, 2] = wrap_angle(transformed[:, 2])
        log_density = np.asarray(log_interpolator(transformed), dtype=float)
        local_count = np.asarray(count_interpolator(transformed), dtype=float)
        valid = np.isfinite(log_density) & (local_count >= 0.05)
        log_density[~valid] = np.nan
        return log_density, valid

    return evaluate


def crossfit_log_ratio(blocks: Blocks, bins: tuple[int, int, int],
                       sigma: float) -> tuple[np.ndarray, np.ndarray]:
    ratio = np.full(blocks.stream.size, np.nan)
    valid = np.zeros(blocks.stream.size, dtype=bool)
    for fold in (0, 1):
        test = blocks.stream % 2 == fold
        train = ~test
        # Consecutive blocks share an endpoint, so using both starts and ends
        # would duplicate nearly every training state.  Use starts only.
        train_points = blocks.start[train]
        evaluate = density_evaluator(train_points, bins, sigma)
        log_start, valid_start = evaluate(blocks.start[test])
        log_end, valid_end = evaluate(blocks.end[test])
        # q is the density with respect to du1 du2 dtheta.  Cartesian volume
        # contributes I1 I2 = exp(u1+u2), so log rho = log q-u1-u2+const.
        result = (
            log_start - log_end
            - blocks.start[test, 0] - blocks.start[test, 1]
            + blocks.end[test, 0] + blocks.end[test, 1]
        )
        fold_valid = valid_start & valid_end & np.isfinite(result)
        indices = np.flatnonzero(test)
        ratio[indices[fold_valid]] = result[fold_valid]
        valid[indices[fold_valid]] = True
    return ratio, valid


def subset(blocks: Blocks, keep: np.ndarray) -> Blocks:
    return Blocks(
        stream=blocks.stream[keep],
        block=blocks.block[keep],
        start=blocks.start[keep],
        end=blocks.end[keep],
        energy_start=blocks.energy_start[keep],
        energy_end=blocks.energy_end[keep],
        entropy_medium=blocks.entropy_medium[keep],
        balance_error=blocks.balance_error[keep],
    )


def trained_log_ratio(train_points: np.ndarray, test: Blocks,
                      bins: tuple[int, int, int],
                      sigma: float) -> tuple[np.ndarray, np.ndarray]:
    evaluate = density_evaluator(train_points, bins, sigma)
    log_start, valid_start = evaluate(test.start)
    log_end, valid_end = evaluate(test.end)
    ratio = (
        log_start - log_end
        - test.start[:, 0] - test.start[:, 1]
        + test.end[:, 0] + test.end[:, 1]
    )
    valid = valid_start & valid_end & np.isfinite(ratio)
    ratio[~valid] = np.nan
    return ratio, valid


def stable_log_mean_exp(values: np.ndarray) -> float:
    return float(logsumexp(values) - math.log(values.size))


def weight_ess(log_weights: np.ndarray) -> float:
    shifted = log_weights - float(np.max(log_weights))
    weights = np.exp(shifted)
    return float(weights.sum() ** 2 / np.square(weights).sum())


def stream_bootstrap_log_ift(total_entropy: np.ndarray, stream: np.ndarray,
                             seed: int = 20260827,
                             replicates: int = 1000) -> tuple[float, float]:
    unique = np.unique(stream)
    by_stream = {value: total_entropy[stream == value] for value in unique}
    rng = np.random.default_rng(seed)
    estimates = np.empty(replicates)
    for replicate in range(replicates):
        chosen = rng.choice(unique, size=unique.size, replace=True)
        values = np.concatenate([by_stream[value] for value in chosen])
        estimates[replicate] = stable_log_mean_exp(-values)
    low, high = np.quantile(estimates, [0.025, 0.975])
    return float(low), float(high)


def entropy_metrics(name: str, blocks: Blocks, log_ratio: np.ndarray,
                    valid: np.ndarray) -> dict[str, object]:
    entropy = blocks.entropy_medium[valid] + log_ratio[valid]
    stream = blocks.stream[valid]
    log_ift = stable_log_mean_exp(-entropy)
    ci_low, ci_high = stream_bootstrap_log_ift(entropy, stream)
    return {
        "dataset": name,
        "samples_total": int(valid.size),
        "samples_supported": int(valid.sum()),
        "support_fraction": float(valid.mean()),
        "mean_total_entropy": float(entropy.mean()),
        "std_total_entropy": float(entropy.std(ddof=1)),
        "log_mean_exp_minus_total_entropy": log_ift,
        "log_ift_stream_bootstrap_low": ci_low,
        "log_ift_stream_bootstrap_high": ci_high,
        "exponential_weight_ess": weight_ess(-entropy),
        "unique_streams": int(np.unique(stream).size),
        "rms_energy_balance_error": float(
            np.sqrt(np.mean(np.square(blocks.balance_error)))
        ),
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def symmetry_diagnostic(total_entropy: np.ndarray, bins: int = 60,
                        min_count: int = 20) -> tuple[list[dict[str, object]],
                                                      dict[str, object]]:
    extent = float(np.quantile(np.abs(total_entropy), 0.995))
    edges = np.linspace(-extent, extent, 2 * bins + 1)
    counts, _ = np.histogram(total_entropy, bins=edges)
    rows: list[dict[str, object]] = []
    for index in range(bins, 2 * bins):
        mirror = 2 * bins - 1 - index
        positive = int(counts[index])
        negative = int(counts[mirror])
        center = float(0.5 * (edges[index] + edges[index + 1]))
        usable = positive >= min_count and negative >= min_count
        log_ratio = math.log(positive / negative) if positive and negative else math.nan
        standard_error = (
            math.sqrt(1.0 / positive + 1.0 / negative)
            if positive and negative else math.nan
        )
        rows.append({
            "entropy_center": center,
            "positive_count": positive,
            "negative_count": negative,
            "log_probability_ratio": log_ratio,
            "poisson_standard_error": standard_error,
            "usable": int(usable),
        })
    usable_rows = [row for row in rows if int(row["usable"]) == 1]
    if len(usable_rows) >= 3:
        x = np.asarray([float(row["entropy_center"]) for row in usable_rows])
        y = np.asarray([float(row["log_probability_ratio"]) for row in usable_rows])
        se = np.asarray([float(row["poisson_standard_error"]) for row in usable_rows])
        design = np.column_stack([x, np.ones_like(x)])
        weighted = design / se[:, None]
        coefficients, *_ = np.linalg.lstsq(weighted, y / se, rcond=None)
        prediction = design @ coefficients
        residual = y - prediction
        dof = max(1, x.size - 2)
        covariance = np.linalg.inv(weighted.T @ weighted) * float(
            np.sum(np.square(residual / se)) / dof
        )
        slope = float(coefficients[0])
        intercept = float(coefficients[1])
        slope_se = float(math.sqrt(covariance[0, 0]))
    else:
        slope = intercept = slope_se = math.nan
    summary = {
        "usable_bins": len(usable_rows),
        "weighted_slope": slope,
        "weighted_slope_se": slope_se,
        "weighted_intercept": intercept,
        "reference_slope": 1.0,
    }
    return rows, summary


def stationarity_rows(name: str, blocks: Blocks) -> list[dict[str, object]]:
    observables = {
        "log_action_1": blocks.start[:, 0],
        "log_action_2": blocks.start[:, 1],
        "cos_theta": np.cos(blocks.start[:, 2]),
        "sin_theta": np.sin(blocks.start[:, 2]),
        "energy": blocks.energy_start,
    }
    result: list[dict[str, object]] = []
    unique = np.unique(blocks.stream)
    blocks_per_stream = int(blocks.block.max()) + 1
    quarter = blocks_per_stream // 4
    if quarter < 1:
        raise ValueError("stationarity check requires at least four blocks")
    first = blocks.block < quarter
    last = blocks.block >= blocks_per_stream - quarter
    for observable, values in observables.items():
        differences: list[float] = []
        first_values: list[float] = []
        last_values: list[float] = []
        for value in unique:
            in_stream = blocks.stream == value
            first_mean = float(np.mean(values[in_stream & first]))
            last_mean = float(np.mean(values[in_stream & last]))
            differences.append(last_mean - first_mean)
            first_values.append(first_mean)
            last_values.append(last_mean)
        stream_differences = np.asarray(differences)
        difference = float(stream_differences.mean())
        se = float(stream_differences.std(ddof=1) /
                   math.sqrt(stream_differences.size))
        result.append({
            "dataset": name,
            "observable": observable,
            "first_quarter_mean": float(np.mean(first_values)),
            "last_quarter_mean": float(np.mean(last_values)),
            "last_minus_first": difference,
            "stream_standard_error": se,
            "z_score": difference / se if se > 0.0 else math.nan,
            "streams": int(unique.size),
        })
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("equilibrium_blocks", type=Path)
    parser.add_argument("driven_blocks", type=Path)
    parser.add_argument("--equilibrium-temperature", type=float, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bins", type=int, nargs=3, default=(24, 24, 48))
    parser.add_argument("--sigma-grid", type=float, nargs="+",
                        default=(0.75, 1.0, 1.5, 2.0))
    args = parser.parse_args()
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)

    equilibrium = load_blocks(args.equilibrium_blocks)
    driven = load_blocks(args.driven_blocks)
    bins = tuple(args.bins)
    validation_rows: list[dict[str, object]] = []
    unique_streams = np.unique(equilibrium.stream)
    if unique_streams.size < 64:
        raise ValueError("equilibrium pilot requires at least 64 streams")
    split = unique_streams.size // 2
    tuning_ids = unique_streams[:split]
    tuning = subset(equilibrium, np.isin(equilibrium.stream, tuning_ids))
    validation = subset(equilibrium, ~np.isin(equilibrium.stream, tuning_ids))
    tuning_target = (
        tuning.energy_end - tuning.energy_start
    ) / args.equilibrium_temperature
    for sigma in args.sigma_grid:
        ratio, valid = crossfit_log_ratio(tuning, bins, sigma)
        error = ratio[valid] - tuning_target[valid]
        slope = float(np.polyfit(tuning_target[valid], ratio[valid], 1)[0])
        correlation = float(
            np.corrcoef(tuning_target[valid], ratio[valid])[0, 1]
        )
        validation_rows.append(
            {
                "stage": "tuning_crossfit",
                "sigma_bins": sigma,
                "samples_supported": int(valid.sum()),
                "support_fraction": float(valid.mean()),
                "density_log_ratio_rmse": float(
                    np.sqrt(np.mean(np.square(error)))
                ),
                "density_log_ratio_mae": float(np.mean(np.abs(error))),
                "density_log_ratio_slope": slope,
                "density_log_ratio_correlation": correlation,
            }
        )
    best = min(validation_rows, key=lambda row: row["density_log_ratio_rmse"])
    sigma = float(best["sigma_bins"])
    equilibrium_ratio, equilibrium_valid = trained_log_ratio(
        tuning.start, validation, bins, sigma
    )
    validation_target = (
        validation.energy_end - validation.energy_start
    ) / args.equilibrium_temperature
    validation_error = (
        equilibrium_ratio[equilibrium_valid]
        - validation_target[equilibrium_valid]
    )
    validation_rows.append(
        {
            "stage": "heldout_validation",
            "sigma_bins": sigma,
            "samples_supported": int(equilibrium_valid.sum()),
            "support_fraction": float(equilibrium_valid.mean()),
            "density_log_ratio_rmse": float(
                np.sqrt(np.mean(np.square(validation_error)))
            ),
            "density_log_ratio_mae": float(np.mean(np.abs(validation_error))),
            "density_log_ratio_slope": float(np.polyfit(
                validation_target[equilibrium_valid],
                equilibrium_ratio[equilibrium_valid], 1,
            )[0]),
            "density_log_ratio_correlation": float(np.corrcoef(
                validation_target[equilibrium_valid],
                equilibrium_ratio[equilibrium_valid],
            )[0, 1]),
        }
    )
    write_csv(output / "equilibrium_density_validation.csv", validation_rows)
    driven_ratio, driven_valid = crossfit_log_ratio(driven, bins, sigma)

    metrics = [
        entropy_metrics("equilibrium_learned", validation,
                        equilibrium_ratio, equilibrium_valid),
        entropy_metrics("driven_learned", driven, driven_ratio, driven_valid),
    ]
    exact_equilibrium_total = (
        validation.entropy_medium
        + (validation.energy_end - validation.energy_start) /
        args.equilibrium_temperature
    )
    exact_valid = np.ones(exact_equilibrium_total.size, dtype=bool)
    metrics.append(
        entropy_metrics(
            "equilibrium_exact_gibbs",
            Blocks(
                validation.stream, validation.block,
                validation.start, validation.end,
                validation.energy_start, validation.energy_end,
                exact_equilibrium_total, validation.balance_error,
            ),
            np.zeros(exact_equilibrium_total.size), exact_valid,
        )
    )
    write_csv(output / "total_entropy_metrics.csv", metrics)
    write_csv(
        output / "endpoint_stationarity.csv",
        stationarity_rows("equilibrium", equilibrium)
        + stationarity_rows("driven", driven),
    )

    driven_total = driven.entropy_medium[driven_valid] + driven_ratio[driven_valid]
    symmetry_rows, symmetry_summary = symmetry_diagnostic(driven_total)
    write_csv(output / "total_entropy_symmetry_bins.csv", symmetry_rows)
    write_csv(output / "total_entropy_symmetry_summary.csv", [symmetry_summary])

    eq_error = (
        equilibrium_ratio[equilibrium_valid]
        - validation_target[equilibrium_valid]
    )
    figure, axes = plt.subplots(1, 2, figsize=(10.0, 4.2))
    axes[0].hexbin(validation_target[equilibrium_valid],
                   equilibrium_ratio[equilibrium_valid], gridsize=50,
                   mincnt=1, bins="log")
    limits = np.quantile(validation_target[equilibrium_valid], [0.01, 0.99])
    axes[0].plot(limits, limits, "k--", lw=1)
    axes[0].set_xlim(limits)
    axes[0].set_ylim(limits)
    axes[0].set_xlabel(r"exact $\log\rho(X_0)/\rho(X_t)=\Delta E/T$")
    axes[0].set_ylabel("cross-fitted density log ratio")
    axes[0].set_title("equal-temperature validation")
    axes[1].hist(eq_error, bins=80, density=True, histtype="step")
    axes[1].set_xlabel("learned minus exact log ratio")
    axes[1].set_ylabel("density")
    axes[1].set_title(f"selected smoothing $\\sigma={sigma:g}$ bins")
    figure.tight_layout()
    figure.savefig(output / "equilibrium_density_validation.png", dpi=220)
    figure.savefig(output / "equilibrium_density_validation.pdf")
    plt.close(figure)

    usable = [row for row in symmetry_rows if int(row["usable"]) == 1]
    figure, axis = plt.subplots(figsize=(6.2, 4.4))
    if usable:
        x = np.asarray([float(row["entropy_center"]) for row in usable])
        y = np.asarray([float(row["log_probability_ratio"]) for row in usable])
        yerr = np.asarray([float(row["poisson_standard_error"]) for row in usable])
        axis.errorbar(x, y, yerr=yerr, fmt="o", ms=3, capsize=2,
                      label="cross-fitted n=2 estimate")
        axis.plot(x, x, "k--", label="detailed-FT reference")
    axis.set_xlabel(r"total entropy $s$")
    axis.set_ylabel(r"$\log[p(s)/p(-s)]$")
    axis.set_title("Exploratory total-entropy symmetry")
    axis.grid(alpha=0.25)
    if usable:
        axis.legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(output / "total_entropy_symmetry.png", dpi=220)
    figure.savefig(output / "total_entropy_symmetry.pdf")
    plt.close(figure)

    metadata = {
        "equilibrium_blocks": str(args.equilibrium_blocks),
        "driven_blocks": str(args.driven_blocks),
        "bins": bins,
        "sigma_grid": args.sigma_grid,
        "selected_sigma": sigma,
        "selection_rule": "minimum equal-temperature density-log-ratio RMSE",
        "interpretation": (
            "Driven total entropy is exploratory unless equilibrium density "
            "validation, support, and exponential-weight ESS are adequate."
        ),
    }
    (output / "analysis_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    print(f"Selected sigma={sigma:g}; wrote total-entropy pilot to {output}")


if __name__ == "__main__":
    main()
