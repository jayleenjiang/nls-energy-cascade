#!/usr/bin/env python3
"""Physics-constrained n=2 NESS total-entropy fluctuation-relation test.

The reduced stationary density q(u1,u2,theta) is represented by a normalized
action--Fourier exponential family.  The fixed D=2, K=1 feature family contains
the exact n=2 Gibbs density and is validated on an independent equal-
temperature data set before being cross-fitted to the driven blocks.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from analyze_total_entropy_n2 import (
    Blocks,
    entropy_metrics,
    load_blocks,
    stationarity_rows,
    subset,
    symmetry_diagnostic,
    wrap_angle,
    write_csv,
)


ACTION_DEGREE = 2
ANGULAR_HARMONICS = 1
QUADRATURE_GRID = (42, 42, 64)
RIDGE = 1.0e-8
MAX_ITERATIONS = 1000


def raw_features(points: np.ndarray) -> np.ndarray:
    u1 = points[:, 0]
    u2 = points[:, 1]
    theta = points[:, 2]
    action1 = np.exp(u1)
    action2 = np.exp(u2)
    columns = [u1, u2]
    monomials: list[np.ndarray] = []
    for total_degree in range(1, ACTION_DEGREE + 1):
        for power1 in range(total_degree + 1):
            power2 = total_degree - power1
            monomials.append(action1 ** power1 * action2 ** power2)
    columns.extend(monomials)
    harmonic_monomials = [np.ones(points.shape[0]), *monomials]
    for harmonic in range(1, ANGULAR_HARMONICS + 1):
        cosine = np.cos(harmonic * theta)
        sine = np.sin(harmonic * theta)
        for monomial in harmonic_monomials:
            columns.append(monomial * cosine)
            columns.append(monomial * sine)
    return np.column_stack(columns)


@dataclass
class DensityModel:
    coefficient: np.ndarray
    center: np.ndarray
    scale: np.ndarray
    low: np.ndarray
    high: np.ndarray
    optimizer_success: bool
    optimizer_iterations: int
    optimizer_gradient_norm: float
    training_support: float

    def evaluate(self, points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        transformed = points.copy()
        transformed[:, 2] = wrap_angle(transformed[:, 2])
        valid = (
            (transformed[:, 0] >= self.low[0])
            & (transformed[:, 0] <= self.high[0])
            & (transformed[:, 1] >= self.low[1])
            & (transformed[:, 1] <= self.high[1])
        )
        features = (raw_features(transformed) - self.center) / self.scale
        log_density = features @ self.coefficient
        log_density[~valid] = np.nan
        return log_density, valid


def feature_sum(points: np.ndarray, center: np.ndarray,
                scale: np.ndarray, chunk: int = 100_000) -> np.ndarray:
    total = np.zeros(center.size)
    for start in range(0, points.shape[0], chunk):
        features = raw_features(points[start:start + chunk])
        total += np.sum((features - center) / scale, axis=0)
    return total


def fit_density(points: np.ndarray) -> DensityModel:
    points = points.copy()
    points[:, 2] = wrap_angle(points[:, 2])
    low = np.quantile(points[:, :2], 0.0001, axis=0) - 0.5
    high = np.quantile(points[:, :2], 0.9999, axis=0) + 0.5
    inside = (
        (points[:, 0] >= low[0]) & (points[:, 0] <= high[0])
        & (points[:, 1] >= low[1]) & (points[:, 1] <= high[1])
    )
    training = points[inside]

    axes = [
        np.linspace(low[0], high[0], QUADRATURE_GRID[0]),
        np.linspace(low[1], high[1], QUADRATURE_GRID[1]),
        np.linspace(-np.pi, np.pi, QUADRATURE_GRID[2], endpoint=False),
    ]
    mesh = np.meshgrid(*axes, indexing="ij")
    grid = np.column_stack([value.ravel() for value in mesh])
    grid_features = raw_features(grid)
    center = np.mean(grid_features, axis=0)
    scale = np.std(grid_features, axis=0)
    scale[scale < 1.0e-12] = 1.0
    grid_features = (grid_features - center) / scale
    empirical_mean = feature_sum(training, center, scale) / training.shape[0]

    def objective(coefficient: np.ndarray) -> tuple[float, np.ndarray]:
        grid_log_density = grid_features @ coefficient
        normalization = logsumexp(grid_log_density)
        probability = np.exp(grid_log_density - normalization)
        value = (
            -float(empirical_mean @ coefficient) + float(normalization)
            + 0.5 * RIDGE * float(coefficient @ coefficient)
        )
        gradient = (
            -empirical_mean + grid_features.T @ probability
            + RIDGE * coefficient
        )
        return value, gradient

    result = minimize(
        objective,
        np.zeros(grid_features.shape[1]),
        jac=True,
        method="L-BFGS-B",
        options={
            "maxiter": MAX_ITERATIONS,
            "ftol": 1.0e-13,
            "gtol": 5.0e-7,
            "maxls": 50,
        },
    )
    return DensityModel(
        coefficient=np.asarray(result.x),
        center=center,
        scale=scale,
        low=low,
        high=high,
        optimizer_success=bool(result.success),
        optimizer_iterations=int(result.nit),
        optimizer_gradient_norm=float(np.linalg.norm(result.jac)),
        training_support=float(inside.mean()),
    )


def endpoint_log_ratio(model: DensityModel,
                       blocks: Blocks) -> tuple[np.ndarray, np.ndarray]:
    log_start, valid_start = model.evaluate(blocks.start)
    log_end, valid_end = model.evaluate(blocks.end)
    ratio = (
        log_start - log_end
        - blocks.start[:, 0] - blocks.start[:, 1]
        + blocks.end[:, 0] + blocks.end[:, 1]
    )
    valid = valid_start & valid_end & np.isfinite(ratio)
    ratio[~valid] = np.nan
    return ratio, valid


def crossfit_driven(blocks: Blocks) -> tuple[np.ndarray, np.ndarray,
                                                list[dict[str, object]]]:
    ratio = np.full(blocks.stream.size, np.nan)
    valid = np.zeros(blocks.stream.size, dtype=bool)
    diagnostics: list[dict[str, object]] = []
    for fold in (0, 1):
        test = blocks.stream % 2 == fold
        train = ~test
        model = fit_density(blocks.start[train])
        fold_blocks = subset(blocks, test)
        fold_ratio, fold_valid = endpoint_log_ratio(model, fold_blocks)
        indices = np.flatnonzero(test)
        ratio[indices[fold_valid]] = fold_ratio[fold_valid]
        valid[indices[fold_valid]] = True
        diagnostics.append({
            "dataset": "driven",
            "fold": fold,
            "training_samples": int(train.sum()),
            "evaluation_samples": int(test.sum()),
            "evaluation_support": float(fold_valid.mean()),
            "optimizer_success": int(model.optimizer_success),
            "optimizer_iterations": model.optimizer_iterations,
            "optimizer_gradient_norm": model.optimizer_gradient_norm,
            "training_support": model.training_support,
        })
    return ratio, valid, diagnostics


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("equilibrium_blocks", type=Path)
    parser.add_argument("driven_blocks", type=Path)
    parser.add_argument("--equilibrium-temperature", type=float, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)

    equilibrium = load_blocks(args.equilibrium_blocks)
    driven = load_blocks(args.driven_blocks)
    equilibrium_ids = np.unique(equilibrium.stream)
    if equilibrium_ids.size < 64 or np.unique(driven.stream).size < 64:
        raise ValueError("parametric endpoint control requires 64 streams")
    split = equilibrium_ids.size // 2
    training = subset(
        equilibrium, np.isin(equilibrium.stream, equilibrium_ids[:split]))
    validation = subset(
        equilibrium, np.isin(equilibrium.stream, equilibrium_ids[split:]))

    equilibrium_model = fit_density(training.start)
    equilibrium_ratio, equilibrium_valid = endpoint_log_ratio(
        equilibrium_model, validation)
    exact_ratio = (
        validation.energy_end - validation.energy_start
    ) / args.equilibrium_temperature
    error = equilibrium_ratio[equilibrium_valid] - exact_ratio[equilibrium_valid]
    slope = float(np.polyfit(
        exact_ratio[equilibrium_valid], equilibrium_ratio[equilibrium_valid], 1
    )[0])
    correlation = float(np.corrcoef(
        exact_ratio[equilibrium_valid], equilibrium_ratio[equilibrium_valid]
    )[0, 1])
    validation_row = {
        "samples_total": int(equilibrium_valid.size),
        "samples_supported": int(equilibrium_valid.sum()),
        "support_fraction": float(equilibrium_valid.mean()),
        "density_log_ratio_rmse": float(np.sqrt(np.mean(np.square(error)))),
        "density_log_ratio_mae": float(np.mean(np.abs(error))),
        "density_log_ratio_slope": slope,
        "density_log_ratio_correlation": correlation,
    }
    write_csv(output / "equilibrium_density_validation.csv", [validation_row])

    driven_ratio, driven_valid, model_diagnostics = crossfit_driven(driven)
    model_diagnostics.insert(0, {
        "dataset": "equilibrium",
        "fold": "heldout",
        "training_samples": int(training.stream.size),
        "evaluation_samples": int(validation.stream.size),
        "evaluation_support": float(equilibrium_valid.mean()),
        "optimizer_success": int(equilibrium_model.optimizer_success),
        "optimizer_iterations": equilibrium_model.optimizer_iterations,
        "optimizer_gradient_norm": equilibrium_model.optimizer_gradient_norm,
        "training_support": equilibrium_model.training_support,
    })
    write_csv(output / "density_model_diagnostics.csv", model_diagnostics)

    metrics = [
        entropy_metrics("equilibrium_parametric", validation,
                        equilibrium_ratio, equilibrium_valid),
        entropy_metrics("driven_parametric", driven,
                        driven_ratio, driven_valid),
    ]
    exact_total = (
        validation.entropy_medium
        + (validation.energy_end - validation.energy_start)
        / args.equilibrium_temperature
    )
    metrics.append(entropy_metrics(
        "equilibrium_exact_gibbs",
        Blocks(
            validation.stream, validation.block, validation.start,
            validation.end, validation.energy_start, validation.energy_end,
            exact_total, validation.balance_error,
        ),
        np.zeros(exact_total.size), np.ones(exact_total.size, dtype=bool),
    ))
    write_csv(output / "total_entropy_metrics.csv", metrics)
    write_csv(
        output / "endpoint_stationarity.csv",
        stationarity_rows("equilibrium", equilibrium)
        + stationarity_rows("driven", driven),
    )

    driven_total = (
        driven.entropy_medium[driven_valid] + driven_ratio[driven_valid]
    )
    symmetry_rows, symmetry_summary = symmetry_diagnostic(
        driven_total, bins=80, min_count=20)
    write_csv(output / "total_entropy_symmetry_bins.csv", symmetry_rows)
    write_csv(output / "total_entropy_symmetry_summary.csv", [symmetry_summary])

    figure, axes = plt.subplots(1, 2, figsize=(10.6, 4.2))
    limits = np.quantile(exact_ratio[equilibrium_valid], [0.0025, 0.9975])
    axes[0].hexbin(
        exact_ratio[equilibrium_valid], equilibrium_ratio[equilibrium_valid],
        gridsize=60, mincnt=1, bins="log")
    axes[0].plot(limits, limits, "k--")
    axes[0].set_xlim(limits)
    axes[0].set_ylim(limits)
    axes[0].set_xlabel("exact equilibrium endpoint ratio")
    axes[0].set_ylabel("parametric endpoint ratio")
    axes[0].set_title(
        f"RMSE={validation_row['density_log_ratio_rmse']:.4f}, "
        f"slope={slope:.4f}")

    usable = [row for row in symmetry_rows if int(row["usable"]) == 1]
    x = np.asarray([float(row["entropy_center"]) for row in usable])
    y = np.asarray([float(row["log_probability_ratio"]) for row in usable])
    yerr = np.asarray([float(row["poisson_standard_error"]) for row in usable])
    axes[1].errorbar(x, y, yerr=yerr, fmt="o", ms=3, capsize=2,
                     label="driven n=2 total entropy")
    if x.size:
        axes[1].plot(x, x, "k--", label="detailed-FT reference")
    axes[1].set_xlabel(r"total entropy $s$")
    axes[1].set_ylabel(r"$\log[p(s)/p(-s)]$")
    axes[1].legend(frameon=False)
    figure.tight_layout()
    figure.savefig(output / "parametric_total_entropy_ft.png", dpi=220)
    figure.savefig(output / "parametric_total_entropy_ft.pdf")
    plt.close(figure)

    metadata = {
        "equilibrium_blocks": str(args.equilibrium_blocks),
        "driven_blocks": str(args.driven_blocks),
        "feature_family": {
            "action_degree": ACTION_DEGREE,
            "angular_harmonics": ANGULAR_HARMONICS,
            "quadrature_grid": QUADRATURE_GRID,
            "ridge": RIDGE,
        },
        "selection_boundary": (
            "All feature and optimizer settings were fixed from equilibrium "
            "validation before applying the model to driven endpoints."
        ),
    }
    (output / "analysis_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n")
    print(f"Wrote parametric n=2 total-entropy analysis to {output}")


if __name__ == "__main__":
    main()
