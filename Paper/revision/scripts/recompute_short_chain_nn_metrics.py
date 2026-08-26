#!/usr/bin/env python3
"""Recompute short-chain neural-network diagnostics from saved models/data.

This is the non-notebook rerun companion for the three-mode Fokker--Planck
section.  It loads the archived Keras models under ``KDE/`` and recomputes the
diagnostics that previously lived only in notebook cells:

* equilibrium-vs-Gibbs slice errors;
* nonequilibrium angular symmetry breaking;
* angular-width diagnostic;
* middle-mode current-balance check from MC/KDE boxes;
* eigenfunction surrogate data-fit error.

The script does not retrain any model.  It verifies the saved model artifacts
against the local source data and writes a JSON file suitable for manuscript
auditing.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import tensorflow as tf


ROOT = Path(__file__).resolve().parents[3]
REVISION = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = REVISION / "short_chain_nn_rerun_metrics.json"


@tf.keras.utils.register_keras_serializable()
def periodic_encode(x: tf.Tensor) -> tf.Tensor:
    """Periodic angle encoding used by the saved Keras models."""

    i_vars = x[:, :3]
    cos_th1 = tf.math.cos(x[:, 3:4])
    sin_th1 = tf.math.sin(x[:, 3:4])
    cos_th3 = tf.math.cos(x[:, 4:5])
    sin_th3 = tf.math.sin(x[:, 4:5])
    return tf.concat([i_vars, cos_th1, sin_th1, cos_th3, sin_th3], axis=-1)


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def file_record(path: Path) -> dict[str, Any]:
    record: dict[str, Any] = {
        "path": rel(path),
        "exists": path.exists(),
    }
    if path.exists():
        record["bytes"] = path.stat().st_size
        record["sha256"] = sha256_file(path)
    return record


def load_model(path: Path) -> tf.keras.Model:
    return tf.keras.models.load_model(path, custom_objects={"periodic_encode": periodic_encode})


def log_normalization(density_path: Path) -> tuple[float, float]:
    dens = np.loadtxt(density_path, dtype=np.float32)
    log_dens = np.log(dens[dens > 0])
    return float(log_dens.mean()), float(log_dens.std())


def predict_density(model: tf.keras.Model, x: np.ndarray, mu: float, sigma: float, batch_size: int) -> np.ndarray:
    pred = model.predict(x.astype(np.float32), batch_size=batch_size, verbose=0).reshape(-1)
    return np.exp(sigma * pred + mu)


def equilibrium_validation(model_eq: tf.keras.Model, mu_eq: float, sigma_eq: float, n_grid: int, batch_size: int) -> list[dict[str, Any]]:
    theta = np.linspace(-np.pi, np.pi, n_grid, dtype=np.float32)
    th1, th3 = np.meshgrid(theta, theta)
    rows: list[dict[str, Any]] = []
    for i_val in [0.5, 1.0, 2.0]:
        sl = np.column_stack(
            [
                np.full(n_grid * n_grid, i_val, dtype=np.float32),
                np.full(n_grid * n_grid, i_val, dtype=np.float32),
                np.full(n_grid * n_grid, i_val, dtype=np.float32),
                th1.ravel(),
                th3.ravel(),
            ]
        )
        p_nn = predict_density(model_eq, sl, mu_eq, sigma_eq, batch_size).reshape(n_grid, n_grid)
        total_action = 3.0 * i_val
        # Formula used in the archived equilibrium-validation notebook cell.
        hamiltonian = (
            total_action**2
            - 1.5 * i_val**2
            + 2.0 * i_val**2 * np.cos(th1)
            + 2.0 * i_val**2 * np.cos(th3)
        )
        p_gibbs = np.exp(-hamiltonian / 10.0)
        p_gibbs_norm = p_gibbs / p_gibbs.sum() * p_nn.sum()
        rel_err = np.abs((p_nn - p_gibbs_norm) / p_gibbs_norm)
        log_ratio = np.log(p_nn) - np.log(p_gibbs_norm)
        rows.append(
            {
                "I": i_val,
                "mean_relative_error_percent": float(100.0 * rel_err.mean()),
                "median_relative_error_percent": float(100.0 * np.median(rel_err)),
                "max_relative_error_percent": float(100.0 * rel_err.max()),
                "log_ratio_min": float(log_ratio.min()),
                "log_ratio_max": float(log_ratio.max()),
                "grid": n_grid,
            }
        )
    return rows


def symmetry_breaking(model_neq: tf.keras.Model, mu_neq: float, sigma_neq: float, n_grid: int, batch_size: int) -> dict[str, Any]:
    theta = np.linspace(-np.pi, np.pi, n_grid, dtype=np.float32)
    th1, th3 = np.meshgrid(theta, theta)
    i_val = 2.0
    sl = np.column_stack(
        [
            np.full(n_grid * n_grid, i_val, dtype=np.float32),
            np.full(n_grid * n_grid, i_val, dtype=np.float32),
            np.full(n_grid * n_grid, i_val, dtype=np.float32),
            th1.ravel(),
            th3.ravel(),
        ]
    )
    p = predict_density(model_neq, sl, mu_neq, sigma_neq, batch_size).reshape(n_grid, n_grid)
    asym = (p - p.T) / (p + p.T + 1e-15)
    high_density_mask = (p + p.T) >= 0.20 * np.max(p + p.T)
    return {
        "I": i_val,
        "grid": n_grid,
        "unmasked_max_abs_asymmetry_percent": float(100.0 * np.abs(asym).max()),
        "unmasked_mean_abs_asymmetry_percent": float(100.0 * np.abs(asym).mean()),
        "high_density_mask_threshold": "p+p.T >= 20% of max(p+p.T)",
        "high_density_bins": int(high_density_mask.sum()),
        "masked_max_abs_asymmetry_percent": float(100.0 * np.abs(asym[high_density_mask]).max()),
        "masked_mean_abs_asymmetry_percent": float(100.0 * np.abs(asym[high_density_mask]).mean()),
    }


def circular_std(theta: np.ndarray, weights: np.ndarray) -> float:
    weights = np.asarray(weights, dtype=np.float64)
    weights = weights / weights.sum()
    resultant = np.abs(np.sum(weights * np.exp(1j * theta)))
    resultant = float(np.clip(resultant, 1e-15, 1.0))
    return float(np.sqrt(-2.0 * np.log(resultant)))


def angular_width_diagnostic(
    model_neq: tf.keras.Model,
    mu_neq: float,
    sigma_neq: float,
    n_theta: int,
    marginal_grid: int,
    batch_size: int,
) -> list[dict[str, Any]]:
    theta = np.linspace(-np.pi, np.pi, n_theta, dtype=np.float32)
    other = np.linspace(-np.pi, np.pi, marginal_grid, dtype=np.float32)
    rows: list[dict[str, Any]] = []
    for i_val in [0.5, 1.0, 2.0, 3.0, 4.0]:
        # Marginal in theta_1 by summing over theta_3.
        th1 = np.tile(theta, marginal_grid)
        th3 = np.repeat(other, n_theta)
        sl1 = np.column_stack(
            [
                np.full(th1.size, i_val, dtype=np.float32),
                np.full(th1.size, i_val, dtype=np.float32),
                np.full(th1.size, i_val, dtype=np.float32),
                th1,
                th3,
            ]
        )
        p1_grid = predict_density(model_neq, sl1, mu_neq, sigma_neq, batch_size).reshape(marginal_grid, n_theta)
        p_th1 = p1_grid.sum(axis=0)

        # Marginal in theta_3 by summing over theta_1.
        th1b = np.repeat(other, n_theta)
        th3b = np.tile(theta, marginal_grid)
        sl3 = np.column_stack(
            [
                np.full(th3b.size, i_val, dtype=np.float32),
                np.full(th3b.size, i_val, dtype=np.float32),
                np.full(th3b.size, i_val, dtype=np.float32),
                th1b,
                th3b,
            ]
        )
        p3_grid = predict_density(model_neq, sl3, mu_neq, sigma_neq, batch_size).reshape(marginal_grid, n_theta)
        p_th3 = p3_grid.sum(axis=0)

        sigma_1 = circular_std(theta.astype(np.float64), p_th1)
        sigma_3 = circular_std(theta.astype(np.float64), p_th3)
        rows.append(
            {
                "I": i_val,
                "sigma_theta1": sigma_1,
                "sigma_theta3": sigma_3,
                "sigma3_over_sigma1": sigma_3 / sigma_1,
                "sqrt_T3_over_T1": float(np.sqrt(8.0 / 2.0)),
                "theta_grid": n_theta,
                "marginal_grid": marginal_grid,
            }
        )
    return rows


def phase_locking_diagnostic(model_neq: tf.keras.Model, mu_neq: float, sigma_neq: float, n_theta: int, batch_size: int) -> dict[str, Any]:
    theta = np.linspace(-np.pi, np.pi, n_theta, dtype=np.float32)
    gamma = 0.1
    stable_branch = np.pi - np.arcsin((np.sqrt(4.0 * gamma**2 + 3.0) - gamma) / (2.0 * (gamma**2 + 1.0)))
    rows: list[dict[str, Any]] = []
    for i_val in [0.5, 1.0, 2.0, 3.0, 4.0, 5.0]:
        sl = np.column_stack(
            [
                np.full(n_theta, i_val, dtype=np.float32),
                np.full(n_theta, i_val, dtype=np.float32),
                np.full(n_theta, i_val, dtype=np.float32),
                theta,
                np.zeros(n_theta, dtype=np.float32),
            ]
        )
        p_nn = predict_density(model_neq, sl, mu_neq, sigma_neq, batch_size)
        rows.append(
            {
                "I": i_val,
                "argmax_density_theta_rad": float(theta[int(np.argmax(p_nn))]),
                "argmax_density_theta_deg": float(np.degrees(theta[int(np.argmax(p_nn))])),
                "legacy_notebook_argmin_density_theta_deg": float(np.degrees(theta[int(np.argmin(p_nn))])),
                "note": "argmax is the density peak; legacy notebook printed argmin and is retained only for provenance comparison.",
            }
        )
    return {
        "gamma": gamma,
        "stable_branch_theta0_rad": float(stable_branch),
        "stable_branch_theta0_deg": float(np.degrees(stable_branch)),
        "rows": rows,
    }


def mc_current_balance(boxes_path: Path, density_path: Path) -> dict[str, Any]:
    x_mc = np.loadtxt(boxes_path, dtype=np.float32)
    y_mc = np.loadtxt(density_path, dtype=np.float32)
    mask = y_mc > 0
    i1, i2, i3, th1, th3 = x_mc[mask, 0], x_mc[mask, 1], x_mc[mask, 2], x_mc[mask, 3], x_mc[mask, 4]
    weights = y_mc[mask]
    f12 = float(np.sum(i2 * i1 * np.sin(th1) * weights) / np.sum(weights))
    f23 = float(np.sum(i2 * i3 * np.sin(th3) * weights) / np.sum(weights))
    return {
        "nonzero_density_boxes": int(mask.sum()),
        "E_I2_I1_sin_theta1": f12,
        "E_I2_I3_sin_theta3": f23,
        "sum": f12 + f23,
        "relative_imbalance_vs_first_current": abs(f12 + f23) / (abs(f12) + 1e-10),
    }


def eigen_data_fit(model_eigen: tf.keras.Model, x_path: Path, q_path: Path, batch_size: int) -> dict[str, Any]:
    x_data = np.loadtxt(x_path, dtype=np.float32)
    q_data = np.loadtxt(q_path, dtype=np.float32).reshape(-1)
    q_mean = float(q_data.mean())
    q_std = float(q_data.std())
    pred_scaled = model_eigen.predict(x_data, batch_size=batch_size, verbose=0).reshape(-1)
    q_pred = q_std * pred_scaled + q_mean
    rmse = float(np.sqrt(np.mean((q_data - q_pred) ** 2)))
    return {
        "data_points": int(len(q_data)),
        "Q1_min": float(q_data.min()),
        "Q1_max": float(q_data.max()),
        "Q1_mean": q_mean,
        "Q1_std": q_std,
        "rmse": rmse,
        "relative_rmse": float(rmse / q_std),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--grid", type=int, default=100)
    parser.add_argument("--phase-grid", type=int, default=500)
    parser.add_argument("--marginal-grid", type=int, default=20)
    parser.add_argument("--batch-size", type=int, default=8192)
    args = parser.parse_args()

    neq_model_path = ROOT / "KDE/4:15_NN/h5_files/final.keras"
    eq_model_path = ROOT / "KDE/4:15_NN/h5_files_eq/final.keras"
    eigen_model_path = ROOT / "KDE/h5_files_eigen/final.keras"
    neq_density_path = ROOT / "KDE/4:15_NN/NLS_FP_density.txt"
    eq_density_path = ROOT / "KDE/4:15_NN/eq/NLS_FP_density.txt"
    boxes_path = ROOT / "KDE/4:15_NN/NLS_FP_boxes.txt"
    eigen_x_path = ROOT / "KDE/backward_NLS_X.txt"
    eigen_q_path = ROOT / "KDE/backward_NLS_Q1.txt"

    model_neq = load_model(neq_model_path)
    model_eq = load_model(eq_model_path)
    model_eigen = load_model(eigen_model_path)
    mu_neq, sigma_neq = log_normalization(neq_density_path)
    mu_eq, sigma_eq = log_normalization(eq_density_path)

    payload = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "script": rel(Path(__file__).resolve()),
        "runtime": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "tensorflow": tf.__version__,
            "keras": tf.keras.__version__,
        },
        "source_files": {
            "nonequilibrium_model": file_record(neq_model_path),
            "equilibrium_model": file_record(eq_model_path),
            "eigen_model": file_record(eigen_model_path),
            "nonequilibrium_boxes": file_record(boxes_path),
            "nonequilibrium_density": file_record(neq_density_path),
            "equilibrium_density": file_record(eq_density_path),
            "eigen_X": file_record(eigen_x_path),
            "eigen_Q1": file_record(eigen_q_path),
        },
        "normalization": {
            "nonequilibrium_log_density_mean": mu_neq,
            "nonequilibrium_log_density_std": sigma_neq,
            "equilibrium_log_density_mean": mu_eq,
            "equilibrium_log_density_std": sigma_eq,
        },
        "equilibrium_validation": equilibrium_validation(model_eq, mu_eq, sigma_eq, args.grid, args.batch_size),
        "symmetry_breaking": symmetry_breaking(model_neq, mu_neq, sigma_neq, args.grid, args.batch_size),
        "angular_width_diagnostic": angular_width_diagnostic(
            model_neq,
            mu_neq,
            sigma_neq,
            args.grid,
            args.marginal_grid,
            args.batch_size,
        ),
        "phase_locking_diagnostic": phase_locking_diagnostic(model_neq, mu_neq, sigma_neq, args.phase_grid, args.batch_size),
        "mc_current_balance": mc_current_balance(boxes_path, neq_density_path),
        "eigen_surrogate_data_fit": eigen_data_fit(model_eigen, eigen_x_path, eigen_q_path, args.batch_size),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({"wrote": rel(args.output), "bytes": args.output.stat().st_size}, indent=2))


if __name__ == "__main__":
    main()
