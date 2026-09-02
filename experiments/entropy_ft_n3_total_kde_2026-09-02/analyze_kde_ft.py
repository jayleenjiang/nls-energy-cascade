#!/usr/bin/env python3
"""Frozen held-out periodic-KDE audit and n=3 total-entropy FT analysis."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter


EXPECTED_HEADER = [
    "stream_id", "block_id", "q_left", "q_right", "delta_energy",
    "entropy_medium", "entropy_rate", "action_current",
    "energy_balance_error", "start_I1", "start_I2", "start_I3",
    "start_theta1", "start_theta3", "end_I1", "end_I2", "end_I3",
    "end_theta1", "end_theta3",
]
EXPECTED_STREAMS = 128
EXPECTED_BLOCKS_PER_STREAM = 7813
EXPECTED_ROWS = EXPECTED_STREAMS * EXPECTED_BLOCKS_PER_STREAM
BLOCK_TIME = 20.0
GRID_SHAPE = (48, 48, 48, 32, 32)
GRID_TRUNCATE = 4.0
BOOTSTRAPS = 2000
BOOTSTRAP_SEED = 2_026_090_291
MIN_PAIR_COUNT = 20
MIN_NEGATIVE = 1000
MIN_PAIRS = 8
KDE_ENDPOINT_RMSE_MAX = 0.15
KDE_INCREMENT_RMSE_MAX = 0.10
KDE_TAIL_INCREMENT_RMSE_MAX = 0.25
KDE_INCREMENT_Q99_MAX = 0.50
IFT_ESS_MIN = 1000.0
IFT_MAX_SHARE_MAX = 0.01
IFT_JACKKNIFE_CHANGE_MAX = 0.10
EXPECTED_DECOMPRESSED_HASHES = {
    "driven": "4f728b3d0e007d704d90734b0888c00ec05b60f09385c6cfd079f3417d7a088f",
    "T6": "dad1506fad5edda62620905a233f00f01f5522c90b55dd22724d861850b2f4eb",
    "T10": "83b5658c6b0a7639ed2ce57d4b6ff1761fe603babdda2dbd3b7450760128d1cc",
}


@dataclass
class Dataset:
    stream: np.ndarray
    block: np.ndarray
    entropy_medium: np.ndarray
    start: np.ndarray
    end: np.ndarray
    balance: np.ndarray


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]),
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256_decompressed(path: Path) -> str:
    process = subprocess.Popen(["zstd", "-dc", str(path)],
                               stdout=subprocess.PIPE)
    assert process.stdout is not None
    digest = hashlib.sha256()
    for chunk in iter(lambda: process.stdout.read(1024 * 1024), b""):
        digest.update(chunk)
    process.stdout.close()
    if process.wait() != 0:
        raise RuntimeError(f"zstd hash stream failed: {path}")
    return digest.hexdigest()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_archive(path: Path) -> Dataset:
    process = subprocess.Popen(["zstd", "-dc", str(path)],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               text=True,
                               bufsize=1024 * 1024)
    assert process.stdout is not None
    header = next(csv.reader([process.stdout.readline().rstrip("\n")]))
    if header != EXPECTED_HEADER:
        raise RuntimeError(f"unexpected header in {path}: {header}")
    matrix = np.loadtxt(process.stdout, delimiter=",", dtype=np.float64)
    process.stdout.close()
    stderr = process.stderr.read() if process.stderr is not None else ""
    if process.wait() != 0:
        raise RuntimeError(f"zstd failed for {path}: {stderr}")
    if matrix.shape != (EXPECTED_ROWS, len(EXPECTED_HEADER)):
        raise RuntimeError(f"unexpected matrix shape {matrix.shape} in {path}")
    if not np.isfinite(matrix).all():
        raise RuntimeError(f"nonfinite raw value in {path}")

    stream = matrix[:, 0].astype(np.int16)
    block = matrix[:, 1].astype(np.int32)
    expected_stream = np.repeat(np.arange(EXPECTED_STREAMS),
                                EXPECTED_BLOCKS_PER_STREAM)
    expected_block = np.tile(np.arange(EXPECTED_BLOCKS_PER_STREAM),
                             EXPECTED_STREAMS)
    if not np.array_equal(stream, expected_stream):
        raise RuntimeError(f"stream order failure in {path}")
    if not np.array_equal(block, expected_block):
        raise RuntimeError(f"block order failure in {path}")
    entropy = matrix[:, 5].copy()
    balance = matrix[:, 8].copy()
    start = matrix[:, 9:14].copy()
    end = matrix[:, 14:19].copy()

    identity_error = np.max(np.abs(
        entropy + matrix[:, 2] / 10.0 + matrix[:, 3] / 2.0
    )) if "driven" in path.name.lower() else float("nan")
    del matrix

    if np.min(start[:, :3]) <= 0.0 or np.min(end[:, :3]) <= 0.0:
        raise RuntimeError(f"nonpositive action in {path}")
    if (np.min(start[:, 3:]) < -math.pi or
            np.max(start[:, 3:]) >= math.pi or
            np.min(end[:, 3:]) < -math.pi or
            np.max(end[:, 3:]) >= math.pi):
        raise RuntimeError(f"angle outside [-pi,pi) in {path}")
    linear_error = float(np.max(np.abs(
        end.reshape(EXPECTED_STREAMS, EXPECTED_BLOCKS_PER_STREAM, 5)[:, :-1, :3]
        - start.reshape(EXPECTED_STREAMS, EXPECTED_BLOCKS_PER_STREAM, 5)[:, 1:, :3]
    )))
    angle_delta = (
        end.reshape(EXPECTED_STREAMS, EXPECTED_BLOCKS_PER_STREAM, 5)[:, :-1, 3:]
        - start.reshape(EXPECTED_STREAMS, EXPECTED_BLOCKS_PER_STREAM, 5)[:, 1:, 3:]
        + math.pi
    ) % (2.0 * math.pi) - math.pi
    angular_error = float(np.max(np.abs(angle_delta)))
    if linear_error > 1e-10 or angular_error > 1e-10:
        raise RuntimeError(f"endpoint continuity failure in {path}")
    print(json.dumps({
        "event": "read_complete", "path": str(path),
        "rows": EXPECTED_ROWS, "entropy_identity_error_if_driven": identity_error,
        "linear_endpoint_error": linear_error,
        "angular_endpoint_error": angular_error,
    }), flush=True)
    return Dataset(stream, block, entropy, start, end, balance)


def transform(state: np.ndarray) -> np.ndarray:
    result = np.empty_like(state, dtype=np.float64)
    result[:, :3] = np.log(state[:, :3])
    result[:, 3:] = (state[:, 3:] + math.pi) % (2.0 * math.pi) - math.pi
    return result


def energy(state: np.ndarray) -> np.ndarray:
    actions = state[:, :3]
    total = actions.sum(axis=1)
    return (0.5 * total * total
            - 0.25 * np.sum(actions * actions, axis=1)
            + actions[:, 0] * actions[:, 1] * np.cos(state[:, 3])
            + actions[:, 1] * actions[:, 2] * np.cos(state[:, 4]))


def circular_scale(values: np.ndarray) -> float:
    resultant = math.hypot(float(np.mean(np.cos(values))),
                           float(np.mean(np.sin(values))))
    resultant = min(1.0, max(resultant, 1e-15))
    return min(math.sqrt(max(0.0, -2.0 * math.log(resultant))),
               math.pi / math.sqrt(3.0))


def statistical_inefficiencies(points: np.ndarray,
                               stream_ids: np.ndarray) -> tuple[float, dict]:
    unique = np.unique(stream_ids)
    if unique.size < 2:
        raise RuntimeError("at least two streams required for autocorrelation")
    counts = np.array([np.count_nonzero(stream_ids == value) for value in unique])
    if not np.all(counts == counts[0]):
        raise RuntimeError("unequal stream lengths in bandwidth estimator")
    length = int(counts[0])
    transformed = [
        ("logI1", points[:, 0]), ("logI2", points[:, 1]),
        ("logI3", points[:, 2]), ("cos_theta1", np.cos(points[:, 3])),
        ("sin_theta1", np.sin(points[:, 3])),
        ("cos_theta3", np.cos(points[:, 4])),
        ("sin_theta3", np.sin(points[:, 4])),
    ]
    nfft = 1 << (2 * length - 1).bit_length()
    factors = {}
    for name, flat in transformed:
        series = np.stack([flat[stream_ids == value] for value in unique])
        series = series - series.mean(axis=1, keepdims=True)
        spectrum = np.fft.rfft(series, n=nfft, axis=1)
        autocov = np.fft.irfft(spectrum * np.conjugate(spectrum),
                               n=nfft, axis=1)[:, :length].real
        autocov /= np.arange(length, 0, -1, dtype=np.float64)[None, :]
        pooled = autocov.mean(axis=0)
        if not (pooled[0] > 0.0):
            raise RuntimeError(f"zero variance in autocorrelation series {name}")
        rho = pooled / pooled[0]
        cutoff = 1
        pair = 0
        while 2 * pair + 2 < length:
            left = 2 * pair + 1
            right = 2 * pair + 2
            if rho[left] + rho[right] <= 0.0:
                break
            cutoff = right + 1
            pair += 1
        g = max(1.0, 1.0 + 2.0 * float(np.sum(rho[1:cutoff])))
        factors[name] = g
    return max(factors.values()), factors


def bandwidth(points: np.ndarray,
              stream_ids: np.ndarray) -> tuple[np.ndarray, dict]:
    maximum_g, factors = statistical_inefficiencies(points, stream_ids)
    effective_n = points.shape[0] / maximum_g
    factor = effective_n ** (-1.0 / 9.0)
    scales = np.empty(5, dtype=np.float64)
    scales[:3] = np.std(points[:, :3], axis=0, ddof=1)
    scales[3] = circular_scale(points[:, 3])
    scales[4] = circular_scale(points[:, 4])
    values = factor * scales
    if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
        raise RuntimeError(f"invalid Scott bandwidth {values}")
    audit = {
        "raw_training_points": int(points.shape[0]),
        "maximum_statistical_inefficiency": float(maximum_g),
        "effective_sample_size": float(effective_n),
        "scott_factor": float(factor),
        **{f"g_{name}": float(value) for name, value in factors.items()},
    }
    return values, audit


def global_grid_bounds(end_z: np.ndarray,
                       stream_ids: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    global_bw, _ = bandwidth(end_z, stream_ids)
    lower = np.empty(5, dtype=np.float64)
    upper = np.empty(5, dtype=np.float64)
    lower[:3] = np.min(end_z[:, :3], axis=0) - 4.0 * global_bw[:3]
    upper[:3] = np.max(end_z[:, :3], axis=0) + 4.0 * global_bw[:3]
    lower[3:] = -math.pi
    upper[3:] = math.pi
    return lower, upper


def grid_indices(points: np.ndarray, lower: np.ndarray,
                 spacing: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    coordinates = (points - lower) / spacing
    indices = np.floor(coordinates).astype(np.int64)
    indices[:, 3] %= GRID_SHAPE[3]
    indices[:, 4] %= GRID_SHAPE[4]
    valid = np.ones(points.shape[0], dtype=bool)
    for axis in range(3):
        valid &= (indices[:, axis] >= 0) & (indices[:, axis] < GRID_SHAPE[axis])
    return indices, valid


def fit_grid(training: np.ndarray, training_streams: np.ndarray,
             lower: np.ndarray, upper: np.ndarray,
             label: str) -> tuple[np.ndarray, dict]:
    bw, correlation_audit = bandwidth(training, training_streams)
    spacing = (upper - lower) / np.asarray(GRID_SHAPE, dtype=np.float64)
    ratio = bw / spacing
    corrected_variance = ratio * ratio - 1.0 / 12.0
    if np.any(corrected_variance <= 0.0):
        raise RuntimeError(
            f"frozen grid cannot resolve bandwidth for {label}: "
            f"bandwidth/spacing={ratio.tolist()}")
    sigma = np.sqrt(corrected_variance)
    indices, valid = grid_indices(training, lower, spacing)
    if not np.all(valid):
        raise RuntimeError(f"training point outside global grid for {label}")
    flat_index = np.ravel_multi_index(indices.T, GRID_SHAPE)
    grid = np.zeros(int(np.prod(GRID_SHAPE)), dtype=np.float32)
    np.add.at(grid, flat_index, 1.0)
    grid = grid.reshape(GRID_SHAPE)
    gaussian_filter(
        grid, sigma=sigma, output=grid,
        mode=("constant", "constant", "constant", "wrap", "wrap"),
        cval=0.0, truncate=GRID_TRUNCATE)
    mass_before_normalization = float(grid.sum(dtype=np.float64))
    if not (mass_before_normalization > 0.0):
        raise RuntimeError(f"zero KDE mass for {label}")
    cell_volume = float(np.prod(spacing))
    grid /= np.float32(mass_before_normalization * cell_volume)
    metadata = {
        "label": label,
        "training_points": int(training.shape[0]),
        "bandwidth_logI1": float(bw[0]),
        "bandwidth_logI2": float(bw[1]),
        "bandwidth_logI3": float(bw[2]),
        "bandwidth_theta1": float(bw[3]),
        "bandwidth_theta3": float(bw[4]),
        "spacing_logI1": float(spacing[0]),
        "spacing_logI2": float(spacing[1]),
        "spacing_logI3": float(spacing[2]),
        "spacing_theta1": float(spacing[3]),
        "spacing_theta3": float(spacing[4]),
        "bw_over_spacing_min": float(np.min(ratio)),
        "raw_smoothed_mass": mass_before_normalization,
        "cell_volume": cell_volume,
        **correlation_audit,
    }
    return grid, metadata


def interpolate_density(grid: np.ndarray, points: np.ndarray,
                        lower: np.ndarray, upper: np.ndarray,
                        chunk_size: int = 100_000) -> tuple[np.ndarray, np.ndarray]:
    spacing = (upper - lower) / np.asarray(GRID_SHAPE, dtype=np.float64)
    result = np.full(points.shape[0], np.nan, dtype=np.float64)
    supported = np.zeros(points.shape[0], dtype=bool)
    flat_grid = grid.reshape(-1)
    strides = np.array([
        GRID_SHAPE[1] * GRID_SHAPE[2] * GRID_SHAPE[3] * GRID_SHAPE[4],
        GRID_SHAPE[2] * GRID_SHAPE[3] * GRID_SHAPE[4],
        GRID_SHAPE[3] * GRID_SHAPE[4], GRID_SHAPE[4], 1,
    ], dtype=np.int64)
    for first in range(0, points.shape[0], chunk_size):
        last = min(points.shape[0], first + chunk_size)
        current = points[first:last]
        coordinate = (current - lower) / spacing - 0.5
        base = np.floor(coordinate).astype(np.int64)
        fraction = coordinate - base
        valid = np.ones(current.shape[0], dtype=bool)
        for axis in range(3):
            valid &= (base[:, axis] >= 0) & (base[:, axis] + 1 < GRID_SHAPE[axis])
        base[:, 3] %= GRID_SHAPE[3]
        base[:, 4] %= GRID_SHAPE[4]
        density = np.zeros(current.shape[0], dtype=np.float64)
        for corner in range(32):
            index = np.zeros(current.shape[0], dtype=np.int64)
            weight = np.ones(current.shape[0], dtype=np.float64)
            for axis in range(5):
                upper_corner = (corner >> axis) & 1
                if upper_corner:
                    axis_index = base[:, axis] + 1
                    weight *= fraction[:, axis]
                else:
                    axis_index = base[:, axis]
                    weight *= 1.0 - fraction[:, axis]
                if axis >= 3:
                    axis_index %= GRID_SHAPE[axis]
                index += axis_index * strides[axis]
            density += weight * flat_grid[index]
        valid &= np.isfinite(density) & (density > 0.0)
        result[first:last][valid] = density[valid]
        supported[first:last] = valid
    return result, supported


def crossfit_density(dataset: Dataset, label: str) -> tuple[np.ndarray, np.ndarray, list[dict], dict]:
    start_z = transform(dataset.start)
    end_z = transform(dataset.end)
    lower, upper = global_grid_bounds(end_z, dataset.stream)
    log_start = np.full(EXPECTED_ROWS, np.nan, dtype=np.float64)
    log_end = np.full(EXPECTED_ROWS, np.nan, dtype=np.float64)
    metadata: list[dict] = []
    for evaluation_parity in (0, 1):
        evaluate = dataset.stream % 2 == evaluation_parity
        train = ~evaluate
        grid, row = fit_grid(
            end_z[train], dataset.stream[train], lower, upper,
            f"{label}_train_parity_{1-evaluation_parity}")
        density_start, support_start = interpolate_density(
            grid, start_z[evaluate], lower, upper)
        density_end, support_end = interpolate_density(
            grid, end_z[evaluate], lower, upper)
        positions = np.flatnonzero(evaluate)
        valid = support_start & support_end
        log_start[positions[valid]] = (
            np.log(density_start[valid]) - np.sum(start_z[evaluate][valid, :3], axis=1))
        log_end[positions[valid]] = (
            np.log(density_end[valid]) - np.sum(end_z[evaluate][valid, :3], axis=1))
        row.update({
            "evaluation_parity": evaluation_parity,
            "evaluation_blocks": int(evaluate.sum()),
            "unsupported_start": int((~support_start).sum()),
            "unsupported_end": int((~support_end).sum()),
            "unsupported_pair": int((~valid).sum()),
        })
        metadata.append(row)
        del grid
    support = np.isfinite(log_start) & np.isfinite(log_end)
    bounds = {
        "lower": lower.tolist(), "upper": upper.tolist(),
        "grid_shape": list(GRID_SHAPE),
        "supported_pairs": int(support.sum()),
        "unsupported_pairs": int((~support).sum()),
    }
    return log_start, log_end, metadata, bounds


def driven_kde_stability(dataset: Dataset) -> tuple[list[dict], dict, list[dict]]:
    start_z = transform(dataset.start)
    end_z = transform(dataset.end)
    lower, upper = global_grid_bounds(end_z, dataset.stream)
    audit_mask = dataset.stream % 3 == 2
    positions = np.flatnonzero(audit_mask)
    estimates = []
    support_masks = []
    metadata = []
    for training_group in (0, 1):
        train = dataset.stream % 3 == training_group
        grid, row = fit_grid(
            end_z[train], dataset.stream[train], lower, upper,
            f"driven_stability_train_mod3_{training_group}")
        density_start, support_start = interpolate_density(
            grid, start_z[audit_mask], lower, upper)
        density_end, support_end = interpolate_density(
            grid, end_z[audit_mask], lower, upper)
        valid = support_start & support_end
        ls = np.full(audit_mask.sum(), np.nan, dtype=np.float64)
        le = np.full(audit_mask.sum(), np.nan, dtype=np.float64)
        ls[valid] = (np.log(density_start[valid])
                     - np.sum(start_z[audit_mask][valid, :3], axis=1))
        le[valid] = (np.log(density_end[valid])
                     - np.sum(end_z[audit_mask][valid, :3], axis=1))
        estimates.append((ls, le))
        support_masks.append(valid)
        row.update({"audit_blocks": int(audit_mask.sum()),
                    "unsupported_audit_pairs": int((~valid).sum())})
        metadata.append(row)
        del grid

    common = support_masks[0] & support_masks[1]
    if not np.any(common):
        raise RuntimeError("no common support in driven KDE stability audit")
    endpoint_difference = estimates[0][1][common] - estimates[1][1][common]
    endpoint_difference -= np.mean(endpoint_difference)
    delta_a = -estimates[0][1][common] + estimates[0][0][common]
    delta_b = -estimates[1][1][common] + estimates[1][0][common]
    increment_difference = delta_a - delta_b

    def metrics(values: np.ndarray) -> dict:
        absolute = np.abs(values)
        return {
            "n": int(values.size), "mean": float(np.mean(values)),
            "rmse": float(np.sqrt(np.mean(values * values))),
            "mae": float(np.mean(absolute)),
            "q90_abs": float(np.quantile(absolute, 0.90)),
            "q95_abs": float(np.quantile(absolute, 0.95)),
            "q99_abs": float(np.quantile(absolute, 0.99)),
            "max_abs": float(np.max(absolute)),
        }

    rows = []
    endpoint_metrics = metrics(endpoint_difference)
    increment_metrics = metrics(increment_difference)
    for quantity, values in (
        ("centered_endpoint_log_density_disagreement", endpoint_difference),
        ("system_entropy_increment_disagreement", increment_difference),
    ):
        row = {"quantity": quantity, "region": "all"}
        row.update(metrics(values))
        rows.append(row)

    average_log_start = 0.5 * (estimates[0][0][common] + estimates[1][0][common])
    average_log_end = 0.5 * (estimates[0][1][common] + estimates[1][1][common])
    lower_endpoint_density = np.minimum(average_log_start, average_log_end)
    cutoff = float(np.quantile(lower_endpoint_density, 0.01))
    tail = lower_endpoint_density <= cutoff
    tail_metrics = metrics(increment_difference[tail])
    tail_row = {"quantity": "system_entropy_increment_disagreement",
                "region": "lowest_density_1_percent",
                "log_density_cutoff": cutoff}
    tail_row.update(tail_metrics)
    rows.append(tail_row)
    passed = bool(
        np.all(common)
        and endpoint_metrics["rmse"] <= KDE_ENDPOINT_RMSE_MAX
        and increment_metrics["rmse"] <= KDE_INCREMENT_RMSE_MAX
        and tail_metrics["rmse"] <= KDE_TAIL_INCREMENT_RMSE_MAX
        and increment_metrics["q99_abs"] <= KDE_INCREMENT_Q99_MAX
    )
    audit = {
        "audit_stream_group": "stream_id mod 3 == 2",
        "audit_blocks": int(positions.size),
        "common_supported_blocks": int(common.sum()),
        "unsupported_by_either_estimator": int((~common).sum()),
        "endpoint_disagreement_rmse": endpoint_metrics["rmse"],
        "increment_disagreement_rmse": increment_metrics["rmse"],
        "tail_increment_disagreement_rmse": tail_metrics["rmse"],
        "increment_disagreement_q99_abs": increment_metrics["q99_abs"],
        "pass": passed,
    }
    return rows, audit, metadata


def summarize_errors(label: str, temperature: float, dataset: Dataset,
                     log_start: np.ndarray, log_end: np.ndarray) -> tuple[list[dict], list[dict], dict]:
    supported = np.isfinite(log_start) & np.isfinite(log_end)
    if not np.any(supported):
        raise RuntimeError(f"no supported equilibrium endpoints for {label}")
    start_energy = energy(dataset.start)
    end_energy = energy(dataset.end)
    endpoint_error = np.full(EXPECTED_ROWS, np.nan, dtype=np.float64)
    fold_offsets = {}
    for parity in (0, 1):
        mask = (dataset.stream % 2 == parity) & supported
        raw = log_end[mask] + end_energy[mask] / temperature
        offset = float(np.mean(raw))
        endpoint_error[mask] = raw - offset
        fold_offsets[str(parity)] = offset
    exact_increment = (end_energy - start_energy) / temperature
    kde_increment = -log_end + log_start
    increment_error = kde_increment - exact_increment

    def metrics(values: np.ndarray) -> dict:
        absolute = np.abs(values)
        return {
            "n": int(values.size),
            "mean": float(np.mean(values)),
            "rmse": float(np.sqrt(np.mean(values * values))),
            "mae": float(np.mean(absolute)),
            "median_abs": float(np.quantile(absolute, 0.5)),
            "q90_abs": float(np.quantile(absolute, 0.9)),
            "q95_abs": float(np.quantile(absolute, 0.95)),
            "q99_abs": float(np.quantile(absolute, 0.99)),
            "max_abs": float(np.max(absolute)),
        }

    overall_rows = []
    for quantity, values in (("centered_endpoint_log_density", endpoint_error),
                             ("system_entropy_increment", increment_error)):
        finite = np.isfinite(values)
        row = {"dataset": label, "temperature": temperature,
               "quantity": quantity, "region": "all",
               "total_n": int(values.size),
               "unsupported_n": int((~finite).sum())}
        row.update(metrics(values[finite]))
        overall_rows.append(row)
        for parity in (0, 1):
            mask = (dataset.stream % 2 == parity) & finite
            fold_row = {"dataset": label, "temperature": temperature,
                        "quantity": quantity,
                        "region": f"evaluation_parity_{parity}",
                        "total_n": int(np.count_nonzero(dataset.stream % 2 == parity)),
                        "unsupported_n": int(np.count_nonzero(
                            (dataset.stream % 2 == parity) & ~finite))}
            fold_row.update(metrics(values[mask]))
            overall_rows.append(fold_row)

    maximum_energy = np.maximum(start_energy, end_energy)
    cut = np.quantile(maximum_energy, [0.0, 0.8, 0.95, 0.99, 1.0])
    tail_rows = []
    labels = ["0-80", "80-95", "95-99", "99-100"]
    for index, region in enumerate(labels):
        if index == 0:
            mask = (maximum_energy >= cut[index]) & (maximum_energy <= cut[index + 1])
        else:
            mask = (maximum_energy > cut[index]) & (maximum_energy <= cut[index + 1])
        supported_mask = mask & supported
        row = {"dataset": label, "temperature": temperature,
               "quantity": "system_entropy_increment", "region": region,
               "total_n": int(mask.sum()),
               "unsupported_n": int((mask & ~supported).sum()),
               "energy_low": float(cut[index]),
               "energy_high": float(cut[index + 1])}
        row.update(metrics(increment_error[supported_mask]))
        tail_rows.append(row)

    endpoint_rmse = metrics(endpoint_error[supported])["rmse"]
    increment_metrics = metrics(increment_error[supported])
    tail_rmse = tail_rows[-1]["rmse"]
    passed = (
        np.all(supported)
        and endpoint_rmse <= KDE_ENDPOINT_RMSE_MAX
        and increment_metrics["rmse"] <= KDE_INCREMENT_RMSE_MAX
        and tail_rmse <= KDE_TAIL_INCREMENT_RMSE_MAX
        and increment_metrics["q99_abs"] <= KDE_INCREMENT_Q99_MAX
    )
    audit = {
        "dataset": label, "temperature": temperature,
        "fold_additive_offsets": fold_offsets,
        "supported_pairs": int(supported.sum()),
        "unsupported_pairs": int((~supported).sum()),
        "support_gate_pass": bool(np.all(supported)),
        "endpoint_rmse": endpoint_rmse,
        "increment_rmse": increment_metrics["rmse"],
        "tail_99_100_increment_rmse": tail_rmse,
        "increment_q99_abs": increment_metrics["q99_abs"],
        "thresholds": {
            "endpoint_rmse_max": KDE_ENDPOINT_RMSE_MAX,
            "increment_rmse_max": KDE_INCREMENT_RMSE_MAX,
            "tail_increment_rmse_max": KDE_TAIL_INCREMENT_RMSE_MAX,
            "increment_q99_abs_max": KDE_INCREMENT_Q99_MAX,
        },
        "pass": bool(passed),
    }
    return overall_rows, tail_rows, audit


def weighted_fit(x: np.ndarray, y: np.ndarray, positive: np.ndarray,
                 negative: np.ndarray) -> tuple[float, float, float]:
    weights = 1.0 / (1.0 / positive + 1.0 / negative)
    design = np.column_stack([np.ones_like(x), x])
    normal = design.T @ (weights[:, None] * design)
    beta = np.linalg.solve(normal, design.T @ (weights * y))
    fitted = design @ beta
    mean = float(np.average(y, weights=weights))
    ss_res = float(np.sum(weights * (y - fitted) ** 2))
    ss_tot = float(np.sum(weights * (y - mean) ** 2))
    r2 = float("nan") if ss_tot == 0.0 else 1.0 - ss_res / ss_tot
    return float(beta[0]), float(beta[1]), r2


def observable_diagnostics(label: str, values: np.ndarray,
                           stream: np.ndarray, rng: np.random.Generator):
    if not np.isfinite(values).all():
        raise RuntimeError(f"nonfinite {label} value")
    shaped = values.reshape(EXPECTED_STREAMS, EXPECTED_BLOCKS_PER_STREAM)
    negative = values < 0.0
    positive = values > 0.0
    zero = values == 0.0
    stream_negative = negative.reshape(EXPECTED_STREAMS,
                                       EXPECTED_BLOCKS_PER_STREAM).mean(axis=1)
    draws = rng.integers(0, EXPECTED_STREAMS,
                         size=(BOOTSTRAPS, EXPECTED_STREAMS))
    probability_draws = stream_negative[draws].mean(axis=1)
    count_row = {
        "observable": label,
        "n_blocks": int(values.size),
        "negative_count": int(negative.sum()),
        "positive_count": int(positive.sum()),
        "zero_count": int(zero.sum()),
        "negative_probability": float(negative.mean()),
        "negative_probability_ci_low": float(np.quantile(probability_draws, 0.025)),
        "negative_probability_ci_high": float(np.quantile(probability_draws, 0.975)),
    }

    q25, q75 = np.quantile(values, [0.25, 0.75])
    width = 2.0 * float(q75 - q25) / np.cbrt(values.size)
    if not (width > 0.0 and np.isfinite(width)):
        raise RuntimeError(f"invalid FD width for {label}")
    indices = np.floor(np.abs(values) / width).astype(np.int64)
    n_bins = int(indices.max()) + 1
    positive_by_stream = np.zeros((EXPECTED_STREAMS, n_bins), dtype=np.int64)
    negative_by_stream = np.zeros((EXPECTED_STREAMS, n_bins), dtype=np.int64)
    np.add.at(positive_by_stream, (stream[positive], indices[positive]), 1)
    np.add.at(negative_by_stream, (stream[negative], indices[negative]), 1)
    n_plus = positive_by_stream.sum(axis=0)
    n_minus = negative_by_stream.sum(axis=0)
    support = (n_plus >= MIN_PAIR_COUNT) & (n_minus >= MIN_PAIR_COUNT)
    centers = (np.arange(n_bins, dtype=np.float64) + 0.5) * width
    bin_rows = []
    for index in range(n_bins):
        bin_rows.append({
            "observable": label, "bin_index": index,
            "abs_left": index * width, "abs_right": (index + 1) * width,
            "a_center": centers[index],
            "positive_count": int(n_plus[index]),
            "negative_count": int(n_minus[index]),
            "qualifies_ge_20_each": int(support[index]),
            "log_count_ratio": (float(math.log(n_plus[index] / n_minus[index]))
                                if n_plus[index] > 0 and n_minus[index] > 0
                                else float("nan")),
        })

    fit_row = {
        "observable": label, "fd_width": width,
        "qualifying_pairs": int(support.sum()),
        "support_gate_negative_min": MIN_NEGATIVE,
        "support_gate_pairs_min": MIN_PAIRS,
        "support_gate_pass": int(negative.sum() >= MIN_NEGATIVE
                                 and support.sum() >= MIN_PAIRS),
        "fit_available": 0,
        "intercept": float("nan"),
        "intercept_ci_low": float("nan"),
        "intercept_ci_high": float("nan"),
        "slope": float("nan"),
        "slope_ci_low": float("nan"),
        "slope_ci_high": float("nan"),
        "weighted_r_squared": float("nan"),
        "bootstrap_valid": 0,
    }
    if support.sum() >= 2:
        x = centers[support]
        p = n_plus[support].astype(np.float64)
        m = n_minus[support].astype(np.float64)
        intercept, slope, r2 = weighted_fit(x, np.log(p / m), p, m)
        boot_intercept = []
        boot_slope = []
        for draw in draws:
            bp = positive_by_stream[draw].sum(axis=0)[support]
            bm = negative_by_stream[draw].sum(axis=0)[support]
            if np.any(bp == 0) or np.any(bm == 0):
                continue
            bi, bs, _ = weighted_fit(x, np.log(bp / bm),
                                     bp.astype(float), bm.astype(float))
            boot_intercept.append(bi)
            boot_slope.append(bs)
        if len(boot_slope) >= int(0.95 * BOOTSTRAPS):
            fit_row.update({
                "fit_available": 1, "intercept": intercept,
                "intercept_ci_low": float(np.quantile(boot_intercept, 0.025)),
                "intercept_ci_high": float(np.quantile(boot_intercept, 0.975)),
                "slope": slope,
                "slope_ci_low": float(np.quantile(boot_slope, 0.025)),
                "slope_ci_high": float(np.quantile(boot_slope, 0.975)),
                "weighted_r_squared": r2,
                "bootstrap_valid": len(boot_slope),
            })

    exponent = -values
    maximum = float(np.max(exponent))
    scaled = np.exp(exponent - maximum)
    total = float(np.sum(scaled))
    log_ift = maximum + math.log(total / values.size)
    ess = total * total / float(np.sum(scaled * scaled))
    max_share = float(np.max(scaled) / total)
    stream_sum = scaled.reshape(EXPECTED_STREAMS,
                                EXPECTED_BLOCKS_PER_STREAM).sum(axis=1)
    boot_sum = stream_sum[draws].sum(axis=1)
    boot_log = maximum + np.log(boot_sum / values.size)
    jackknife = np.empty(EXPECTED_STREAMS, dtype=np.float64)
    for removed in range(EXPECTED_STREAMS):
        remaining_sum = total - stream_sum[removed]
        remaining_n = values.size - EXPECTED_BLOCKS_PER_STREAM
        jackknife[removed] = maximum + math.log(remaining_sum / remaining_n)
    jackknife_change = float(np.max(np.abs(jackknife - log_ift)))
    ift_row = {
        "observable": label, "log_mean_exp_minus": log_ift,
        "bootstrap_ci_low": float(np.quantile(boot_log, 0.025)),
        "bootstrap_ci_high": float(np.quantile(boot_log, 0.975)),
        "exponential_weight_ess": ess,
        "ess_fraction": ess / values.size,
        "maximum_single_weight_share": max_share,
        "maximum_leave_one_stream_change": jackknife_change,
        "resolution_pass": int(ess >= IFT_ESS_MIN
                               and max_share <= IFT_MAX_SHARE_MAX
                               and jackknife_change <= IFT_JACKKNIFE_CHANGE_MAX),
    }
    return count_row, fit_row, bin_rows, ift_row


def write_derived(path: Path, dataset: Dataset, log_start: np.ndarray,
                  log_end: np.ndarray, delta_sys: np.ndarray,
                  delta_total: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    process = subprocess.Popen(["zstd", "-T0", "-12", "-q", "-o", str(path)],
                               stdin=subprocess.PIPE, text=True)
    assert process.stdin is not None
    handle = process.stdin
    handle.write("stream_id,block_id,entropy_medium,log_rho_start,log_rho_end,delta_s_sys,delta_s_tot\n")
    for index in range(EXPECTED_ROWS):
        handle.write(
            f"{int(dataset.stream[index])},{int(dataset.block[index])},"
            f"{dataset.entropy_medium[index]:.17g},{log_start[index]:.17g},"
            f"{log_end[index]:.17g},{delta_sys[index]:.17g},"
            f"{delta_total[index]:.17g}\n")
    handle.close()
    if process.wait() != 0:
        raise RuntimeError("failed to compress derived block output")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--driven", type=Path, required=True)
    parser.add_argument("--equilibrium-t6", type=Path, required=True)
    parser.add_argument("--equilibrium-t10", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--write-derived", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    input_rows = []
    for label, path in (("driven", args.driven), ("T6", args.equilibrium_t6),
                        ("T10", args.equilibrium_t10)):
        subprocess.run(["zstd", "-t", "-q", str(path)], check=True)
        decompressed_hash = sha256_decompressed(path)
        if decompressed_hash != EXPECTED_DECOMPRESSED_HASHES[label]:
            raise RuntimeError(
                f"{label} decompressed SHA-256 mismatch: {decompressed_hash} != "
                f"{EXPECTED_DECOMPRESSED_HASHES[label]}")
        input_rows.append({
            "label": label, "path": str(path.resolve()),
            "compressed_sha256": sha256(path),
            "decompressed_sha256": decompressed_hash,
            "expected_decompressed_sha256": EXPECTED_DECOMPRESSED_HASHES[label],
            "hash_match": 1,
            "compressed_bytes": path.stat().st_size,
        })
    write_csv(args.output / "input_hashes.csv", input_rows)

    kde_rows = []
    accuracy_rows = []
    tail_rows = []
    accuracy_audits = []
    for label, temperature, archive in (
        ("T6", 6.0, args.equilibrium_t6),
        ("T10", 10.0, args.equilibrium_t10),
    ):
        dataset = read_archive(archive)
        log_start, log_end, metadata, bounds = crossfit_density(dataset, label)
        kde_rows.extend(metadata)
        rows, tails, audit = summarize_errors(
            label, temperature, dataset, log_start, log_end)
        accuracy_rows.extend(rows)
        tail_rows.extend(tails)
        audit["bounds"] = bounds
        accuracy_audits.append(audit)
        del dataset, log_start, log_end

    write_csv(args.output / "kde_bandwidths.csv", kde_rows)
    write_csv(args.output / "equilibrium_kde_accuracy.csv", accuracy_rows)
    write_csv(args.output / "equilibrium_kde_tail_accuracy.csv", tail_rows)
    kde_gate_pass = all(row["pass"] for row in accuracy_audits)
    with (args.output / "equilibrium_kde_gate.json").open("w") as handle:
        json.dump({"all_pass": kde_gate_pass, "cases": accuracy_audits},
                  handle, indent=2)
        handle.write("\n")

    driven = read_archive(args.driven)
    stability_rows, stability_audit, stability_metadata = driven_kde_stability(driven)
    kde_rows.extend(stability_metadata)
    write_csv(args.output / "driven_kde_stability.csv", stability_rows)
    with (args.output / "driven_kde_stability_gate.json").open("w") as handle:
        json.dump(stability_audit, handle, indent=2)
        handle.write("\n")
    log_start, log_end, metadata, driven_bounds = crossfit_density(driven, "driven")
    kde_rows.extend(metadata)
    write_csv(args.output / "kde_bandwidths.csv", kde_rows)
    supported = np.isfinite(log_start) & np.isfinite(log_end)
    support_row = {
        "n_blocks": EXPECTED_ROWS, "supported_pairs": int(supported.sum()),
        "unsupported_pairs": int((~supported).sum()),
        "supported_fraction": float(supported.mean()),
    }
    write_csv(args.output / "driven_kde_support.csv", [support_row])
    delta_sys = -log_end + log_start
    delta_total = driven.entropy_medium + delta_sys
    system_rows = []
    for region, mask in (
        ("all", np.ones(EXPECTED_ROWS, dtype=bool)),
        ("evaluation_parity_0", driven.stream % 2 == 0),
        ("evaluation_parity_1", driven.stream % 2 == 1),
    ):
        total_count = int(mask.sum())
        valid = mask & supported
        values = delta_sys[valid]
        if values.size == 0:
            continue
        system_rows.append({
            "region": region, "n": int(values.size),
            "total_n": total_count,
            "unsupported_n": total_count - int(values.size),
            "mean": float(np.mean(values)),
            "std": float(np.std(values, ddof=1)),
            "minimum": float(np.min(values)),
            "q01": float(np.quantile(values, 0.01)),
            "q05": float(np.quantile(values, 0.05)),
            "median": float(np.quantile(values, 0.5)),
            "q95": float(np.quantile(values, 0.95)),
            "q99": float(np.quantile(values, 0.99)),
            "maximum": float(np.max(values)),
        })
    write_csv(args.output / "driven_system_entropy_summary.csv", system_rows)
    if args.write_derived:
        write_derived(args.output / "derived_blocks.csv.zst", driven,
                      log_start, log_end, delta_sys, delta_total)

    rng = np.random.default_rng(BOOTSTRAP_SEED)
    count_rows = []
    fit_rows = []
    bin_rows = []
    ift_rows = []
    count, fit, bins, ift = observable_diagnostics(
        "medium_entropy", driven.entropy_medium, driven.stream, rng)
    count["analysis_status"] = "computed"
    count["kde_unsupported_pairs"] = 0
    fit["analysis_status"] = "computed"
    ift["analysis_status"] = "computed"
    count_rows.append(count)
    fit_rows.append(fit)
    bin_rows.extend(bins)
    ift_rows.append(ift)

    if np.all(supported):
        count, fit, bins, ift = observable_diagnostics(
            "total_entropy_kde", delta_total, driven.stream, rng)
        count["analysis_status"] = "computed"
        count["kde_unsupported_pairs"] = 0
        fit["analysis_status"] = "computed"
        ift["analysis_status"] = "computed"
        count_rows.append(count)
        fit_rows.append(fit)
        bin_rows.extend(bins)
        ift_rows.append(ift)
    else:
        count_rows.append({
            "observable": "total_entropy_kde", "n_blocks": EXPECTED_ROWS,
            "negative_count": float("nan"), "positive_count": float("nan"),
            "zero_count": float("nan"), "negative_probability": float("nan"),
            "negative_probability_ci_low": float("nan"),
            "negative_probability_ci_high": float("nan"),
            "analysis_status": "not_computed_no_kde_extrapolation",
            "kde_unsupported_pairs": int((~supported).sum()),
        })
        fit_rows.append({
            "observable": "total_entropy_kde", "fd_width": float("nan"),
            "qualifying_pairs": 0, "support_gate_negative_min": MIN_NEGATIVE,
            "support_gate_pairs_min": MIN_PAIRS, "support_gate_pass": 0,
            "fit_available": 0, "intercept": float("nan"),
            "intercept_ci_low": float("nan"),
            "intercept_ci_high": float("nan"), "slope": float("nan"),
            "slope_ci_low": float("nan"), "slope_ci_high": float("nan"),
            "weighted_r_squared": float("nan"), "bootstrap_valid": 0,
            "analysis_status": "not_computed_no_kde_extrapolation",
        })
        ift_rows.append({
            "observable": "total_entropy_kde",
            "log_mean_exp_minus": float("nan"),
            "bootstrap_ci_low": float("nan"),
            "bootstrap_ci_high": float("nan"),
            "exponential_weight_ess": float("nan"),
            "ess_fraction": float("nan"),
            "maximum_single_weight_share": float("nan"),
            "maximum_leave_one_stream_change": float("nan"),
            "resolution_pass": 0,
            "analysis_status": "not_computed_no_kde_extrapolation",
        })
    write_csv(args.output / "negative_tail_counts.csv", count_rows)
    write_csv(args.output / "detailed_ft_fits.csv", fit_rows)
    write_csv(args.output / "symmetric_bin_counts.csv", bin_rows)
    write_csv(args.output / "integral_ft.csv", ift_rows)

    total_fit = next(row for row in fit_rows
                     if row["observable"] == "total_entropy_kde")
    total_ift = next(row for row in ift_rows
                     if row["observable"] == "total_entropy_kde")
    detailed_reference_pass = bool(
        total_fit.get("fit_available", 0)
        and total_fit["slope_ci_low"] <= 1.0 <= total_fit["slope_ci_high"]
    )
    integral_reference_pass = bool(
        total_ift["bootstrap_ci_low"] <= 0.0 <= total_ift["bootstrap_ci_high"]
    )
    final_validated = bool(
        kde_gate_pass
        and stability_audit["pass"]
        and total_fit["support_gate_pass"]
        and total_fit.get("fit_available", 0)
        and detailed_reference_pass
        and total_ift["resolution_pass"]
        and integral_reference_pass
    )
    audit = {
        "grid_shape": list(GRID_SHAPE), "kernel": "product Gaussian",
        "angle_boundary": "wrapped periodic", "action_transform": "log",
        "action_log_jacobian_included": True,
        "crossfit": "stream parity two-fold",
        "bootstrap_replicates": BOOTSTRAPS,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "equilibrium_kde_gate_pass": kde_gate_pass,
        "driven_kde_stability_gate_pass": bool(stability_audit["pass"]),
        "driven_bounds": driven_bounds,
        "driven_support": support_row,
        "total_two_sided_support_pass": bool(total_fit["support_gate_pass"]),
        "detailed_reference_slope_ci_contains_one": detailed_reference_pass,
        "integral_resolution_pass": bool(total_ift["resolution_pass"]),
        "integral_reference_ci_contains_zero": integral_reference_pass,
        "total_entropy_ft_validated": final_validated,
        "claim": (
            "validated finite-time n=3 total-entropy FT diagnostic"
            if final_validated else
            "n=3 KDE total-entropy FT not reliably resolved under frozen gates"
        ),
    }
    with (args.output / "FINAL_AUDIT.json").open("w") as handle:
        json.dump(audit, handle, indent=2)
        handle.write("\n")
    print(json.dumps({"event": "analysis_complete", **audit}), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
