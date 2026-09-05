#!/usr/bin/env python3
"""Frozen integrity and decorrelation audit for n=3 NEEP trajectory data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from pathlib import Path

import numpy as np


RAW_HEADER = [
    "stream_id", "block_id", "q_left", "q_right", "delta_energy",
    "entropy_medium", "entropy_rate", "action_current",
    "energy_balance_error", "start_I1", "start_I2", "start_I3",
    "start_theta1", "start_theta3", "end_I1", "end_I2", "end_I3",
    "end_theta1", "end_theta3",
]
TRIG_HEADER = [
    "start_theta1_cos", "start_theta1_sin",
    "start_theta3_cos", "start_theta3_sin",
    "end_theta1_cos", "end_theta1_sin",
    "end_theta3_cos", "end_theta3_sin",
]
STREAMS = 128
TRANSITIONS_PER_STREAM = 39063
ROWS = STREAMS * TRANSITIONS_PER_STREAM
SNAPSHOTS = STREAMS * (TRANSITIONS_PER_STREAM + 1)
DELTA_T = 0.1
MAX_LAG = 10000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stream_metadata(path: Path) -> tuple[str, int, bytes]:
    process = subprocess.Popen(["zstd", "-dc", str(path)],
                               stdout=subprocess.PIPE)
    assert process.stdout is not None
    digest = hashlib.sha256()
    line_count = 0
    first = b""
    carry = b""
    for chunk in iter(lambda: process.stdout.read(4 * 1024 * 1024), b""):
        digest.update(chunk)
        if not first:
            combined = carry + chunk
            index = combined.find(b"\n")
            if index >= 0:
                first = combined[:index]
                carry = b""
            else:
                carry = combined
        line_count += chunk.count(b"\n")
    process.stdout.close()
    if process.wait() != 0:
        raise RuntimeError(f"zstd decompression failed: {path}")
    return digest.hexdigest(), line_count, first


def read_raw(path: Path) -> np.ndarray:
    process = subprocess.Popen(["zstd", "-dc", str(path)],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               text=True, bufsize=4 * 1024 * 1024)
    assert process.stdout is not None
    header = next(csv.reader([process.stdout.readline().rstrip("\n")]))
    if header != RAW_HEADER:
        raise RuntimeError(f"unexpected raw header: {header}")
    matrix = np.loadtxt(process.stdout, delimiter=",", dtype=np.float64)
    process.stdout.close()
    stderr = process.stderr.read() if process.stderr is not None else ""
    if process.wait() != 0:
        raise RuntimeError(f"zstd failed for {path}: {stderr}")
    if matrix.shape != (ROWS, len(RAW_HEADER)):
        raise RuntimeError(f"unexpected raw shape {matrix.shape}")
    if not np.isfinite(matrix).all():
        raise RuntimeError("nonfinite raw value")
    return matrix


def verify_neep_values(path: Path) -> dict:
    process = subprocess.Popen(["zstd", "-dc", str(path)],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               text=True, bufsize=4 * 1024 * 1024)
    assert process.stdout is not None
    header = next(csv.reader([process.stdout.readline().rstrip("\n")]))
    if header != RAW_HEADER + TRIG_HEADER:
        raise RuntimeError("NEEP header mismatch during value audit")
    columns = (0, 1, 12, 13, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26)
    values = np.loadtxt(process.stdout, delimiter=",", dtype=np.float64,
                        usecols=columns)
    process.stdout.close()
    stderr = process.stderr.read() if process.stderr is not None else ""
    if process.wait() != 0:
        raise RuntimeError(f"zstd failed for {path}: {stderr}")
    if values.shape != (ROWS, len(columns)):
        raise RuntimeError(f"unexpected NEEP value shape {values.shape}")
    if not np.isfinite(values).all():
        raise RuntimeError("nonfinite NEEP value")
    streams = values[:, 0].astype(np.int16)
    intervals = values[:, 1].astype(np.int32)
    if not np.array_equal(streams, np.repeat(np.arange(STREAMS),
                                              TRANSITIONS_PER_STREAM)):
        raise RuntimeError("NEEP stream ordering failure")
    if not np.array_equal(intervals, np.tile(np.arange(TRANSITIONS_PER_STREAM),
                                              STREAMS)):
        raise RuntimeError("NEEP interval ordering failure")
    angles = values[:, 2:6]
    encoded = values[:, 6:14]
    expected = np.column_stack((
        np.cos(angles[:, 0]), np.sin(angles[:, 0]),
        np.cos(angles[:, 1]), np.sin(angles[:, 1]),
        np.cos(angles[:, 2]), np.sin(angles[:, 2]),
        np.cos(angles[:, 3]), np.sin(angles[:, 3]),
    ))
    encoding_error = float(np.max(np.abs(encoded - expected)))
    unit_circle_error = float(np.max(np.abs(
        encoded[:, 0::2] * encoded[:, 0::2]
        + encoded[:, 1::2] * encoded[:, 1::2] - 1.0)))
    return {
        "trigonometric_encoding_max_abs_error": encoding_error,
        "unit_circle_identity_max_abs_error": unit_circle_error,
    }


def autocovariance_matrix(values: np.ndarray) -> np.ndarray:
    centered = values - values.mean(axis=1, keepdims=True)
    length = centered.shape[1]
    nfft = 1 << (2 * length - 1).bit_length()
    spectrum = np.fft.rfft(centered, n=nfft, axis=1)
    autocov = np.fft.irfft(spectrum * np.conjugate(spectrum),
                           n=nfft, axis=1)[:, :MAX_LAG + 1].real
    autocov /= np.arange(length, length - MAX_LAG - 1, -1,
                         dtype=np.float64)[None, :]
    return autocov


def first_crossing(rho: np.ndarray, level: float = math.exp(-1.0)) -> float:
    below = np.flatnonzero(rho <= level)
    if below.size == 0:
        return float("nan")
    index = int(below[0])
    if index == 0:
        return 0.0
    y0 = rho[index - 1]
    y1 = rho[index]
    if y1 == y0:
        return float(index)
    return (index - 1) + (level - y0) / (y1 - y0)


def statistical_inefficiency(rho: np.ndarray) -> tuple[float, int, bool]:
    cutoff = 1
    stopped = False
    pair = 0
    while 2 * pair + 2 < rho.size:
        left = 2 * pair + 1
        right = 2 * pair + 2
        if rho[left] + rho[right] <= 0.0:
            stopped = True
            break
        cutoff = right + 1
        pair += 1
    g = max(1.0, 1.0 + 2.0 * float(np.sum(rho[1:cutoff])))
    return g, cutoff, stopped


def correlation_summary(values: np.ndarray, periodic: bool) -> tuple[dict, np.ndarray]:
    if periodic:
        autocov = (autocovariance_matrix(np.cos(values))
                   + autocovariance_matrix(np.sin(values)))
    else:
        autocov = autocovariance_matrix(values)
    if np.any(autocov[:, 0] <= 0.0):
        raise RuntimeError("nonpositive stream variance")
    pooled = autocov.mean(axis=0)
    rho = pooled / pooled[0]
    stream_rho = autocov / autocov[:, :1]
    g, cutoff, stopped = statistical_inefficiency(rho)
    stream_g = np.array([statistical_inefficiency(row)[0]
                         for row in stream_rho])
    stream_cross = np.array([first_crossing(row) for row in stream_rho])
    cross = first_crossing(rho)
    result = {
        "lag1_correlation": float(rho[1]),
        "one_over_e_lag": float(cross),
        "one_over_e_time": float(cross * DELTA_T),
        "statistical_inefficiency_g": float(g),
        "tau_integrated": float(0.5 * g * DELTA_T),
        "ips_cutoff_lag": int(cutoff),
        "ips_found_nonpositive_pair": bool(stopped),
        "delta_t_over_one_over_e_time": float(1.0 / cross),
        "stream_one_over_e_time_q05": float(np.nanquantile(stream_cross * DELTA_T, 0.05)),
        "stream_one_over_e_time_median": float(np.nanmedian(stream_cross * DELTA_T)),
        "stream_one_over_e_time_q95": float(np.nanquantile(stream_cross * DELTA_T, 0.95)),
        "stream_tau_integrated_q05": float(np.quantile(0.5 * stream_g * DELTA_T, 0.05)),
        "stream_tau_integrated_median": float(np.median(0.5 * stream_g * DELTA_T)),
        "stream_tau_integrated_q95": float(np.quantile(0.5 * stream_g * DELTA_T, 0.95)),
    }
    return result, rho


def verify_neep_header(path: Path) -> dict:
    decompressed_sha, lines, first = stream_metadata(path)
    header = next(csv.reader([first.decode("utf-8")]))
    if header != RAW_HEADER + TRIG_HEADER:
        raise RuntimeError(f"unexpected NEEP header: {header}")
    if lines != ROWS + 1:
        raise RuntimeError(f"unexpected NEEP line count: {lines}")
    return {
        "path": str(path),
        "compressed_bytes": path.stat().st_size,
        "compressed_sha256": sha256(path),
        "decompressed_sha256": decompressed_sha,
        "line_count_including_header": lines,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", required=True, type=Path)
    parser.add_argument("--neep", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--tl", required=True, type=float)
    parser.add_argument("--tr", required=True, type=float)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    raw_meta = verify_neep_header(args.neep)
    raw_sha, raw_lines, raw_first = stream_metadata(args.raw)
    if next(csv.reader([raw_first.decode("utf-8")])) != RAW_HEADER:
        raise RuntimeError("raw header mismatch")
    if raw_lines != ROWS + 1:
        raise RuntimeError(f"unexpected raw line count: {raw_lines}")
    matrix = read_raw(args.raw)

    stream = matrix[:, 0].astype(np.int16)
    interval = matrix[:, 1].astype(np.int32)
    if not np.array_equal(stream, np.repeat(np.arange(STREAMS),
                                             TRANSITIONS_PER_STREAM)):
        raise RuntimeError("stream ordering failure")
    if not np.array_equal(interval, np.tile(np.arange(TRANSITIONS_PER_STREAM),
                                             STREAMS)):
        raise RuntimeError("interval ordering failure")

    start = matrix[:, 9:14].reshape(STREAMS, TRANSITIONS_PER_STREAM, 5)
    end = matrix[:, 14:19].reshape(STREAMS, TRANSITIONS_PER_STREAM, 5)
    if np.min(start[:, :, :3]) <= 0.0 or np.min(end[:, :, :3]) <= 0.0:
        raise RuntimeError("nonpositive action")
    if (np.min(start[:, :, 3:]) < -math.pi or
            np.max(start[:, :, 3:]) >= math.pi or
            np.min(end[:, :, 3:]) < -math.pi or
            np.max(end[:, :, 3:]) >= math.pi):
        raise RuntimeError("angle outside [-pi,pi)")
    linear_continuity = float(np.max(np.abs(end[:, :-1, :3] -
                                                  start[:, 1:, :3])))
    angle_delta = (end[:, :-1, 3:] - start[:, 1:, 3:] + math.pi) % (
        2.0 * math.pi) - math.pi
    angular_continuity = float(np.max(np.abs(angle_delta)))
    if linear_continuity > 1e-10 or angular_continuity > 1e-10:
        raise RuntimeError("endpoint continuity failure")

    residual = matrix[:, 8]
    recomputed_residual = matrix[:, 2] + matrix[:, 3] - matrix[:, 4]
    entropy_error = matrix[:, 5] + matrix[:, 2] / args.tl + matrix[:, 3] / args.tr
    state_i2 = np.concatenate((start[:, :1, 1], end[:, :, 1]), axis=1)
    state_theta1 = np.concatenate((start[:, :1, 3], end[:, :, 3]), axis=1)

    i2_summary, i2_rho = correlation_summary(state_i2, periodic=False)
    theta_summary, theta_rho = correlation_summary(state_theta1, periodic=True)
    quantiles = np.quantile(residual, [0.0, 0.001, 0.01, 0.5, 0.99, 0.999, 1.0])
    neep_value_audit = verify_neep_values(args.neep)

    summary = {
        "label": args.label,
        "temperatures": [args.tl, args.tr],
        "delta_t": DELTA_T,
        "integration_steps_per_transition": 200,
        "streams": STREAMS,
        "transitions_per_stream": TRANSITIONS_PER_STREAM,
        "transition_rows": ROWS,
        "state_snapshots": SNAPSHOTS,
        "raw_archive": {
            "path": str(args.raw),
            "compressed_bytes": args.raw.stat().st_size,
            "compressed_sha256": sha256(args.raw),
            "decompressed_sha256": raw_sha,
            "line_count_including_header": raw_lines,
        },
        "neep_archive": raw_meta,
        "integrity": {
            "all_values_finite": True,
            "positive_actions": True,
            "angles_in_minus_pi_pi": True,
            "linear_endpoint_continuity_max_abs": linear_continuity,
            "angular_endpoint_continuity_max_abs": angular_continuity,
            "stored_vs_recomputed_balance_max_abs": float(np.max(np.abs(
                residual - recomputed_residual))),
            "entropy_identity_max_abs": float(np.max(np.abs(entropy_error))),
            **neep_value_audit,
        },
        "first_law_residual": {
            "mean": float(np.mean(residual)),
            "std": float(np.std(residual, ddof=1)),
            "rms": float(np.sqrt(np.mean(residual * residual))),
            "rms_rate": float(np.sqrt(np.mean(residual * residual)) / DELTA_T),
            "min": float(quantiles[0]),
            "q001": float(quantiles[1]),
            "q01": float(quantiles[2]),
            "median": float(quantiles[3]),
            "q99": float(quantiles[4]),
            "q999": float(quantiles[5]),
            "max": float(quantiles[6]),
            "max_abs": float(np.max(np.abs(residual))),
        },
        "correlation": {"I2": i2_summary, "theta1_periodic": theta_summary},
    }

    with (args.output / f"{args.label}_sanity.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with (args.output / f"{args.label}_correlation.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["lag", "time", "rho_I2", "rho_theta1_periodic"])
        for lag, (rho_i2, rho_theta) in enumerate(zip(i2_rho, theta_rho)):
            writer.writerow([lag, lag * DELTA_T, rho_i2, rho_theta])
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
