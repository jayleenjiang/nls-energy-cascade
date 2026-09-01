#!/usr/bin/env python3
"""Predeclared direct-sampling feasibility and integrity audit for n=3."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np


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
BOOTSTRAP_REPLICATES = 2000
BOOTSTRAP_SEED = 2026083193


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stream_bootstrap_interval(indicator: np.ndarray, rng: np.random.Generator):
    per_stream = indicator.mean(axis=1)
    draws = rng.integers(0, len(per_stream),
                         size=(BOOTSTRAP_REPLICATES, len(per_stream)))
    values = per_stream[draws].mean(axis=1)
    return np.quantile(values, [0.025, 0.975])


def symmetric_counts(values: np.ndarray):
    q25, q75 = np.quantile(values, [0.25, 0.75])
    width = 2.0 * (q75 - q25) / np.cbrt(values.size)
    if not np.isfinite(width) or width <= 0.0:
        width = 3.5 * np.std(values, ddof=1) / np.cbrt(values.size)
    if not np.isfinite(width) or width <= 0.0:
        raise RuntimeError("cannot define a positive histogram width")
    indices = np.floor(np.abs(values) / width).astype(np.int64)
    size = int(indices.max()) + 1
    positive = np.bincount(indices[values >= 0.0], minlength=size)
    negative = np.bincount(indices[values < 0.0], minlength=size)
    return width, positive, negative


def weighted_line_fit(x: np.ndarray,
                      y: np.ndarray,
                      positive: np.ndarray,
                      negative: np.ndarray):
    weights = 1.0 / (1.0 / positive + 1.0 / negative)
    design = np.column_stack([np.ones_like(x), x])
    root_weight = np.sqrt(weights)
    coefficients, *_ = np.linalg.lstsq(
        design * root_weight[:, None], y * root_weight, rcond=None)
    fitted = design @ coefficients
    mean = np.average(y, weights=weights)
    ss_residual = np.sum(weights * (y - fitted) ** 2)
    ss_total = np.sum(weights * (y - mean) ** 2)
    r_squared = 1.0 - ss_residual / ss_total if ss_total > 0.0 else np.nan
    return float(coefficients[0]), float(coefficients[1]), float(r_squared)


def per_stream_symmetric_counts(values: np.ndarray,
                                width: float,
                                size: int):
    positive = np.zeros((values.shape[0], size), dtype=np.int64)
    negative = np.zeros((values.shape[0], size), dtype=np.int64)
    for stream in range(values.shape[0]):
        indices = np.floor(np.abs(values[stream]) / width).astype(np.int64)
        positive[stream] = np.bincount(
            indices[values[stream] >= 0.0], minlength=size)[:size]
        negative[stream] = np.bincount(
            indices[values[stream] < 0.0], minlength=size)[:size]
    return positive, negative


def log_mean_exponential_minus(values: np.ndarray):
    exponent = -values
    maximum = float(np.max(exponent))
    scaled = np.exp(exponent - maximum)
    total = float(np.sum(scaled))
    log_mean = maximum + math.log(total / values.size)
    ess = total * total / float(np.sum(scaled * scaled))
    maximum_share = float(np.max(scaled) / total)
    return log_mean, ess, maximum_share, maximum, scaled


def write_csv(path: Path, rows: list[dict], fields: list[str]):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("blocks", type=Path)
    parser.add_argument("summary", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    with args.blocks.open(newline="") as handle:
        header = next(csv.reader(handle))
    if header != EXPECTED_HEADER:
        raise RuntimeError(f"unexpected block header: {header}")

    data = np.loadtxt(args.blocks, delimiter=",", skiprows=1)
    if data.shape != (EXPECTED_ROWS, len(EXPECTED_HEADER)):
        raise RuntimeError(f"unexpected block matrix shape: {data.shape}")
    if not np.isfinite(data).all():
        raise RuntimeError("nonfinite value in production blocks")

    shaped = data.reshape(EXPECTED_STREAMS, EXPECTED_BLOCKS_PER_STREAM, -1)
    expected_stream = np.arange(EXPECTED_STREAMS)[:, None]
    expected_block = np.arange(EXPECTED_BLOCKS_PER_STREAM)[None, :]
    if not np.all(shaped[:, :, 0] == expected_stream):
        raise RuntimeError("stream ordering or identifiers are invalid")
    if not np.all(shaped[:, :, 1] == expected_block):
        raise RuntimeError("block ordering or identifiers are invalid")

    entropy_identity_error = np.max(np.abs(
        shaped[:, :, 5] + shaped[:, :, 2] / 10.0 + shaped[:, :, 3] / 2.0))
    balance_identity_error = np.max(np.abs(
        shaped[:, :, 8] -
        (shaped[:, :, 2] + shaped[:, :, 3] - shaped[:, :, 4])))

    linear_endpoint_error = np.max(np.abs(
        shaped[:, :-1, 14:17] - shaped[:, 1:, 9:12]))
    angle_delta = shaped[:, :-1, 17:19] - shaped[:, 1:, 12:14]
    angle_delta = (angle_delta + np.pi) % (2.0 * np.pi) - np.pi
    angular_endpoint_error = np.max(np.abs(angle_delta))
    action_minimum = float(np.min(shaped[:, :, [9, 10, 11, 14, 15, 16]]))
    angle_minimum = float(np.min(shaped[:, :, [12, 13, 17, 18]]))
    angle_maximum = float(np.max(shaped[:, :, [12, 13, 17, 18]]))

    rng = np.random.default_rng(BOOTSTRAP_SEED)
    count_rows: list[dict] = []
    first_law_rows: list[dict] = []
    bin_rows: list[dict] = []
    symmetry_rows: list[dict] = []
    fit_bin_rows: list[dict] = []
    ift_rows: list[dict] = []

    for multiplier in range(1, 11):
        groups = EXPECTED_BLOCKS_PER_STREAM // multiplier
        usable = groups * multiplier
        grouped = shaped[:, :usable, :].reshape(
            EXPECTED_STREAMS, groups, multiplier, -1)
        entropy = grouped[:, :, :, 5].sum(axis=2)
        residual = grouped[:, :, :, 8].sum(axis=2)
        q_left = grouped[:, :, :, 2].sum(axis=2)
        q_right = grouped[:, :, :, 3].sum(axis=2)
        delta_energy = grouped[:, :, :, 4].sum(axis=2)
        recomputed_residual = q_left + q_right - delta_energy
        aggregation_identity_error = float(np.max(np.abs(
            residual - recomputed_residual)))

        negative = entropy < 0.0
        ci_low, ci_high = stream_bootstrap_interval(negative, rng)
        width, positive_counts, negative_counts = symmetric_counts(
            entropy.ravel())
        qualifying = (positive_counts >= 20) & (negative_counts >= 20)
        qualifying_pairs = int(np.sum(qualifying))
        negative_count = int(np.sum(negative))
        resolved = negative_count >= 1000 and qualifying_pairs >= 8
        count_rows.append({
            "time": 20 * multiplier,
            "multiplier": multiplier,
            "blocks_per_stream": groups,
            "n_blocks": entropy.size,
            "negative_count": negative_count,
            "negative_probability": negative.mean(),
            "negative_probability_ci_low": ci_low,
            "negative_probability_ci_high": ci_high,
            "fd_bin_width": width,
            "symmetric_pairs_count_ge_20_each": qualifying_pairs,
            "resolved_by_frozen_gate": int(resolved),
        })
        for index, (pos, neg) in enumerate(
                zip(positive_counts, negative_counts)):
            bin_rows.append({
                "time": 20 * multiplier,
                "bin_index": index,
                "abs_left": index * width,
                "abs_right": (index + 1) * width,
                "positive_count": int(pos),
                "negative_count": int(neg),
                "qualifies_ge_20_each": int(pos >= 20 and neg >= 20),
            })

        fit_indices = np.flatnonzero(qualifying)
        if fit_indices.size >= 2:
            centers = (fit_indices.astype(float) + 0.5) * width
            fit_positive = positive_counts[fit_indices].astype(float)
            fit_negative = negative_counts[fit_indices].astype(float)
            log_ratio = np.log(fit_positive / fit_negative)
            intercept, slope, r_squared = weighted_line_fit(
                centers, log_ratio, fit_positive, fit_negative)

            stream_positive, stream_negative = per_stream_symmetric_counts(
                entropy, width, len(positive_counts))
            draws = rng.integers(
                0, EXPECTED_STREAMS,
                size=(BOOTSTRAP_REPLICATES, EXPECTED_STREAMS))
            bootstrap_intercept = []
            bootstrap_slope = []
            for draw in draws:
                pos = stream_positive[draw].sum(axis=0)[fit_indices]
                neg = stream_negative[draw].sum(axis=0)[fit_indices]
                if np.any(pos == 0) or np.any(neg == 0):
                    continue
                boot_intercept, boot_slope, _ = weighted_line_fit(
                    centers, np.log(pos / neg), pos.astype(float),
                    neg.astype(float))
                bootstrap_intercept.append(boot_intercept)
                bootstrap_slope.append(boot_slope)
            if len(bootstrap_slope) < int(0.95 * BOOTSTRAP_REPLICATES):
                raise RuntimeError("too many invalid symmetry bootstrap draws")
            slope_ci = np.quantile(bootstrap_slope, [0.025, 0.975])
            intercept_ci = np.quantile(
                bootstrap_intercept, [0.025, 0.975])
            symmetry_rows.append({
                "time": 20 * multiplier,
                "n_blocks": entropy.size,
                "n_fit_bins": fit_indices.size,
                "intercept": intercept,
                "intercept_ci_low": intercept_ci[0],
                "intercept_ci_high": intercept_ci[1],
                "slope": slope,
                "slope_ci_low": slope_ci[0],
                "slope_ci_high": slope_ci[1],
                "weighted_r_squared": r_squared,
                "bootstrap_valid_replicates": len(bootstrap_slope),
            })
            for index, center, pos, neg, ratio in zip(
                    fit_indices, centers, fit_positive, fit_negative,
                    log_ratio):
                fit_bin_rows.append({
                    "time": 20 * multiplier,
                    "bin_index": int(index),
                    "a_center": center,
                    "positive_count": int(pos),
                    "negative_count": int(neg),
                    "log_count_ratio": ratio,
                    "log_count_ratio_se_poisson": math.sqrt(1.0 / pos +
                                                              1.0 / neg),
                })

        log_ift, weight_ess, maximum_weight_share, maximum, scaled = (
            log_mean_exponential_minus(entropy))
        per_stream_scaled_sum = scaled.sum(axis=1)
        draws = rng.integers(
            0, EXPECTED_STREAMS,
            size=(BOOTSTRAP_REPLICATES, EXPECTED_STREAMS))
        bootstrap_scaled_sum = per_stream_scaled_sum[draws].sum(axis=1)
        bootstrap_log_ift = maximum + np.log(
            bootstrap_scaled_sum /
            (EXPECTED_STREAMS * entropy.shape[1]))
        ift_ci = np.quantile(bootstrap_log_ift, [0.025, 0.975])
        ift_rows.append({
            "time": 20 * multiplier,
            "n_blocks": entropy.size,
            "log_mean_exp_minus_entropy_medium": log_ift,
            "ci_low": ift_ci[0],
            "ci_high": ift_ci[1],
            "exponential_weight_ess": weight_ess,
            "exponential_weight_ess_fraction": weight_ess / entropy.size,
            "maximum_single_weight_share": maximum_weight_share,
            "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        })

        flat_residual = residual.ravel()
        first_law_rows.append({
            "time": 20 * multiplier,
            "n_blocks": flat_residual.size,
            "mean": np.mean(flat_residual),
            "std": np.std(flat_residual, ddof=1),
            "rms": np.sqrt(np.mean(flat_residual ** 2)),
            "minimum": np.min(flat_residual),
            "q025": np.quantile(flat_residual, 0.025),
            "median": np.median(flat_residual),
            "q975": np.quantile(flat_residual, 0.975),
            "maximum": np.max(flat_residual),
            "max_abs": np.max(np.abs(flat_residual)),
            "aggregation_identity_max_abs_error": aggregation_identity_error,
        })

    write_csv(args.output / "negative_tail_counts.csv", count_rows,
              list(count_rows[0]))
    write_csv(args.output / "symmetric_bin_counts.csv", bin_rows,
              list(bin_rows[0]))
    write_csv(args.output / "first_law_residuals.csv", first_law_rows,
              list(first_law_rows[0]))
    if symmetry_rows:
        write_csv(args.output / "medium_entropy_symmetry.csv", symmetry_rows,
                  list(symmetry_rows[0]))
        write_csv(args.output / "medium_entropy_fit_bins.csv", fit_bin_rows,
                  list(fit_bin_rows[0]))
    write_csv(args.output / "medium_entropy_ift.csv", ift_rows,
              list(ift_rows[0]))

    with args.summary.open(newline="") as handle:
        summary_rows = list(csv.DictReader(handle))
    if len(summary_rows) != 1:
        raise RuntimeError("summary must contain exactly one row")
    summary = summary_rows[0]
    midpoint_failures = int(summary["midpoint_failure_count"])

    integrity_pass = (
        entropy_identity_error <= 2e-14 and
        balance_identity_error <= 2e-14 and
        linear_endpoint_error <= 2e-14 and
        angular_endpoint_error <= 2e-14 and
        action_minimum >= 0.0 and
        angle_minimum >= -math.pi and angle_maximum < math.pi and
        midpoint_failures == 0
    )
    audit = {
        "blocks": str(args.blocks.resolve()),
        "blocks_sha256": sha256(args.blocks),
        "summary": str(args.summary.resolve()),
        "summary_sha256": sha256(args.summary),
        "rows": int(data.shape[0]),
        "columns": int(data.shape[1]),
        "entropy_identity_max_abs_error": float(entropy_identity_error),
        "balance_identity_max_abs_error": float(balance_identity_error),
        "linear_endpoint_continuity_max_abs_error": float(linear_endpoint_error),
        "angular_endpoint_continuity_max_abs_error": float(angular_endpoint_error),
        "minimum_saved_action": action_minimum,
        "minimum_saved_angle": angle_minimum,
        "maximum_saved_angle": angle_maximum,
        "midpoint_failure_count": midpoint_failures,
        "integrity_pass": bool(integrity_pass),
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "feasibility_thresholds": {
            "minimum_negative_blocks": 1000,
            "minimum_symmetric_pairs": 8,
            "minimum_count_each_side_per_pair": 20,
        },
    }
    (args.output / "analysis_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n")

    t20 = count_rows[0]
    lines = [
        "# Direct-sampling feasibility verdict",
        "",
        f"Integrity gate: {'PASS' if integrity_pass else 'FAIL'}.",
        "",
        "The table below is copied from the raw-count analysis; no fit-window "
        "selection enters this decision.",
        "",
        "| t | blocks | negative Sigma_m | P(Sigma_m<0) | 95% stream-bootstrap CI | symmetric pairs >=20 each | resolved |",
        "|---:|---:|---:|---:|---:|---:|:---:|",
    ]
    for row in count_rows:
        lines.append(
            f"| {row['time']} | {row['n_blocks']} | {row['negative_count']} | "
            f"{row['negative_probability']:.9g} | "
            f"[{row['negative_probability_ci_low']:.9g},"
            f"{row['negative_probability_ci_high']:.9g}] | "
            f"{row['symmetric_pairs_count_ge_20_each']} | "
            f"{'yes' if row['resolved_by_frozen_gate'] else 'no'} |")
    lines.extend(["", (
        "At t=20, direct sampling resolves the two-sided medium-entropy tail "
        "under the frozen gate.  Cloning is not activated."
        if t20["resolved_by_frozen_gate"] else
        "At t=20, direct sampling does not resolve the two-sided medium-entropy "
        "tail under the frozen gate.  No cloning has been started; it may only "
        "be considered explicitly after reporting this failure."
    ), ""])
    (args.output / "FEASIBILITY_VERDICT.md").write_text("\n".join(lines))

    if not integrity_pass:
        raise RuntimeError("integrity gate failed; see analysis_audit.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
