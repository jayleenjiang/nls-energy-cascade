#!/usr/bin/env python3
"""Frozen stream-bootstrap equilibrium audit for compressed sampler output."""

from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np


BLOCK_TIME = 20.0
EXPECTED_ROWS = 1_000_064
EXPECTED_STREAMS = 128
EXPECTED_PER_STREAM = 7_813
BOOTSTRAPS = 5_000
BOOTSTRAP_SEED = 2_026_090_191
MIN_PAIR_COUNT = 200


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def weighted_fit(x: np.ndarray, y: np.ndarray, w: np.ndarray):
    design = np.column_stack([np.ones_like(x), x])
    normal = design.T @ (w[:, None] * design)
    beta = np.linalg.solve(normal, design.T @ (w * y))
    fitted = design @ beta
    ybar = np.sum(w * y) / np.sum(w)
    ss_res = np.sum(w * (y - fitted) ** 2)
    ss_tot = np.sum(w * (y - ybar) ** 2)
    r2 = float("nan") if ss_tot == 0.0 else 1.0 - ss_res / ss_tot
    return float(beta[0]), float(beta[1]), float(r2)


def stream_bootstrap_ci(stream_values: np.ndarray, rng: np.random.Generator):
    indices = rng.integers(0, stream_values.size,
                           size=(BOOTSTRAPS, stream_values.size))
    draws = stream_values[indices].mean(axis=1)
    return float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))


def read_case(archive: Path):
    entropy_rate = np.empty(EXPECTED_ROWS, dtype=np.float64)
    stream_ids = np.empty(EXPECTED_ROWS, dtype=np.int16)
    metric_names = ["q_left_rate", "q_right_rate", "entropy_rate",
                    "energy_current", "action_current"]
    sums = {name: np.zeros(EXPECTED_STREAMS, dtype=np.float64)
            for name in metric_names}
    counts = np.zeros(EXPECTED_STREAMS, dtype=np.int64)
    minimum_balance = math.inf
    maximum_balance = -math.inf
    sum_balance = 0.0
    sumsq_balance = 0.0

    process = subprocess.Popen(
        ["zstd", "-dc", str(archive)], stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, text=True, bufsize=1024 * 1024)
    assert process.stdout is not None
    reader = csv.DictReader(process.stdout)
    required = {"stream_id", "q_left", "q_right", "entropy_rate",
                "action_current", "energy_balance_error"}
    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        raise RuntimeError(f"missing columns in {archive}: {reader.fieldnames}")

    row_count = 0
    for row in reader:
        if row_count >= EXPECTED_ROWS:
            raise RuntimeError(f"too many rows in {archive}")
        stream = int(row["stream_id"])
        q_left = float(row["q_left"])
        q_right = float(row["q_right"])
        erate = float(row["entropy_rate"])
        action = float(row["action_current"])
        balance = float(row["energy_balance_error"])
        values = {
            "q_left_rate": q_left / BLOCK_TIME,
            "q_right_rate": q_right / BLOCK_TIME,
            "entropy_rate": erate,
            "energy_current": (q_left - q_right) / (2.0 * BLOCK_TIME),
            "action_current": action,
        }
        for name, value in values.items():
            sums[name][stream] += value
        counts[stream] += 1
        entropy_rate[row_count] = erate
        stream_ids[row_count] = stream
        minimum_balance = min(minimum_balance, balance)
        maximum_balance = max(maximum_balance, balance)
        sum_balance += balance
        sumsq_balance += balance * balance
        row_count += 1

    process.stdout.close()
    stderr = process.stderr.read() if process.stderr is not None else ""
    return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(f"zstd failed for {archive}: {stderr}")
    if row_count != EXPECTED_ROWS:
        raise RuntimeError(f"{archive}: {row_count} rows, expected {EXPECTED_ROWS}")
    if not np.all(counts == EXPECTED_PER_STREAM):
        raise RuntimeError(f"bad per-stream counts: {counts.tolist()}")
    if not np.all(np.isfinite(entropy_rate)):
        raise RuntimeError("non-finite entropy rate")

    stream_means = {name: values / counts for name, values in sums.items()}
    balance_audit = {
        "minimum": minimum_balance,
        "maximum": maximum_balance,
        "mean": sum_balance / row_count,
        "rms": math.sqrt(sumsq_balance / row_count),
    }
    return entropy_rate, stream_ids, stream_means, balance_audit


def analyze_case(label: str, temperature: float, archive: Path, output: Path):
    entropy_rate, stream_ids, stream_means, balance = read_case(archive)
    rng = np.random.default_rng(BOOTSTRAP_SEED + int(temperature))
    mean_rows = []
    for name, values in stream_means.items():
        estimate = float(np.mean(values))
        se = float(np.std(values, ddof=1) / math.sqrt(EXPECTED_STREAMS))
        ci_low, ci_high = stream_bootstrap_ci(values, rng)
        mean_rows.append({
            "label": label,
            "temperature": temperature,
            "observable": name,
            "estimate": estimate,
            "stream_se": se,
            "sigma_from_zero": estimate / se,
            "bootstrap_ci_low": ci_low,
            "bootstrap_ci_high": ci_high,
            "n_streams": EXPECTED_STREAMS,
            "n_blocks": EXPECTED_ROWS,
        })

    n_positive = int(np.count_nonzero(entropy_rate > 0.0))
    n_negative = int(np.count_nonzero(entropy_rate < 0.0))
    n_zero = int(np.count_nonzero(entropy_rate == 0.0))
    n_nonzero = n_positive + n_negative
    naive_z = (n_positive - n_negative) / math.sqrt(n_nonzero)
    signs = np.sign(entropy_rate)
    stream_sign_mean = np.bincount(stream_ids, weights=signs,
                                   minlength=EXPECTED_STREAMS) / EXPECTED_PER_STREAM
    sign_mean = float(np.mean(stream_sign_mean))
    sign_se = float(np.std(stream_sign_mean, ddof=1) /
                    math.sqrt(EXPECTED_STREAMS))
    sign_row = {
        "label": label,
        "temperature": temperature,
        "positive_count": n_positive,
        "negative_count": n_negative,
        "zero_count": n_zero,
        "positive_fraction_nonzero": n_positive / n_nonzero,
        "naive_counting_z": naive_z,
        "stream_sign_mean": sign_mean,
        "stream_sign_se": sign_se,
        "stream_sign_sigma_from_zero": sign_mean / sign_se,
    }

    absolute = np.abs(entropy_rate)
    q99 = float(np.quantile(absolute, 0.99))
    q25, q75 = np.quantile(entropy_rate, [0.25, 0.75])
    iqr = float(q75 - q25)
    fd_width = 2.0 * iqr / (EXPECTED_ROWS ** (1.0 / 3.0))
    if not (fd_width > 0.0 and q99 > 0.0):
        raise RuntimeError("degenerate entropy distribution")
    n_bins = int(np.clip(math.ceil(q99 / fd_width), 20, 80))
    edges = np.linspace(0.0, q99, n_bins + 1)
    positive_by_stream = np.zeros((EXPECTED_STREAMS, n_bins), dtype=np.int64)
    negative_by_stream = np.zeros((EXPECTED_STREAMS, n_bins), dtype=np.int64)

    positive_mask = entropy_rate > 0.0
    negative_mask = entropy_rate < 0.0
    positive_bins = np.searchsorted(edges, entropy_rate[positive_mask],
                                    side="right") - 1
    negative_bins = np.searchsorted(edges, -entropy_rate[negative_mask],
                                    side="right") - 1
    pos_valid = (positive_bins >= 0) & (positive_bins < n_bins)
    neg_valid = (negative_bins >= 0) & (negative_bins < n_bins)
    np.add.at(positive_by_stream,
              (stream_ids[positive_mask][pos_valid], positive_bins[pos_valid]), 1)
    np.add.at(negative_by_stream,
              (stream_ids[negative_mask][neg_valid], negative_bins[neg_valid]), 1)
    n_plus = positive_by_stream.sum(axis=0)
    n_minus = negative_by_stream.sum(axis=0)
    support = (n_plus >= MIN_PAIR_COUNT) & (n_minus >= MIN_PAIR_COUNT)
    if int(np.count_nonzero(support)) < 3:
        raise RuntimeError("fewer than three supported symmetric bin pairs")
    centers = 0.5 * (edges[:-1] + edges[1:])
    y = np.log(n_plus[support] / n_minus[support]) / BLOCK_TIME
    variance = (
        1.0 / n_plus[support] + 1.0 / n_minus[support]
    ) / (BLOCK_TIME ** 2)
    weights = 1.0 / variance
    intercept, slope, r2 = weighted_fit(centers[support], y, weights)

    bootstrap_slopes = np.empty(BOOTSTRAPS, dtype=np.float64)
    bootstrap_intercepts = np.empty(BOOTSTRAPS, dtype=np.float64)
    for bootstrap in range(BOOTSTRAPS):
        selected = rng.integers(0, EXPECTED_STREAMS, size=EXPECTED_STREAMS)
        bp = positive_by_stream[selected].sum(axis=0)[support]
        bm = negative_by_stream[selected].sum(axis=0)[support]
        by = np.log(bp / bm) / BLOCK_TIME
        bvar = (1.0 / bp + 1.0 / bm) / (BLOCK_TIME ** 2)
        bintercept, bslope, _ = weighted_fit(
            centers[support], by, 1.0 / bvar)
        bootstrap_intercepts[bootstrap] = bintercept
        bootstrap_slopes[bootstrap] = bslope

    slope_se = float(np.std(bootstrap_slopes, ddof=1))
    intercept_se = float(np.std(bootstrap_intercepts, ddof=1))
    fit_row = {
        "label": label,
        "temperature": temperature,
        "slope": slope,
        "slope_bootstrap_se": slope_se,
        "slope_sigma_from_zero": slope / slope_se,
        "slope_ci_low": float(np.quantile(bootstrap_slopes, 0.025)),
        "slope_ci_high": float(np.quantile(bootstrap_slopes, 0.975)),
        "intercept": intercept,
        "intercept_bootstrap_se": intercept_se,
        "intercept_ci_low": float(np.quantile(bootstrap_intercepts, 0.025)),
        "intercept_ci_high": float(np.quantile(bootstrap_intercepts, 0.975)),
        "weighted_r2": r2,
        "supported_bin_pairs": int(np.count_nonzero(support)),
        "total_symmetric_bins": n_bins,
        "abs_rate_q99": q99,
        "fd_width_unclamped": fd_width,
        "minimum_count_per_side": MIN_PAIR_COUNT,
    }
    bin_rows = []
    for index in range(n_bins):
        bin_rows.append({
            "label": label,
            "temperature": temperature,
            "bin_index": index,
            "abs_rate_left": edges[index],
            "abs_rate_right": edges[index + 1],
            "abs_rate_center": centers[index],
            "positive_count": int(n_plus[index]),
            "negative_count": int(n_minus[index]),
            "supported": int(support[index]),
            "log_ratio_over_t": (math.log(n_plus[index] / n_minus[index]) /
                                 BLOCK_TIME
                                 if n_plus[index] > 0 and n_minus[index] > 0
                                 else float("nan")),
        })

    return mean_rows, sign_row, fit_row, bin_rows, balance


def main() -> int:
    if len(sys.argv) != 3:
        print(f"usage: {sys.argv[0]} EXPERIMENT_DIR OUTPUT_DIR", file=sys.stderr)
        return 2
    experiment = Path(sys.argv[1]).resolve()
    output = Path(sys.argv[2]).resolve()
    output.mkdir(parents=True, exist_ok=True)
    cases = [("T6", 6.0), ("T10", 10.0)]
    all_means = []
    all_signs = []
    all_fits = []
    all_bins = []
    audits = {}
    for label, temperature in cases:
        archive = experiment / "raw" / f"{label}_blocks.csv.zst"
        means, signs, fit, bins, balance = analyze_case(
            label, temperature, archive, output)
        all_means.extend(means)
        all_signs.append(signs)
        all_fits.append(fit)
        all_bins.extend(bins)
        audits[label] = balance

    write_csv(output / "means_raw.csv", list(all_means[0]), all_means)
    write_csv(output / "sign_counts.csv", list(all_signs[0]), all_signs)
    write_csv(output / "symmetry_fit.csv", list(all_fits[0]), all_fits)
    write_csv(output / "symmetry_bins.csv", list(all_bins[0]), all_bins)
    with (output / "integrity_audit.json").open("w") as handle:
        json.dump({
            "expected_rows_per_case": EXPECTED_ROWS,
            "expected_streams": EXPECTED_STREAMS,
            "expected_blocks_per_stream": EXPECTED_PER_STREAM,
            "block_time": BLOCK_TIME,
            "bootstrap_replicates": BOOTSTRAPS,
            "bootstrap_seed_base": BOOTSTRAP_SEED,
            "balance_error": audits,
        }, handle, indent=2)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
