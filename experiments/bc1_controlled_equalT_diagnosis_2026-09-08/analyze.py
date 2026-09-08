#!/usr/bin/env python3
"""Frozen analysis for the projection-free BC1 equal-T diagnostic."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from statistics import NormalDist

import numpy as np


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
ANALYSIS = ROOT / "analysis"
SIMD_ROOT = ROOT.parent / "boundary_profiles_2026-09-05" / "curated_results" / "profiles"

CONFIG = {
    25: {
        "seeds": [2026090825, 2026090826],
        "simd_mean": 0.494670165647086,
        "simd_relative_range": 0.00893622658,
        "simd_max_abs_z": 8.678126289279529,
        "simd_pointwise_failures": 19,
    },
    100: {
        "seeds": [2026090900, 2026090901],
        "simd_mean": 0.249219499939971,
        "simd_relative_range": 0.0973137457,
        "simd_max_abs_z": 12.92852477395248,
        "simd_pointwise_failures": 71,
    },
}
OFFSETS = list(range(0, 256, 16))
EXPECTED_TRAJECTORIES = 512
POINTWISE_CRITICAL = 1.9599639845


def data_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(line for line in handle if not line.startswith("#")))


def read_simd(n: int) -> list[dict[str, str]]:
    return data_rows(SIMD_ROOT / f"BC1_T6_T6_n{n}_profile.csv")


def prefix(n: int, seed: int, offset: int) -> Path:
    return RAW / f"n{n}_seed{seed}_offset{offset:03d}"


def mean_se(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return values.mean(axis=0), values.std(axis=0, ddof=1) / math.sqrt(values.shape[0])


def analyze_n(n: int, require_complete: bool) -> dict:
    cfg = CONFIG[n]
    trajectory_records: dict[tuple[int, int], dict[int, tuple[float, float]]] = {}
    summaries: list[dict[str, str]] = []
    missing: list[str] = []
    integrity_errors: list[str] = []

    for seed in cfg["seeds"]:
        for offset in OFFSETS:
            p = prefix(n, seed, offset)
            trajectory_path = Path(str(p) + "_trajectory_profiles.csv")
            summary_path = Path(str(p) + "_summary.csv")
            metadata_path = Path(str(p) + "_metadata.csv")
            for needed in (trajectory_path, summary_path, metadata_path):
                if not needed.exists():
                    missing.append(str(needed))
            if not trajectory_path.exists() or not summary_path.exists():
                continue

            rows = data_rows(trajectory_path)
            if len(rows) != 16 * n:
                integrity_errors.append(f"{trajectory_path}: rows={len(rows)} expected={16*n}")
            for row in rows:
                stream = int(row["stream_id"])
                if not offset <= stream < offset + 16:
                    integrity_errors.append(f"{trajectory_path}: unexpected stream {stream}")
                key = (seed, stream)
                site = int(row["j"])
                action = float(row["mean_I"])
                sine = float(row["mean_sin_theta"])
                trajectory_records.setdefault(key, {})[site] = (action, sine)

            summary_rows = data_rows(summary_path)
            if len(summary_rows) != 1:
                integrity_errors.append(f"{summary_path}: summary rows={len(summary_rows)}")
            else:
                summaries.append(summary_rows[0])

    complete = not missing and not integrity_errors
    if require_complete and not complete:
        raise RuntimeError(
            f"n={n} incomplete: missing={len(missing)} integrity_errors={integrity_errors}"
        )
    if not trajectory_records:
        result = {"n": n, "complete": False, "missing": missing, "integrity_errors": integrity_errors}
        ANALYSIS.mkdir(parents=True, exist_ok=True)
        (ANALYSIS / f"n{n}_summary.json").write_text(json.dumps(result, indent=2) + "\n")
        return result

    keys = sorted(trajectory_records)
    action = np.full((len(keys), n), np.nan)
    sine = np.full((len(keys), n - 1), np.nan)
    for row_index, key in enumerate(keys):
        sites = trajectory_records[key]
        if len(sites) != n:
            integrity_errors.append(f"trajectory {key}: sites={len(sites)} expected={n}")
        for site, (a, s) in sites.items():
            action[row_index, site - 1] = a
            if site < n:
                sine[row_index, site - 1] = s

    if not np.isfinite(action).all() or not np.isfinite(sine).all():
        integrity_errors.append("non-finite trajectory profile value")

    action_mean, action_se = mean_se(action)
    sine_mean, sine_se = mean_se(sine)
    z = sine_mean / sine_se

    simd_rows = read_simd(n)
    simd_action_mean = np.array([float(row["mean_I"]) for row in simd_rows])
    simd_action_se = np.array([float(row["se_mean_I"]) for row in simd_rows])
    simd_sine_mean = np.array([float(row["mean_sin_theta"]) for row in simd_rows[:-1]])
    simd_sine_se = np.array([float(row["se_mean_sin_theta"]) for row in simd_rows[:-1]])
    simd_z = simd_sine_mean / simd_sine_se

    ANALYSIS.mkdir(parents=True, exist_ok=True)
    with (ANALYSIS / f"n{n}_action_comparison.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "j", "controlled_mean_I", "controlled_se_mean_I",
            "simd_mean_I", "simd_se_mean_I", "difference",
        ])
        for j in range(n):
            writer.writerow([
                j + 1, action_mean[j], action_se[j], simd_action_mean[j],
                simd_action_se[j], action_mean[j] - simd_action_mean[j],
            ])

    with (ANALYSIS / f"n{n}_bond_sine_comparison.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "j", "controlled_mean_sin", "controlled_se", "controlled_z",
            "simd_mean_sin", "simd_se", "simd_z", "controlled_over_simd_se",
        ])
        for j in range(n - 1):
            writer.writerow([
                j + 1, sine_mean[j], sine_se[j], z[j], simd_sine_mean[j],
                simd_sine_se[j], simd_z[j], sine_se[j] / simd_sine_se[j],
            ])

    controlled_spatial_mean = float(action_mean.mean())
    controlled_relative_range = float((action_mean.max() - action_mean.min()) / controlled_spatial_mean)
    mean_relative_difference = abs(controlled_spatial_mean - cfg["simd_mean"]) / cfg["simd_mean"]
    port_mean_pass = mean_relative_difference <= 0.05
    port_flatness_pass = controlled_relative_range <= cfg["simd_relative_range"]
    porting_pass = bool(port_mean_pass and port_flatness_pass and not integrity_errors and len(keys) == EXPECTED_TRAJECTORIES)

    bonferroni_critical = NormalDist().inv_cdf(1.0 - 0.05 / (2.0 * (n - 1)))
    max_index = int(np.argmax(np.abs(z)))
    max_abs_z = float(abs(z[max_index]))
    pointwise_failures = int(np.count_nonzero(np.abs(z) > POINTWISE_CRITICAL))
    bonferroni_pass = bool(max_abs_z <= bonferroni_critical)

    projection_count = sum(int(row["projection_count"]) for row in summaries)
    floor_count = sum(int(row["floor_count"]) for row in summaries)
    zero_radius_events = sum(int(row["zero_radius_events"]) for row in summaries)
    midpoint_failure_steps = sum(int(row["midpoint_failure_steps"]) for row in summaries)
    nonfinite_batches = sum(int(row["nonfinite"]) for row in summaries)
    min_sampled_action = min(float(row["min_sampled_action"]) for row in summaries)
    profile_samples = sorted({int(row["profile_samples_per_trajectory"]) for row in summaries})

    result = {
        "n": n,
        "complete": bool(complete and not integrity_errors and len(keys) == EXPECTED_TRAJECTORIES),
        "trajectory_count": len(keys),
        "expected_trajectory_count": EXPECTED_TRAJECTORIES,
        "base_seeds": cfg["seeds"],
        "profile_samples_per_trajectory": profile_samples,
        "missing": missing,
        "integrity_errors": integrity_errors,
        "porting": {
            "controlled_spatial_mean": controlled_spatial_mean,
            "simd_spatial_mean": cfg["simd_mean"],
            "mean_relative_difference": mean_relative_difference,
            "mean_within_5_percent": bool(port_mean_pass),
            "controlled_relative_range": controlled_relative_range,
            "simd_relative_range": cfg["simd_relative_range"],
            "relative_range_pass": bool(port_flatness_pass),
            "pass": porting_pass,
        },
        "sine_control": {
            "max_abs_z": max_abs_z,
            "max_abs_z_bond": max_index + 1,
            "pointwise_critical": POINTWISE_CRITICAL,
            "pointwise_failures": pointwise_failures,
            "bond_count": n - 1,
            "bonferroni_critical": bonferroni_critical,
            "bonferroni_pass": bonferroni_pass,
            "controlled_se_min": float(sine_se.min()),
            "controlled_se_mean": float(sine_se.mean()),
            "controlled_se_max": float(sine_se.max()),
            "simd_se_min": float(simd_sine_se.min()),
            "simd_se_mean": float(simd_sine_se.mean()),
            "simd_se_max": float(simd_sine_se.max()),
            "median_controlled_over_simd_se": float(np.median(sine_se / simd_sine_se)),
            "simd_max_abs_z": cfg["simd_max_abs_z"],
            "simd_pointwise_failures": cfg["simd_pointwise_failures"],
        },
        "numerics": {
            "projection_count": projection_count,
            "floor_count": floor_count,
            "zero_radius_events": zero_radius_events,
            "min_sampled_action": min_sampled_action,
            "midpoint_failure_steps": midpoint_failure_steps,
            "nonfinite_batches": nonfinite_batches,
        },
    }
    result["numerical_pass"] = bool(
        projection_count == 0 and floor_count == 0 and midpoint_failure_steps == 0
        and nonfinite_batches == 0 and result["complete"]
    )
    (ANALYSIS / f"n{n}_summary.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", choices=["25", "100", "all"], default="all")
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()
    ns = [25, 100] if args.n == "all" else [int(args.n)]
    results = [analyze_n(n, args.require_complete) for n in ns]
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
